version 1.0

# WORKFLOW DEFINITION
workflow PreProcessingForVariantDiscovery_GATK4 {
  input {
    String sample_name
    String ref_name
    Array[File] fastq_files

    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File? ref_alt
    File ref_sa
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_amb

    File dbSNP_vcf
    File dbSNP_vcf_index

    Array[File] known_indels_sites_VCFs
    Array[File] known_indels_sites_indices

    String gatk_path
    String gotc_path

    String bwa_commandline = "bwa mem -K 100000000 -v 3 -t 16 -Y $bash_ref_fasta"
    Int compression_level = 5

    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.0.0"
    String gotc_docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
    String python_docker = "python:2.7"

    Int flowcell_small_disk = 100
    Int flowcell_medium_disk = 200
    Int agg_small_disk = 200
    Int agg_medium_disk = 300
    Int agg_large_disk = 400

    Int preemptible_tries = 3
    Int executorMemory
  }
    String base_file_name = sample_name + "." + ref_name

  # Get the version of BWA to include in the PG record in the header of the BAM produced
  # by MergeBamAlignment.
  call GetBwaVersion {
    input:
      docker_image = gotc_docker,
      bwa_path = gotc_path,
      preemptible_tries = preemptible_tries
  }

  call BwaMem {
      input:
        inFileFastqR1 = fastq_files[0],
        inFileFastqR2 = fastq_files[1],
        sampleName = sample_name,

        bwa_commandline = bwa_commandline,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        ref_alt = ref_alt,
        ref_sa = ref_sa,
        ref_ann = ref_ann,
        ref_bwt = ref_bwt,
        ref_pac = ref_pac,
        ref_amb = ref_amb,
        docker_image = gotc_docker,
        bwa_path = gotc_path,
        disk_size = flowcell_medium_disk,
        preemptible_tries = preemptible_tries,
        compression_level = compression_level
  }

  # Aggregate aligned+merged flowcell BAM files and mark duplicates
  # We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
  # to avoid having to spend time just merging BAM files.
  call MarkDuplicates {
    input:
      input_bam = BwaMem.output_bam,
      output_bam_basename = base_file_name + ".aligned.unsorted.duplicates_marked",
      metrics_filename = base_file_name + ".duplicate_metrics",
      docker_image = gatk_docker,
      gatk_path = gatk_path,
      disk_size = agg_large_disk,
      compression_level = compression_level,
      preemptible_tries = preemptible_tries,
      executorMemory = executorMemory
  }

  # Sort aggregated+deduped BAM file and fix tags
  call SortAndFixTags {
    input:
      input_bam = MarkDuplicates.output_bam,
      output_bam_basename = base_file_name + ".aligned.duplicate_marked.sorted",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      docker_image = gatk_docker,
      gatk_path = gatk_path,
      disk_size = agg_large_disk,
      preemptible_tries = 0,
      compression_level = compression_level,
      executorMemory = executorMemory,
      mem_size_gb = executorMemory - 2
  }

  # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel
  call BaseRecalibrator {
    input:
      input_bam = SortAndFixTags.output_bam,
      input_bam_index = SortAndFixTags.output_bam_index,
      recalibration_report_filename = base_file_name + ".recal_data.csv",
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index,
      known_indels_sites_VCFs = known_indels_sites_VCFs,
      known_indels_sites_indices = known_indels_sites_indices,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      docker_image = gatk_docker,
      gatk_path = gatk_path,
      disk_size = agg_small_disk,
      preemptible_tries = preemptible_tries,
      executorMemory = executorMemory
  }

#  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping_with_unmapped) {

# Apply the recalibration model by interval
  call ApplyBQSR {
      input:
        input_bam = SortAndFixTags.output_bam,
        input_bam_index = SortAndFixTags.output_bam_index,
        output_bam_basename = base_file_name + ".aligned.duplicates_marked.recalibrated",
        recalibration_report = BaseRecalibrator.recalibration_report,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        docker_image = gatk_docker,
        gatk_path = gatk_path,
        disk_size = agg_small_disk,
        preemptible_tries = preemptible_tries,
        executorMemory = executorMemory
  }

  # Outputs that will be retained when execution is complete
  output {
    File duplication_metrics = MarkDuplicates.duplicate_metrics
    File bqsr_report = BaseRecalibrator.recalibration_report
    File analysis_ready_bam = ApplyBQSR.recalibrated_bam
  }
}

# TASK DEFINITIONS

# Get version of BWA
task GetBwaVersion {
  input {
    Float mem_size_gb = 1
    Int preemptible_tries
    String docker_image
    String bwa_path
  }

  command {
    # Not setting "set -o pipefail" here because /bwa has a rc=1 and we don't want to allow rc=1 to succeed
    # because the sed may also fail with that error and that is something we actually want to fail on.
    ~{bwa_path}bwa 2>&1 | \
    grep -e '^Version' | \
    sed 's/Version: //'
  }
  runtime {
    preemptible: preemptible_tries
    docker: docker_image
    memory: "~{mem_size_gb} GiB"
  }
  output {
    String version = "BWA-0.7.17"
  }
}

# BWA MEM for alignment
task BwaMem {
  # This is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit),
  # listing the reference contigs that are "alternative". Leave blank in JSON for legacy
  # references such as b37 and hg19.
  input {
    File inFileFastqR1
    File inFileFastqR2
    String sampleName

    String bwa_commandline
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File? ref_alt
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa

    Float mem_size_gb = 14
    String num_cpu = 16

    Int compression_level
    Int preemptible_tries
    Int disk_size

    String docker_image
    String bwa_path
  }
  Int command_mem_gb = ceil(mem_size_gb/2)

  command {
    set -o pipefail
    set -e

    # set the bash variable needed for the command-line
    bash_ref_fasta=~{ref_fasta}

    # bwa reference preload.
    # ~{bwa_path}/bwa shm ~{ref_fasta}

    ~{bwa_path}~{bwa_commandline} \
      -R "@RG\tID:~{sampleName}\tLB:~{sampleName}\tSM:~{sampleName}\tPL:ILLUMINA" \
      ~{inFileFastqR1} ~{inFileFastqR2} \
    | \
    samtools view -1 - > ~{sampleName}.bam
  }
  runtime {
    preemptible: preemptible_tries
    docker: docker_image
    memory: "~{mem_size_gb} GiB"
    cpu: num_cpu
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bam = "~{sampleName}.bam"
  }
}

# Sort BAM file by coordinate order and fix tag values for NM and UQ
task SortAndFixTags {
  input {
    File input_bam
    String output_bam_basename
    File ref_dict
    File ref_fasta
    File ref_fasta_index

    Int compression_level
    Int preemptible_tries
    Int disk_size
    Float mem_size_gb = 10

    String docker_image
    String gatk_path
    Int executorMemory
  }
    Int command_mem_gb_sort = ceil(mem_size_gb) - 1
    Int command_mem_gb_fix = ceil((mem_size_gb - 1)/10)
    Int maxMem = executorMemory - 2

  command {
    set -o pipefail

    ~{gatk_path} --java-options "-Dsamjdk.compression_level=~{compression_level} -Xmx~{maxMem}G -Xms~{command_mem_gb_sort}G" \
      SortSam \
      --INPUT ~{input_bam} \
      --OUTPUT /dev/stdout \
      --SORT_ORDER "coordinate" \
      --CREATE_INDEX false \
      --CREATE_MD5_FILE false \
    | \
    ~{gatk_path} --java-options "-Dsamjdk.compression_level=~{compression_level} -Xmx~{maxMem}G -Xms~{command_mem_gb_fix}G" \
      SetNmMdAndUqTags \
      --INPUT /dev/stdin \
      --OUTPUT ~{output_bam_basename}.bam \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true \
      --REFERENCE_SEQUENCE ~{ref_fasta}
  }
  runtime {
    preemptible: preemptible_tries
    docker: docker_image
    memory: "~{mem_size_gb} GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File output_bam_index = "~{output_bam_basename}.bai"
    File output_bam_md5 = "~{output_bam_basename}.bam.md5"
  }
}

# Mark duplicate reads to avoid counting non-independent observations
task MarkDuplicates {
  input {
    File input_bam
    String output_bam_basename
    String metrics_filename

    Int compression_level
    Int preemptible_tries
    Int disk_size
    Float mem_size_gb = 7.5

    String docker_image
    String gatk_path
    Int executorMemory
  }
    Int command_mem_gb = ceil(mem_size_gb) - 2
    Int maxMem = executorMemory - 2
 # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly.
 # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
 # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
  command {
    ~{gatk_path} --java-options "-Dsamjdk.compression_level=~{compression_level} -Xmx~{maxMem}G -Xms~{command_mem_gb}G" \
      MarkDuplicates \
      --INPUT ~{input_bam} \
      --OUTPUT ~{output_bam_basename}.bam \
      --METRICS_FILE ~{metrics_filename} \
      --VALIDATION_STRINGENCY SILENT \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
      --ASSUME_SORT_ORDER "queryname" \
      --CREATE_MD5_FILE true \
      --READ_NAME_REGEX null
  }
  runtime {
    preemptible: preemptible_tries
    docker: docker_image
    memory: "~{mem_size_gb}  GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File duplicate_metrics = "~{metrics_filename}"
  }
}

# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
  input {
    File input_bam
    File input_bam_index
    String recalibration_report_filename
    File dbSNP_vcf
    File dbSNP_vcf_index
    Array[File] known_indels_sites_VCFs
    Array[File] known_indels_sites_indices
    File ref_dict
    File ref_fasta
    File ref_fasta_index

    Int preemptible_tries
    Int disk_size
    Float mem_size_gb = 6

    String docker_image
    String gatk_path
    Int executorMemory
  }
  Int command_mem_gb = ceil(mem_size_gb) - 2
  Int maxMem = executorMemory - 2

  command {
    ~{gatk_path} --java-options "-Xmx~{maxMem}G -Xms~{command_mem_gb}G" \
      BaseRecalibrator \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      --use-original-qualities \
      -O ~{recalibration_report_filename} \
      --known-sites ~{dbSNP_vcf} \
      --known-sites ~{sep=" --known-sites " known_indels_sites_VCFs}
  }
  runtime {
    preemptible: preemptible_tries
    docker: docker_image
    memory: "~{mem_size_gb} GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File recalibration_report = "~{recalibration_report_filename}"
  }
}

# Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {
  input {
    File input_bam
    File input_bam_index
    String output_bam_basename
    File recalibration_report
    File ref_dict
    File ref_fasta
    File ref_fasta_index

    Int preemptible_tries
    Int disk_size
    Float mem_size_gb = 4

    String docker_image
    String gatk_path
    Int executorMemory
  }
  Int command_mem_gb = ceil(mem_size_gb) - 1
  Int maxMem = executorMemory - 2

  command {
    ~{gatk_path} --java-options "-Xmx~{maxMem}G -Xms~{command_mem_gb}G" \
      ApplyBQSR \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      -O ~{output_bam_basename}.bam \
      -bqsr ~{recalibration_report} \
      --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
      --add-output-sam-program-record \
      --create-output-bam-md5 \
      --use-original-qualities
  }
  runtime {
    preemptible: preemptible_tries
    docker: docker_image
    memory: "~{mem_size_gb} GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File recalibrated_bam = "~{output_bam_basename}.bam"
  }
}
