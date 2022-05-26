version 1.0

# WORKFLOW DEFINITION 
workflow HaplotypeCallerGvcf_GATK4 {
  input {
    File input_bam
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    String gatk_path

    Boolean make_gvcf = false
    Boolean make_bamout = false
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.0.0"
    String gitc_docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
    String samtools_path = "samtools"
    Int executorMemory
  }  

    #is the input a cram file?
    Boolean is_cram = sub(basename(input_bam), ".*\\.", "") == "cram"

    String sample_basename = if is_cram then  basename(input_bam, ".cram") else basename(input_bam, ".bam")
    String vcf_basename = sample_basename
    String output_suffix = if make_gvcf then ".g.vcf.gz" else ".vcf.gz"
    String output_filename = vcf_basename + output_suffix

    # Generate GVCF by interval
    call HaplotypeCaller {
      input:
        input_bam = input_bam,
        samtools_path = samtools_path,
        output_filename = output_filename,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        make_gvcf = make_gvcf,
        make_bamout = make_bamout,
        docker = gatk_docker,
        gatk_path = gatk_path,
        executorMemory = executorMemory
    }

  # Outputs that will be retained when execution is complete
  output {
    File output_vcf = HaplotypeCaller.output_vcf
    File output_vcf_index = HaplotypeCaller.output_vcf_index
  }
}

# TASK DEFINITIONS

task CramToBamTask {
  input {
    # Command parameters
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File input_cram
    String sample_name

    # Runtime parameters
    String docker
    Int? machine_mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? preemptible_attempts
    String samtools_path
  }
    Float output_bam_size = size(input_cram, "GB") / 0.60
    Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
    Int disk_size = ceil(size(input_cram, "GB") + output_bam_size + ref_size) + 20
  
  command {
    set -e
    set -o pipefail

    ~{samtools_path} view -h -T ~{ref_fasta} ~{input_cram} |
    ~{samtools_path} view -b -o ~{sample_name}.bam -
    ~{samtools_path} index -b ~{sample_name}.bam
    mv ~{sample_name}.bam.bai ~{sample_name}.bai
  }
  runtime {
    docker: docker
    memory: select_first([machine_mem_gb, 15]) + " GB"
    disks: "local-disk " + select_first([disk_space_gb, disk_size]) + if use_ssd then " SSD" else " HDD"
    preemptible: select_first([preemptible_attempts, 3])
 }
  output {
    File output_bam = "~{sample_name}.bam"
    File output_bai = "~{sample_name}.bai"
  }
}

# HaplotypeCaller per-sample in GVCF mode
task HaplotypeCaller {
  input {
    # Command parameters
    File input_bam
    String output_filename
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Float? contamination
    Boolean make_gvcf
    Boolean make_bamout
    String samtools_path

    String? gcs_project_for_requester_pays

    String gatk_path
    String? java_options

    # Runtime parameters
    String docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? preemptible_attempts
    Int executorMemory
  }

  String java_opt = select_first([java_options, "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"]) 

  Int machine_mem_gb = select_first([mem_gb, 7])
  Int command_mem_gb = machine_mem_gb - 1
  Int maxMem = executorMemory - 2

  # Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")

  String vcf_basename = if make_gvcf then  basename(output_filename, ".gvcf") else basename(output_filename, ".vcf")
  String bamout_arg = if make_bamout then "-bamout ~{vcf_basename}.bamout.bam" else ""

  parameter_meta {
    input_bam: {
      description: "a bam file",
      localization_optional: true
    }
  }
  command {
    set -e

    ~{samtools_path} index -b ~{input_bam}

    ~{gatk_path} --java-options "-Xmx~{maxMem}G -Xms~{command_mem_gb}G ~{java_opt}" \
      HaplotypeCaller \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      -O ~{output_filename} \
      -contamination ~{default="0" contamination} \
      -G StandardAnnotation -G StandardHCAnnotation ~{true="-G AS_StandardAnnotation" false="" make_gvcf} \
      -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
      ~{true="-ERC GVCF" false="" make_gvcf} \
      ~{bamout_arg}

    # Cromwell doesn't like optional task outputs, so we have to touch this file.
    touch ~{vcf_basename}.bamout.bam 
  }
  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    preemptible: select_first([preemptible_attempts, 3])
  }
  output {
    File output_vcf = "~{output_filename}"
    File output_vcf_index = "~{output_filename}.tbi"
    File bamout = "~{vcf_basename}.bamout.bam"
  }
}
# Merge GVCFs generated per-interval for the same sample
task MergeGVCFs {
  input {
    # Command parameters
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String output_filename

    String gatk_path

    # Runtime parameters
    String docker
    Int? mem_gb
    Int? disk_space_gb
    Int? preemptible_attempts
  }
    Boolean use_ssd = false
    Int machine_mem_gb = select_first([mem_gb, 3])
    Int command_mem_gb = machine_mem_gb - 1
  
  command {
  set -e

    ~{gatk_path} --java-options "-Xmx~{command_mem_gb}G -Xms~{command_mem_gb}G"  \
      MergeVcfs \
      --INPUT ~{sep=' --INPUT ' input_vcfs} \
      --OUTPUT ~{output_filename}
  }
  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + select_first([disk_space_gb, 100]) + if use_ssd then " SSD" else " HDD"
    preemptible: select_first([preemptible_attempts, 3])
  }
  output {
    File output_vcf = "~{output_filename}"
    File output_vcf_index = "~{output_filename}.tbi"
  }
}

