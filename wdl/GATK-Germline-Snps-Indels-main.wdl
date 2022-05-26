version 1.0

import "subworkflows/processing-for-variant-discovery-gatk4.wdl" as mapper
import "subworkflows/haplotypecaller-gvcf-gatk4.wdl" as caller
import "subworkflows/tasks.wdl" as tasks

workflow GermlineSnpsIndelsGatk4 {
    input {
        String sampleName
        String refName
        Array[File] inFileFqs
        File refFa
        File refFai
        File refDict
        File? refAlt
        File refSa
        File refAnn
        File refBwt
        File refPac
        File refAmb
        
        File dbsnpVCF
        File dbsnpVCFTbi
        
        Array[File] knownIndelsSitesVCFs
        Array[File] knownIndelsSitesIdxs

        Boolean makeGVCF = true
        Boolean makeBamout = false
        
        String gatk_path
        String gotc_path
        String bwaCommandline = "bwa mem -K 100000000 -v 3 -t 16 -Y $bash_ref_fasta"
        Int executorMemory = 9

    }

    call mapper.PreProcessingForVariantDiscovery_GATK4 {
        input:
            sample_name = sampleName,
            ref_name    = refName,
            fastq_files = inFileFqs,

            ref_fasta       = refFa,
            ref_fasta_index = refFai,
            ref_dict        = refDict,
            ref_alt         = refAlt,
            ref_sa          = refSa,
            ref_ann         = refAnn,
            ref_bwt         = refBwt,
            ref_pac         = refPac,
            ref_amb         = refAmb,

            dbSNP_vcf       = dbsnpVCF,
            dbSNP_vcf_index = dbsnpVCFTbi,

            known_indels_sites_VCFs    = knownIndelsSitesVCFs,
            known_indels_sites_indices = knownIndelsSitesIdxs,

            bwa_commandline   = bwaCommandline,
            compression_level = 5,

            gatk_docker   = "us.gcr.io/broad-gatk/gatk:4.2.0.0",
            gatk_path     = gatk_path,
            gotc_docker   = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710",
            gotc_path     = gotc_path,
            python_docker = "python:2.7",

            flowcell_small_disk  = 100,
            flowcell_medium_disk = 200,
            agg_small_disk       = 200,
            agg_medium_disk      = 300,
            agg_large_disk       = 400,

            preemptible_tries    = 3,
            executorMemory  = executorMemory
    }

    call caller.HaplotypeCallerGvcf_GATK4 {
        input:
            input_bam       = PreProcessingForVariantDiscovery_GATK4.analysis_ready_bam,
            ref_dict        = refDict,
            ref_fasta       = refFa,
            ref_fasta_index = refFai,

            make_gvcf       = makeGVCF,
            make_bamout     = makeBamout,
            gatk_docker     = "us.gcr.io/broad-gatk/gatk:4.2.0.0",
            gatk_path       = gatk_path,
            gitc_docker     = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710",
            samtools_path   = "samtools",
            executorMemory  = executorMemory
    }

    call tasks.IndexBam {
        input:
            inFileBam = PreProcessingForVariantDiscovery_GATK4.analysis_ready_bam
    }

    call tasks.GenotypeGVCFs {
        input:
            inFileGVCF = HaplotypeCallerGvcf_GATK4.output_vcf,
            outBasename = basename(PreProcessingForVariantDiscovery_GATK4.analysis_ready_bam, ".bam"),
            refFa = refFa,
            refFai = refFai,
            refDict = refDict,
            dbsnpVCF = dbsnpVCF,
            dbsnpVCFTbi = dbsnpVCFTbi,
            executorMemory  = executorMemory
    }

    output {
        File outFileBam = IndexBam.outFileBam
        File outFileBai = IndexBam.outFileBai
        File outFileGVCF = GenotypeGVCFs.outFileGVCF
        File outFileVCF = GenotypeGVCFs.outFileVCF
        File outFileVCFIdx = GenotypeGVCFs.outFileVCFTbi
        File outFileBqsrReport = PreProcessingForVariantDiscovery_GATK4.bqsr_report
        File outFileDuplicationMetrics = PreProcessingForVariantDiscovery_GATK4.duplication_metrics
    }
}
