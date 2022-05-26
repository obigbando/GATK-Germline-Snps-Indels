version 1.0

task GenotypeGVCFs {
    input {
        File inFileGVCF
        String outBasename
        File refFa
        File refFai
        File refDict
        File dbsnpVCF
        File dbsnpVCFTbi
        Int executorMemory
        String dockerImage = "atgenomix.azurecr.io/atgenomix/seqslab_runtime-1.4_ubuntu-18.04_gatk4-4.2.0.0:latest"
    }
    
    Int maxMem = executorMemory - 2

    command <<<
        set -e -o pipefail
        
        tabix ~{inFileGVCF}

        gatk --java-options "-Xmx~{maxMem}g -Xms~{maxMem}g" \
        GenotypeGVCFs \
        --variant ~{inFileGVCF} \
        --output ~{outBasename}.vcf.gz \
        --reference ~{refFa} \
        --dbsnp ~{dbsnpVCF} \
        --annotation-group StandardAnnotation --annotation-group AS_StandardAnnotation \

        mv ~{inFileGVCF} ~{outBasename}.g.vcf.gz
    >>>

    runtime {
        docker: dockerImage
        cpu: "1"
        memory: "~{maxMem} GB"
    }
    
    output {
        File outFileGVCF = "~{outBasename}.g.vcf.gz"
        File outFileVCF = "~{outBasename}.vcf.gz"
        File outFileVCFTbi = "~{outBasename}.vcf.gz.tbi"
    }
}

task IndexBam {
    input {
        File inFileBam
        String dockerImage = "atgenomix.azurecr.io/atgenomix/seqslab_runtime-1.4_ubuntu-18.04_gatk4-4.2.0.0:latest"
    }

    String outFileBamPrefix = basename(inFileBam, ".bam")

    command <<<
        set -e -o pipefail
        
        samtools index -@ 16 -b ~{inFileBam} ~{outFileBamPrefix}.bai
        mv ~{inFileBam} ~{outFileBamPrefix}.bam
    >>>

    runtime {
        docker: dockerImage
        cpu: "1"
        memory: "7 GB"
    }

    output {
        File outFileBam = "~{outFileBamPrefix}.bam"
        File outFileBai = "~{outFileBamPrefix}.bai"
    }
}

task IndexVcf {
    input {
        File inFileVCF
        String dockerImage = "atgenomix.azurecr.io/atgenomix/seqslab_runtime-1.4_ubuntu-18.04_gatk4-4.2.0.0:latest"
    }

    command <<<
        set -e -o pipefail
        
        bcftools index -t ~{inFileVCF} -o ~{inFileVCF}.tbi
    >>>

    runtime {
        docker: dockerImage
        cpu: "1"
        memory: "7 GB"
    }

    output {
        File outFileVCFTbi = "~{inFileVCF}.tbi"
    }
}


