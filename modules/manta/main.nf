process manta {

        publishDir "${params.outdir}/${indivID}/${sampleID}/Manta", mode: 'copy'

        label 'manta'

        scratch true

        input:
        set val(indivID), val(sampleID), path(bam),path(bai)

        output:
        set val(indivID), val(sampleID), path(vcf)

        script:
        vcf = bam.getBaseName() + ".manta_diploidSV.vcf.gz"

        """
                configManta.py --bam $bam \
                        --referenceFasta $params.fasta \
                        --runDir manta \
                        --callRegions $params.bed \

                manta/runWorkflow.py -j ${task.cpus}

                cp manta/results/variants/diploidSV.vcf.gz $vcf
        """
}

