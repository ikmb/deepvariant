process deepvariant {

        publishDir "${params.outdir}/${indivID}/${sampleID}/DeepVariant", mode: 'copy'

        label 'deepvariant'

        scratch true

        input:
        tuple val(indivID), val(sampleID), path(bam),path(bai)
        path(bed)
	path(fastaGz)
	path(gzFai)
	path(gzi)
	path(fai)

        output:
        tuple val(indivID),val(sampleID),path(gvcf)
        tuple val(indivID),val(sampleID),path(vcf)
	val(sampleID)

        script:
        gvcf = bam.getBaseName() + ".g.vcf.gz"
        vcf = bam.getBaseName() + ".vcf.gz"

	def model = "WGS"
	if (params.pacbio) {
		model = "PACBIO"
	}
        """
		unset TMPDIR

                /opt/deepvariant/bin/run_deepvariant \
                --model_type=${model} \
                --ref=$fastaGz \
                --reads $bam \
                --output_vcf=$vcf \
                --output_gvcf=$gvcf \
                --regions=$bed \
                --num_shards=${task.cpus}
        """
}

