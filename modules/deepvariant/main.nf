process DEEPVARIANT {

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
        tuple val(indivID),val(sampleID),path(gvcf), emit: gvcf
        tuple val(indivID),val(sampleID),path(vcf), emit: vcf
	tuple val(indivID),val(sampleID), emit: sample

        script:
        gvcf = bam.getBaseName() + ".g.vcf.gz"
        vcf = bam.getBaseName() + ".vcf.gz"

	def model = "WGS"
	def options = ""
	if (params.pacbio) {
		model = "PACBIO"
	} else {
		options = "--regions=${bed}"
	}
        """
		set TMPDIR=\$PWD

                /opt/deepvariant/bin/run_deepvariant \
                --model_type=${model} \
                --ref=$fastaGz \
                --reads $bam \
                --output_vcf=$vcf \
                --output_gvcf=$gvcf \
                $options \
                --num_shards=${task.cpus} \
		--intermediate_results_dir tmp_data

        """
}

