process DEEPVARIANT {

	tag "${meta.patient_id}|${meta.sample_id}"

        publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/DeepVariant", mode: 'copy'

        label 'deepvariant'

        scratch true

        input:
        tuple val(meta), path(bam),path(bai)
        path(bed)
	path(fastaGz)
	path(gzFai)
	path(gzi)
	path(fai)

        output:
        tuple val(meta),path(gvcf), emit: gvcf
        tuple val(meta),path(vcf), emit: vcf
	val(meta), emit: sample

        script:
        gvcf = bam.getBaseName() + ".g.vcf.gz"
        vcf = bam.getBaseName() + ".vcf.gz"

	def model = "WGS"
	def options = ""
	if (params.pacbio) {
		model = "PACBIO"
		if (params.intervals) {
			options = "--regions=${bed}"
		}
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

