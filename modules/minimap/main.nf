process PBMM2 {

    container 'quay.io/biocontainers/pbmm2:1.12.0--h9ee0642_0'

	tag "${meta.patient_id}|${meta.sample_id}"

	label 'long_parallel'

	input:
	tuple val(meta),path(fastq)
	val(mmi)

	output:
	tuple val(meta),path(bam),path(bai), emit: bam

	script:
	bam = fastq.getSimpleName() + ".bam"
	bai = bam + ".bai"
	rgid = fastq.getSimpleName()
	
	"""
	pbmm2 align ${mmi} $fastq $bam --preset CCS --sort --rg '@RG\\tID:${rgid}\\tSM:${meta.sample_id}\\tCN:CCGA' -j ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbmm2: \$( pbmm2 --version | head -n1 | sed -e "s/pbmm2 //g" )
    END_VERSIONS
	"""
	
}
