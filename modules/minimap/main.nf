process PBMM2 {

	tag "${meta.patient_id}|${meta.sample_id}"

	label 'pbmm2'

	input:
	tuple val(meta),path(fastq)

	output:
	tuple val(meta),path(bam),path(bai), emit: bam

	script:
	bam = fastq.getSimpleName() + ".bam"
	bai = bam + ".bai"
	rgid = fastq.getSimpleName()
	
	"""
		pbmm2 align ${params.mmi} $fastq $bam --preset CCS --sort --rg '@RG\\tID:${rgid}\\tSM:${meta.sample_id}\\tCN:CCGA' -j ${task.cpus}
	"""
	
}
