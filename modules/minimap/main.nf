process long_read_align {

	label 'pbmm2'

	input:
	tuple val(indivID),val(sampleID),path(fastq)

	output:
	tuple val(indivID),val(sampleID),path(bam),path(bai)

	script:
	bam = fastq.getSimpleName() + ".bam"
	bai = bam + ".bai"
	rgid = fastq.getSimpleName()
	
	"""
		pbmm2 align ${params.mmi} $fastq $bam --preset CCS --sort --rg '@RG\\tID:${rgid}\\tSM:${sampleID}\\tCN:${params.center}' -j ${task.cpus}
	"""
	
}