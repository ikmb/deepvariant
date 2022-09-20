process BEDTOOLS_SLOP {

	label 'bedtools'

	input:
	path(bed)

	output:
	path(bed_padded), emit: bed

	script:
	bed_padded = bed.getBaseName() + ".padded.bed"

	"""
		bedtools slop -b 10 -g $params.fasta -i $bed > $bed_padded
	"""

}
