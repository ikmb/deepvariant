process vcf_compress_and_index {

	//publishDir "${params.outdir}/SVs", mode: 'copy'

	input:
	path(vcf)

	output:
	tuple path(vcf_gz),path(vcf_gz_tbi)

	script:
	
	vcf_gz = vcf.getBaseName() + ".vcf.gz"
	vcf_gz_tbi = vcf_gz + ".tbi"

	"""
		bgzip -c $vcf > $vcf_gz
		tabix $vcf_gz
	"""

}

process bed_compress_and_index {

	input:
	path(bed)

	output:
	tuple path(bed_gz),path(bed_gz_tbi)

	script:
	bed_gz = bed + ".gz"
	bed_gz_tbi = bed_gz + ".tbi"

	"""
		bgzip -c $bed > $bed_gz
		tabix $bed_gz
	"""

}
