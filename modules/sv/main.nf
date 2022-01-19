process pbsv_sig {

	label 'pbsv'

	input:
	tuple val(indivID),val(sampleID),path(bam),path(pbi)
	path(repeat_ref)

	output:
	tuple val(indivID),val(sampleID),path(sig)
	path(sig)

	script:
	sig = bam.getBaseName() + ".svsig.gz"

	"""
		pbsv discover --ccs --tandem-repeats $repeat_ref $bam $sig
	"""
}

process pbsv_call {

	//publishDir "${params.outdir}/SVs", mode: 'copy'

	label 'pbsv'

	input:
	path(sigs)

	output:
	path(vcf)

	script:
	vcf = "SVs.vcf"

	"""
		pbsv call -j ${task.cpus} ${params.fasta} ${sigs} $vcf
	"""

}
