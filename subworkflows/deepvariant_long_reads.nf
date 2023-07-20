include { DEEPVARIANT } from "./../modules/deepvariant/main.nf"
include { TABIX } from "./../modules/htslib/tabix.nf"

ch_versions = Channel.from([])

workflow DEEPVARIANT_LONG_READS {

	take:
		bam
		bed
		ch_fasta

	main:

		DEEPVARIANT(
			bam,
			bed.collect(),
			ch_fasta 
		)

		ch_versions = ch_versions.mix(DEEPVARIANT.out.versions)

		TABIX(
			DEEPVARIANT.out.vcf
		)

		ch_versions = ch_versions.mix(TABIX.out.versions)
	
	emit:
		gvcf = DEEPVARIANT.out.gvcf
		vcf = TABIX.out.vcf
		versions = ch_versions

}
