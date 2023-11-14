include { DEEPVARIANT } from "./../modules/deepvariant/main.nf"
include { TABIX } 		from "./../modules/htslib/tabix.nf"
include { GLNEXUS } 	from "./../modules/glnexus/main"

ch_versions = Channel.from([])
ch_vcfs     = Channel.from([])

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
		ch_vcfs = ch_vcfs.mix(TABIX.out.vcf)

		if (params.joint_calling) {
            GLNEXUS(
                DEEPVARIANT.out.gvcf.map { m,g -> g }.collect(),
                bed
			)

			ch_versions = ch_versions.mix(GLNEXUS.out.versions)
			
		}
	
	emit:
		gvcf = DEEPVARIANT.out.gvcf
		vcf = ch_vcfs
		versions = ch_versions

}
