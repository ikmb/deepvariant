include { DEEPVARIANT } from "./../modules/deepvariant/main.nf"
include { VCF_ADD_DBSNP; VCF_INDEX ; VCF_PASS } from "./../modules/vcf/main.nf"
include { MANTA } from "./../modules/sv/main.nf"
include { WHATSHAP } from "./../modules/whatshap/main"

workflow DEEPVARIANT_SHORT_READS {

	take:
		bam
		bed
		fastaGz
		gzFai
		gzi
		fai
		
	main:
		
		ch_vcfs = Channel.from([])

		DEEPVARIANT(
			bam,
			bed.collect(),
			fastaGz.collect(),
			gzFai.collect(),
			gzi.collect(),
			fai.collect()
		)
		VCF_INDEX(DEEPVARIANT.out.vcf)
		VCF_PASS(VCF_INDEX.out.vcf)
		VCF_ADD_DBSNP(VCF_PASS.out.vcf)
		if (params.phase) {
			WHATSHAP(
				VCF_ADD_DBSNP.out.vcf.join(SAMTOOLS_MERGE_AND_DEDUP.out.bam)
			)
			ch_vcfs = WHATSHAP.out.vcf
		} else {
			ch_vcfs = VCF_ADD_DBSNP.out.vcf
		}

	emit:
		gvcf = DEEPVARIANT.out.gvcf
		vcf = ch_vcfs
		vcf_dv = VCF_ADD_DBSNP.out.vcf
}
