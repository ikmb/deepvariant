include { DEEPVARIANT } from "./../modules/deepvariant/main.nf"
include { SAMTOOLS_MERGE_BAM; SAMTOOLS_INDEX } from "./../modules/samtools/main.nf"
include { VCF_ADD_DBSNP ; VCF_INDEX ; VCF_PASS ; VCF_GET_SAMPLE } from "./../modules/vcf/main.nf"
include { VCF_COMPRESS_AND_INDEX } from "./../modules/htslib/main.nf" 

workflow DEEPVARIANT_LONG_READS {

	take:
		bam
		bed
		fastaGz
		gzFai
		gzi
		fai

	main:

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
	emit:
		gvcf = DEEPVARIANT.out.gvcf
		vcf = VCF_ADD_DBSNP.out.mix(VCF_GET_SAMPLE.out.vcf)
		vcf_dv = VCF_ADD_DBSNP.out.vcf
}
