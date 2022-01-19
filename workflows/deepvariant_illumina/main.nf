include { trim } from "./../../modules/fastp/main.nf" params(params)
include { align } from "./../../modules/bwa/main.nf" params(params)
include { merge_and_dedup } from "./../../modules/samtools/main.nf" params(params)
include { deepvariant } from "./../../modules/deepvariant/main.nf" params(params)
include { vcf_index ; vcf_pass } from "./../../modules/vcf/main.nf" params(params)

workflow DEEPVARIANT_SHORT_READS {

	take:
		reads
		bed
		fastaGz
		gzFai
		gzi
		fai
		
	main:
		trim(reads)
		align(trim.out[0])
		merge_and_dedup(align.out[0])
		deepvariant(merge_and_dedup.out[0],bed.collect(),fastaGz.collect(),gzFai.collect(),gzi.collect(),fai.collect())
		vcf_index(deepvariant.out[1])
		vcf_pass(vcf_index.out)
	emit:
		gvcf = deepvariant.out[0]
		vcf = vcf_pass.out

}

