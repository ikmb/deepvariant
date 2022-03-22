include { trim } from "./../../modules/fastp/main.nf" params(params)
include { align } from "./../../modules/bwa/main.nf" params(params)
include { merge_and_dedup } from "./../../modules/samtools/main.nf" params(params)
include { deepvariant } from "./../../modules/deepvariant/main.nf" params(params)
include { vcf_add_dbsnp; vcf_index ; vcf_pass } from "./../../modules/vcf/main.nf" params(params)
include { manta } from "./../../modules/sv/main.nf" params(params)
include { bed_compress_and_index } from "./../../modules/htslib/main.nf" params(params)

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
		merge_and_dedup(
			align.out.bam.groupTuple(by: [0,1])
		)
		deepvariant(merge_and_dedup.out.bam,bed.collect(),fastaGz.collect(),gzFai.collect(),gzi.collect(),fai.collect())
		vcf_index(deepvariant.out[1])
		vcf_pass(vcf_index.out)
		vcf_add_dbsnp(vcf_pass.out)
		bed_compress_and_index(bed)
		manta(merge_and_dedup.out[0],bed_compress_and_index.out.collect())
	emit:
		gvcf = deepvariant.out[0]
		vcf = vcf_add_dbsnp.out.mix(manta.out[0],manta.out[1],manta.out[2])
		bam = merge_and_dedup.out.bam
		vcf_dv = vcf_add_dbsnp.out.vcf
}

