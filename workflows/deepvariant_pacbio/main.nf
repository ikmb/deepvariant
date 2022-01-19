include { long_read_align } from "./../../modules/minimap/main.nf" params(params)
include { deepvariant } from "./../../modules/deepvariant/main.nf" params(params)
include { merge_bam_files } from "./../../modules/samtools/main.nf" params(params)
include { vcf_add_dbsnp ; vcf_index ; vcf_pass ; vcf_get_sample} from "./../../modules/vcf/main.nf" params(params)
include { pbsv_sig; pbsv_call } from "./../../modules/sv/main.nf" params(params)
include { vcf_compress_and_index } from "./../../modules/htslib/main.nf" params(params)

workflow DEEPVARIANT_PACBIO {

	take:
		reads
		bed
		repeats
		fastaGz
		gzFai
		gzi
		fai

	main:
		long_read_align(reads)
		merge_bam_files(long_read_align.out[0].groupTuple(by: [0,1]) )
		deepvariant(merge_bam_files.out[0],bed.collect(),fastaGz.collect(),gzFai.collect(),gzi.collect(),fai.collect() )
		vcf_index(deepvariant.out[1])
		vcf_pass(vcf_index.out)
		vcf_add_dbsnp(vcf_pass.out)
		pbsv_sig(long_read_align.out[0],repeats.collect())
		pbsv_call(pbsv_sig.out[1].collect())
		vcf_compress_and_index(pbsv_call.out)
		vcf_get_sample(deepvariant.out[2].combine(vcf_compress_and_index.out).collect(),"SVs")
	emit:
		bam = merge_bam_files.out[0]
		gvcf = deepvariant.out[0]
		vcf = vcf_add_dbsnp.out

}
