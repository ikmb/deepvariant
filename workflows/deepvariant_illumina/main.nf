include { FASTP } from "./../../modules/fastp/main.nf"
include { BWA } from "./../../modules/bwa/main.nf"
include { SAMTOOLS_MERGE_AND_DEDUP} from "./../../modules/samtools/main.nf"
include { DEEPVARIANT } from "./../../modules/deepvariant/main.nf"
include { VCF_ADD_DBSNP; VCF_INDEX ; VCF_PASS } from "./../../modules/vcf/main.nf"
include { MANTA } from "./../../modules/sv/main.nf"
include { BED_COMPRESS_AND_INDEX } from "./../../modules/htslib/main.nf"

workflow DEEPVARIANT_SHORT_READS {

	take:
		reads
		bed
		fastaGz
		gzFai
		gzi
		fai
		
	main:
		FASTP(reads)
		BWA(FASTP.out.reads)
		SAMTOOLS_MERGE_AND_DEDUP(
			BWA.out.bam.groupTuple(by: [0,1])
		)
		DEEPVARIANT(
			SAMTOOLS_MERGE_AND_DEDUP.out.bam,
			bed.collect(),
			fastaGz.collect(),
			gzFai.collect(),
			gzi.collect(),
			fai.collect()
		)
		VCF_INDEX(DEEPVARIANT.out.vcf)
		VCF_PASS(VCF_INDEX.out.vcf)
		VCF_ADD_DBSNP(VCF_PASS.out.vcf)
		BED_COMPRESS_AND_INDEX(bed)
		MANTA(
			SAMTOOLS_MERGE_AND_DEDUP.out.bam,
			BED_COMPRESS_AND_INDEX.out.bed.collect()
		)
	emit:
		gvcf = DEEPVARIANT.out.gvcf
		vcf = VCF_ADD_DBSNP.out.vcf.mix(MANTA.out[0],MANTA.out[1],MANTA.out[2])
		bam = SAMTOOLS_MERGE_AND_DEDUP.out.bam
		vcf_dv = VCF_ADD_DBSNP.out.vcf
}

