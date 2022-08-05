include { PBMM2 } from "./../../modules/minimap/main.nf"
include { DEEPVARIANT } from "./../../modules/deepvariant/main.nf"
include { SAMTOOLS_MERGE_BAM; SAMTOOLS_INDEX } from "./../../modules/samtools/main.nf"
include { VCF_ADD_DBSNP ; VCF_INDEX ; VCF_PASS ; VCF_GET_SAMPLE } from "./../../modules/vcf/main.nf"
include { PBSV_SIG; PBSV_CALL } from "./../../modules/sv/main.nf"
include { VCF_COMPRESS_AND_INDEX } from "./../../modules/htslib/main.nf" 

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
		PBMM2(reads)

		bam_mapped = PBMM2.out.bam.map { meta, bam, bai ->
                        new_meta = [:]
                        new_meta.patient_id = meta.patient_id
                        new_meta.sample_id = meta.sample_id
                        def groupKey = meta.sample_id
                        tuple( groupKey, new_meta, bam, bai)
                }.groupTuple(by: [0,1]).map { g ,new_meta ,bam, bai -> [ new_meta, bam, bai ] }

                bam_mapped.branch {
                        single:   it[1].size() == 1
                        multiple: it[1].size() > 1
                }.set { bam_to_merge }

		// Merge and index multi-cell samples
		SAMTOOLS_MERGE_BAM(
			bam_to_merge.multiple.map {  m,b,i -> tuple(m,b) }
		)

		ch_bams = bam_to_merge.single.mix(SAMTOOLS_MERGE_BAM.out.bam)

		DEEPVARIANT(
			ch_bams,
			bed.collect(),
			fastaGz.collect(),
			gzFai.collect(),
			gzi.collect(),
			fai.collect() 
		)
		VCF_INDEX(DEEPVARIANT.out.vcf)
		VCF_PASS(VCF_INDEX.out.vcf)
		VCF_ADD_DBSNP(VCF_PASS.out.vcf)
		PBSV_SIG(
			PBMM2.out.bam,
			repeats.collect()
		)
		PBSV_CALL(PBSV_SIG.out[1].collect())
		VCF_COMPRESS_AND_INDEX(PBSV_CALL.out.vcf)
		VCF_GET_SAMPLE(
			DEEPVARIANT.out.sample.combine(VCF_COMPRESS_AND_INDEX.out).collect(),
			"SVs"
		)
	emit:
		bam = SAMTOOLS_MERGE_BAM.out.bam
		gvcf = DEEPVARIANT.out.gvcf
		vcf = VCF_ADD_DBSNP.out.mix(VCF_GET_SAMPLE.out.vcf)
		vcf_dv = VCF_ADD_DBSNP.out.vcf
}
