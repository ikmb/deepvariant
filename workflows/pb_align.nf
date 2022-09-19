include { PBMM2 } from "./../modules/minimap/main.nf"
include { SAMTOOLS_MERGE_BAM; SAMTOOLS_INDEX } from "./../modules/samtools/main.nf"

workflow PB_ALIGN {

	take:
		reads

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

	emit:
		bam = SAMTOOLS_MERGE_BAM.out.bam
}
