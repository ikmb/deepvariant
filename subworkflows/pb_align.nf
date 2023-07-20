include { PBMM2 }				from "./../modules/minimap/main"
include { SAMTOOLS_MERGE }		from "./../modules/samtools/merge/main"
include { SAMTOOLS_INDEX }		from "./../modules/samtools/index/main"

ch_bams 	= Channel.from([])
ch_versions = Channel.from([])

workflow PB_ALIGN {

	take:
		reads

	main:
		PBMM2(reads)

		ch_versions = ch_versions.mix(PBMM2.out.versions)

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
		SAMTOOLS_MERGE(
			bam_to_merge.multiple.map {  m,b,i -> tuple(m,b) }
		)

		ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

        SAMTOOLS_INDEX(
			SAMTOOLS_MERGE.out.bam
		)

		ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

		ch_bams = SAMTOOLS_INDEX.out.bam.mix(bam_to_merge.single)

	emit:
		bam = ch_bams
		versions = ch_versions
}
