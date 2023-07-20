include { FASTP }				from "./../modules/fastp/main"
include { BWA2_MEM } 			from "./../modules/bwa2/mem/main"
include { SAMTOOLS_MERGE } 		from "./../modules/samtools/merge/main"
include { SAMTOOLS_INDEX } 		from "./../modules/samtools/index/main"
include { SAMTOOLS_MARKDUP } 	from "./../modules/samtools/markdup/main"

ch_versions = Channel.from([])

workflow TRIM_AND_ALIGN {

	take:
		reads
		bwa_index
		ch_fasta

	main:
		
		FASTP(
			reads
		)

		ch_versions = ch_versions.mix(FASTP.out.versions)

		BWA2_MEM(
			FASTP.out.reads,
			ch_fasta,
			bwa_index
		)

		ch_versions = ch_versions.mix(BWA2_MEM.out.versions)

        ch_aligned_bams = BWA2_MEM.out.bam

		bam_mapped = ch_aligned_bams.map { meta, bam ->
            new_meta = [:]
            new_meta.patient_id = meta.patient_id
            new_meta.sample_id = meta.sample_id
            def groupKey = meta.sample_id
            tuple( groupKey, new_meta, bam)
        }.groupTuple(by: [0,1]).map { g ,new_meta ,bam -> [ new_meta, bam ] }
            
        bam_mapped.branch {
            single:   it[1].size() == 1
            multiple: it[1].size() > 1
        }.set { bam_to_merge }

        SAMTOOLS_MERGE( 
			bam_to_merge.multiple 
		)

		ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

        SAMTOOLS_INDEX(
			SAMTOOLS_MERGE.out.bam.mix( bam_to_merge.single )
		)

		ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

        SAMTOOLS_MARKDUP(
			SAMTOOLS_INDEX.out.bam,
			ch_fasta
		)

		ch_versions = ch_versions.mix(SAMTOOLS_MARKDUP.out.versions)

	emit:
		bam = SAMTOOLS_MARKDUP.out.bam
		versions = ch_versions
}

