include { FASTP } from "./../modules/fastp/main.nf"
include { BWA } from "./../modules/bwa/main.nf"
include { SAMTOOLS_MERGE_AND_DEDUP} from "./../modules/samtools/main.nf"

workflow TRIM_AND_ALIGN {

	take:
		reads
		bed
		
	main:
		
		FASTP(reads)
		BWA(FASTP.out.reads)
		SAMTOOLS_MERGE_AND_DEDUP(
			BWA.out.bam.map { m,b ->
				def new_meta = [:]
				new_meta.patient_id = m.patient_id
				new_meta.sample_id = m.sample_id
				tuple(new_meta,b)
			}.groupTuple()
		)

	emit:
		bam = SAMTOOLS_MERGE_AND_DEDUP.out.bam
}

