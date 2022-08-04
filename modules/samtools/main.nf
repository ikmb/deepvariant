process SAMTOOLS_MERGE_AND_DEDUP {

	tag "${meta.patient_id}|${meta.sample_id}"

        publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/Align", mode: 'copy'

        input:
        tuple val(meta), path(aligned_bam_list)

        output:
        tuple val(meta),path(merged_bam),path(merged_bam_index), emit: bam
        tuple path(merged_bam),path(merged_bam_index), emit: bam_simple

        script:
        merged_bam = meta.patient_id + "_" + meta.sample_id + ".merged.md.cram"
        merged_bam_index = merged_bam + ".crai"

        if (aligned_bam_list.size() > 1 && aligned_bam_list.size() < 1000 ) {
                """
                        samtools merge -@ ${task.cpus} tmp.bam ${aligned_bam_list.join(' ')}
                        samtools index tmp.bam
                        samtools markdup -O CRAM --reference $params.fasta -@ ${task.cpus} tmp.bam $merged_bam
                        rm tmp.bam*
                        samtools index $merged_bam

                """
        } else {
                """
                        samtools index $aligned_bam_list
                        samtools markdup -O CRAM --reference $params.fasta -@ ${task.cpus} $aligned_bam_list $merged_bam
                        samtools index $merged_bam
                """
        }
}

process SAMTOOLS_MERGE_BAM {

	tag "${meta.patient_id}|${meta.sample_id}"

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/Align", mode: 'copy'

	input:
        tuple val(meta), path(bams),path(bais)

        output:
        tuple val(meta),path(merged_bam),path(merged_bam_index), emit: bam
	tuple path(merged_bam),path(merged_bam_index)

        script:
        merged_bam = meta.patient_id + "_" + meta.sample_id + ".merged.cram"
        merged_bam_index = merged_bam + ".crai"

	"""
		samtools merge -O CRAM --reference $params.fasta -@ ${task.cpus} $merged_bam $bams
		samtools index $merged_bam
	"""

}

