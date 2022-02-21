process merge_and_dedup {

        publishDir "${params.outdir}/${indivID}/${sampleID}/Align", mode: 'copy'

        input:
        tuple val(indivID), val(sampleID), path(aligned_bam_list)

        output:
        tuple val(indivID),val(sampleID),path(merged_bam),path(merged_bam_index), emit: bam
        tuple path(merged_bam),path(merged_bam_index), emit: bam_simple

        script:
        merged_bam = indivID + "_" + sampleID + ".merged.md.cram"
        merged_bam_index = merged_bam + ".crai"
        sample_name = indivID + "_" + sampleID

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

process merge_bam_files {

	publishDir "${params.outdir}/${indivID}/${sampleID}/Align", mode: 'copy'

	input:
        tuple val(indivID), val(sampleID), path(bams),path(bais)

        output:
        tuple val(indivID),val(sampleID),path(merged_bam),path(merged_bam_index)
	tuple path(merged_bam),path(merged_bam_index)

        script:
        merged_bam = indivID + "_" + sampleID + ".merged.cram"
        merged_bam_index = merged_bam + ".crai"

	"""
		samtools merge -O CRAM --reference $params.fasta -@ ${task.cpus} $merged_bam $bams
		samtools index $merged_bam
	"""

}

