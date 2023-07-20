process BAM_SELECT_READS {

    tag "${meta.patient_id}|${meta.sample_id}"

    container 'quay.io/biocontainers/samtools:1.16.1--h6899075_1'

    label 'medium_serial'
    
    publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/SubsetReads", mode: 'copy'

    input:
    tuple val(meta),path(bam),path(bai)
    path(bed)
    tuple path(fasta),path(fai),path(dict)

    output:
    tuple val(meta),path(R1),path(R2), emit: fastq
    path("versions.yml"), emit: versions

    script:
    R1 = bam.getBaseName() + "_R1_001.fastq.gz"
    R2 = bam.getBaseName() + "_R2_001.fastq.gz"

    """
    samtools view --reference $fasta -hb -o mapped.cram -L $bed $bam
    samtools view --reference $fasta -hb -o unmapped.cram -f 4 $bam
    samtools merge merged.cram mapped.cram unmapped.cram
    samtools collate -u -O merged.cram | samtools fastq -1 $R1 -2 $R2 -0 /dev/null -s /dev/null -n
    rm *mapped.cram merged.cram

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}