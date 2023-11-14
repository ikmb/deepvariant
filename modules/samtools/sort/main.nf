process SAMTOOLS_SORT {

    container 'quay.io/biocontainers/samtools:1.16.1--h6899075_1'

    label 'short_parallel'
    
    tag "${meta.patient_id}|${meta.sample_id}"

    input:
    tuple val(meta), path(b)

    output:
    tuple val(meta), path(sorted_bam), emit: bam
    path("versions.yml"), emit: versions

    script:
    
    sorted_bam = b.getBaseName() + ".sorted.bam"

    """
    samtools sort -@ ${task.cpus} -m 2G -O bam -o $sorted_bam  $b

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}
