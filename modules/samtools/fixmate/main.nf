process SAMTOOLS_FIXMATE {

    container 'quay.io/biocontainers/samtools:1.16.1--h6899075_1'

    label 'short_parallel'
    
    tag "${meta.patient_id}|${meta.sample_id}"

    input:
    tuple val(meta), path(b),path(i)

    output:
    tuple val(meta), path(fixed_bam),path(fixed_bai), emit: bam
    path("versions.yml"), emit: versions

    script:
    fixed_bam = b.getBaseName() + ".fm.bam"
    fixed_bai = fixed_bam + ".bai"
    
    """
    samtools sort -n $b | samtools fixmate -@ ${task.cpus} -m -O bam - $fixed_bam
    samtools index $fixed_bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}
