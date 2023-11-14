process KMC {

    container 'quay.io/biocontainers/kmc:3.2.1--hf1761c0_2'

    label 'medium_parallel'

    tag "${meta.patient_id}|${meta.sample_id}"

    input:
    tuple val(meta),path(fastqR1),path(fastqR2)

    output:
    tuple val(meta),path(fastqR1),path(fastqR2),path(kff), emit: reads
    path("versions.yml"), emit: versions

    script:
    kff = meta.sample_id + ".kff"

    """
    for i in \$(echo *.fastq.gz); do echo \$i >> reads.paths; done;
    kmc -k29 -m${task.memory.toGiga()} -okff -t${task.cpus} @reads.paths ${meta.sample_id} \$TMPDIR

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kmc: 3.2.1
    END_VERSIONS
    """

}
