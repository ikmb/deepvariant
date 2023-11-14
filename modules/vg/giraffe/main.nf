process VG_GIRAFFE {

    container 'quay.io/biocontainers/vg:1.52.0--h9ee0642_0'

    publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/VG", mode: 'copy'

    label 'long_parallel'

    tag "${meta.patient_id}|${meta.sample_id}"

    input:
    tuple val(meta),path(R1),path(R2),path(kff)
    tuple path(gbz),path(hapl),path(paths)

    output:
    tuple val(meta),path(bam), emit: bam
    path("versions.yml"), emit: versions

    script:
    bam = meta.sample_id + "-" + meta.readgroup_id + ".vg.bam"

    """
    vg giraffe --progress \
    --read-group "ID:${meta.readgroup_id} LB:${meta.library_id} SM:${meta.sample_id} PL:illumina PU:${meta.platform_unit}" \
    --sample ${meta.sample_id} \
    -o BAM \
    --ref-paths $paths \
    -P \
    -L  3000 \
    -f $R1 -f $R2 \
    -Z $gbz  \
    --kff-name $kff \
    --haplotype-name $hapl \
    -t ${task.cpus} > $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: 1.52
    END_VERSIONS
    """

}
