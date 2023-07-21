
process PICARD_WGS_METRICS {

    tag "${meta.patient_id}|${meta.sample_id}"

    publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/QC", mode: 'copy'

    container "docker://quay.io/biocontainers/picard:2.26.11--hdfd78af_0"

    input:
    tuple val(meta),path(bam),path(bai)
    tuple path(fasta),path(fai),path(dict)
    path(bed)

    output:
    path(picard_stats), emit: stats
    
    script:
    base_name = bam.getBaseName()
    picard_stats = base_name + "_wgs_metrics.txt"
    intervals = bed.getBaseName() + ".interval_list"

    """
    picard BedToIntervalList I=$bed O=$intervals SD=${dict}
    picard CollectWgsMetrics I=$bam REFERENCE_SEQUENCE=${fasta} O=$picard_stats INTERVALS=$intervals
    """
}

