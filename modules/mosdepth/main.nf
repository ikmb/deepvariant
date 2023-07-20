process MOSDEPTH {

    tag "${meta.patient_id}|${meta.sample_id}"

    container "quay.io/biocontainers/mosdepth:0.2.9--hbeb723e_1"
    
    publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/QC", mode: 'copy'
    
    label 'medium_serial'

    input:
    tuple val(meta),path(bam),path(bai)
    tuple path(fasta),path(fai),path(dict)
    path(bed)

    output:
    tuple path(genome_bed_coverage),path(genome_global_coverage), emit: coverage

    script:
    base_name = bam.getBaseName()
    genome_bed_coverage = base_name + ".mosdepth.region.dist.txt"
    genome_global_coverage = base_name + ".mosdepth.global.dist.txt"

    """
    mosdepth -t ${task.cpus} -n -f $fasta -x -Q 10 -b $bed $base_name $bam
    """
}