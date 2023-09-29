process DEEPVARIANT {

    tag "${meta.patient_id}|${meta.sample_id}"

    container "docker://google/deepvariant:1.5.0"

    publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/DeepVariant", mode: 'copy'

    label 'long_parallel'

    scratch true

    input:
    tuple val(meta), path(bam),path(bai)
    path(bed)
    tuple path(fasta),path(fai),path(dict)

    output:
    tuple val(meta),path(gvcf), emit: gvcf
    tuple val(meta),path(vcf), emit: vcf
    val(meta), emit: sample
    path("versions.yml"), emit: versions

    script:
    gvcf = bam.getBaseName() + ".g.vcf.gz"
    vcf = bam.getBaseName() + ".vcf.gz"

    def model = "WGS"
    def options = ""
    if (params.pacbio) {
        model = "PACBIO"
        if (params.intervals) {
            options = "--regions=${bed}"
        }
    } else {
        if (bed)
            options = "--regions=${bed}"
        end
    }
    """
    set TMPDIR=\$PWD

    /opt/deepvariant/bin/run_deepvariant \
    --model_type=${model} \
    --ref=$fasta \
    --reads $bam \
    --output_vcf=$vcf \
    --output_gvcf=$gvcf \
    $options \
    --num_shards=${task.cpus} \
    --intermediate_results_dir tmp_data

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: \$(echo \$(/opt/deepvariant/bin/run_deepvariant --version) | sed 's/^.*version //; s/ .*\$//' )
    END_VERSIONS
    """
}

