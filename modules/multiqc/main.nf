process MULTIQC {

    container 'quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0'

    publishDir "${params.outdir}/Summary/", mode: 'copy'

    input:
    path('*')

    output:
    path("${params.run_name}_multiqc.html"), emit: report
    path("versions.yml"), emit: versions

    script:

    """
    cp $baseDir/conf/multiqc_config.yaml multiqc_config.yaml
    multiqc -n ${params.run_name}_multiqc *

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS

    """
}


