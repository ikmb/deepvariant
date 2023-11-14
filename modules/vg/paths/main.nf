process VG_PATHS {

    tag "${path_name}"
    
    container 'quay.io/biocontainers/vg:1.52.0--h9ee0642_0'

    label 'medium_parallel'

    input:
    tuple path(gbz),path(hapl)
    val(path_name)
    
    output:
    tuple path(gbz),path(hapl),path(paths), emit: index
    path("versions.yml"), emit: versions

    script:
    paths = path_name + ".path_list.txt"

    """
    vg paths -x $gbz -L -Q $path_name > $paths
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: 1.52
    END_VERSIONS
    """
}
