include { KMC }                 from "./../modules/kmc/main"
include { VG_PATHS }            from "./../modules/vg/paths/main"
include { VG_GIRAFFE }          from "./../modules/vg/giraffe/main"
include { SAMTOOLS_MERGE } 		from "./../modules/samtools/merge/main"
include { SAMTOOLS_INDEX } 		from "./../modules/samtools/index/main"
include { SAMTOOLS_MARKDUP } 	from "./../modules/samtools/markdup/main"
include { SAMTOOLS_SORT }       from "./../modules/samtools/sort/main"
include { SAMTOOLS_FIXMATE }    from "./../modules/samtools/fixmate/main"

ch_versions = Channel.from([])

// To Do:
// Build a new container for VG that includes samtools to run sort and fixmate in one process and avoid
// that the data touches the harddrive 3 or 4 times for no reason. 
workflow VG {

    take:
    reads
    vg_index
    path_name
    ch_fasta
    
    main:

    // Generate Kmer stats
    KMC(
        reads
    )

    // Extract the desired mapping path
    VG_PATHS(
        vg_index.collect(),
        path_name
    )

    ch_versions = ch_versions.mix(VG_PATHS.out.versions)

    // Run the aligner against the mapping path and kmer stats
    VG_GIRAFFE(
        KMC.out.reads,
        VG_PATHS.out.index.collect()
    )

    ch_versions = ch_versions.mix(VG_GIRAFFE.out.versions)
    ch_aligned_bams = VG_GIRAFFE.out.bam

    SAMTOOLS_SORT(
        ch_aligned_bams
    )

	bam_mapped = SAMTOOLS_SORT.out.bam.map { meta, bam ->
        new_meta = [:]
        new_meta.patient_id = meta.patient_id
        new_meta.sample_id = meta.sample_id
        tuple( new_meta, bam)
    }.groupTuple()
            
    bam_mapped.branch {
        single:   it[1].size() == 1
        multiple: it[1].size() > 1
    }.set { bam_to_merge }

    // Merge multi-lane data
    SAMTOOLS_MERGE( 
        bam_to_merge.multiple 
    )

    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)
    
    // Index final alignment file
    SAMTOOLS_INDEX(
        SAMTOOLS_MERGE.out.bam.mix( bam_to_merge.single )
    )
    
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    SAMTOOLS_FIXMATE(
        SAMTOOLS_INDEX.out.bam
    )

    ch_versions = ch_versions.mix(SAMTOOLS_FIXMATE.out.versions)
    
    // Perform duplicate marking
    SAMTOOLS_MARKDUP(
        SAMTOOLS_FIXMATE.out.bam,
        ch_fasta
    )
    
    ch_versions = ch_versions.mix(SAMTOOLS_MARKDUP.out.versions)

    emit:
    bam = SAMTOOLS_MARKDUP.out.bam
    stats = SAMTOOLS_MARKDUP.out.report
    versions = ch_versions
}
