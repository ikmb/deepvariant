include { KMC }                 from "./../modules/kmc/main"
include { VG_PATHS }            from "./../modules/vg/paths/main"
include { VG_GIRAFFE }          from "./../modules/vg/giraffe/main"
include { SAMTOOLS_MERGE } 		from "./../modules/samtools/merge/main"
include { SAMTOOLS_INDEX } 		from "./../modules/samtools/index/main"
include { SAMTOOLS_MARKDUP } 	from "./../modules/samtools/markdup/main"

ch_versions = Channel.from([])

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

	bam_mapped = ch_aligned_bams.map { meta, bam ->
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
    
    // Perform duplicate marking
    SAMTOOLS_MARKDUP(
        SAMTOOLS_INDEX.out.bam,
        ch_fasta
    )
                                                                    
    emit:
    bam = SAMTOOLS_MARKDUP.out.bam
    stats = SAMTOOLS_MARKDUP.out.report
    versions = ch_versions
}