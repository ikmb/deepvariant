include { DEEPVARIANT } from "./../modules/deepvariant/main"
include { TABIX } from "./../modules/htslib/tabix"

ch_vcfs = Channel.from([])
ch_versions = Channel.from([])

workflow DEEPVARIANT_SHORT_READS {

    take:
    bam
    bed
    ch_fasta
    
    main:
    
    DEEPVARIANT(
        bam,
        bed.collect(),
        ch_fasta
    )

    ch_versions = ch_versions.mix(DEEPVARIANT.out.versions)

    TABIX(
        DEEPVARIANT.out.vcf
    )
    
    ch_versions = ch_versions.mix(TABIX.out.versions)
    ch_vcfs = ch_vcfs.mix(TABIX.out.vcf)

    emit:
    gvcf = DEEPVARIANT.out.gvcf
    vcf = ch_vcfs
    versions = ch_versions

}

