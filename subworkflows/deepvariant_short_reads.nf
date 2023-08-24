include { DEEPVARIANT } from "./../modules/deepvariant/main"
include { TABIX } from "./../modules/htslib/tabix"

ch_vcfs = Channel.from([])

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

    TABIX(
        DEEPVARIANT.out.vcf
    )
    
    ch_vcfs = ch_vcfs.mix(TABIX.out.vcf)

    emit:
    gvcf = DEEPVARIANT.out.gvcf
    vcf = ch_vcfs

}

