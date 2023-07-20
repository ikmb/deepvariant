nextflow.enable.dsl=2

fasta       = params.genomes[params.assembly].fasta
fai         = params.genomes[params.assembly].fai
dict        = params.genomes[params.assembly].dict

ch_fasta    = Channel.fromList( [ file(fasta , checkIfExists: true), file(fai, checkIfExits: true), file(dict, checkIfExists: true) ] ).collect()

bwa2_index  = params.genomes[params.assembly].bwa2_index
mmi_index   = params.genomes[params.assembly].mmi

params.bed  = params.intervals ?: params.genomes[params.assembly].bed

tools = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase().replaceAll('-', '').replaceAll('_', '')} : []

if (params.pacbio) {

    params.dv_model = "PACBIO"

    Channel.fromPath(params.samples)
        .splitCsv(sep: ';', header: true)
        .map { create_pacbio_channel(it) }
        .set { reads }

    Channel.fromPath(params.tandem_repeats)
        .ifEmpty { exit 1; "Missing tandem repeats for pacbio SV calling..." }
        .set { tandem_repeats }

} else {

    params.dv_model = "WGS"

    Channel.fromPath(params.samples)
        .splitCsv(sep: ';', header: true)
        .map { create_fastq_channel(it) }
        .set { reads }
}

// Import workflows
include { DEEPVARIANT_SHORT_READS }                     from "../subworkflows/deepvariant_short_reads"
include { DEEPVARIANT_LONG_READS }                      from "../subworkflows/deepvariant_long_reads"
include { PB_ALIGN }                                    from "../subworkflows/pb_align"
include { TRIM_AND_ALIGN }                              from "../subworkflows/trim_and_align"
include { PICARD_WGS_METRICS  }                         from "../modules/picard/collect_wgs_metrics"
include { MULTIQC }                                     from "../modules/multiqc/main"
include { MOSDEPTH }                                    from "../modules/mosdepth/main"
include { BAM_SELECT_READS }                            from "../modules/helper/bam_select_reads"
include { MANTA }                                       from "../modules/manta/main"
include { BGZIP_INDEX }                                 from "../modules/htslib/bgzip_index"
include { PBSV_SIG; PBSV_CALL }                         from "../modules/sv/main"
include { CUSTOM_DUMPSOFTWAREVERSIONS }                 from "./../modules/custom/dumpsoftwareversions/main"

// Initialize channels
Channel
    .fromPath(params.bed)
    .ifEmpty {exit 1; "Could not find a BED file"}
    .set { ch_bed }.collect()

ch_vcf      = Channel.from([])
ch_bam      = Channel.from([])
ch_versions = Channel.from([])

workflow DEEPVARIANT_PIPELINE {

    main:

    if (params.pacbio) {

        PB_ALIGN(
            reads,
            mmi
        )
        
        ch_bam = ch_bam.mix(PB_ALIGN.out.bam)

        ch_versions = ch_versions.mix(PB_ALIGN.out.versions)

        if ('deepvariant' in tools) {
            DEEPVARIANT_LONG_READS(
                ch_bam,
                ch_bed,
                ch_fasta
            )
        }

        ch_versions = ch_versions.mix(DEEPVARIANT_LONG_READS.out.versions)
        
        if ('pbsv' in tools) {
            PBSV_SIG(
                ch_bam,
                tandem_repeats.collect()        
            )
            PBSV_CALL(PBSV_SIG.out[1].collect())
        }

        ch_vcf = ch_vcf.mix(DEEPVARIANT_LONG_READS.out.vcf)

    } else {

        TRIM_AND_ALIGN(
            reads,
            bwa2_index,
            ch_fasta
        )

        ch_versions = ch_versions.mix(TRIM_AND_ALIGN.out.versions)

        ch_bam = ch_bam.mix(TRIM_AND_ALIGN.out.bam)

        if ('intersect' in tools) {
            BAM_SELECT_READS(
                TRIM_AND_ALIGN.out.bam,
                ch_bed.collect(),
                ch_fasta
            )

            ch_versions = ch_versions.mix(BAM_SELECT_READS.out.versions)

        }

        if ('deepvariant' in tools) {
            DEEPVARIANT_SHORT_READS(
                TRIM_AND_ALIGN.out.bam,
                ch_bed.collect(),
                ch_fasta
            )
            ch_vcf = ch_vcf.mix(DEEPVARIANT_SHORT_READS.out.vcf)

            ch_versions = ch_versions.mix(DEEPVARIANT_SHORT_READS.out.versions)
        }

        if ('manta' in tools) {
            BGZIP_INDEX(
                ch_bed
            )
            
            ch_versions = ch_versions.mix(BGZIP_INDEX.out.versions)

            MANTA(
                ch_bam,
                BED_COMPRESS_AND_INDEX.out.bed.collect(),
                ch_fasta.collect()
            )

            ch_versions = ch_versions.mix(MANTA.out.versions)

        }
    }

    MOSDEPTH(
        ch_bam,
        ch_fasta,
        ch_bed.collect()
    )
    
    PICARD_WGS_METRICS(
        ch_bam,
        ch_fasta,
        ch_bed.collect()
    )

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    MULTIQC(
        MOSDEPTH.out.mix(
            PICARD_WGS_METRICS.out,CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml
        ).collect()
    )

    emit:
    qc = MULTIQC.out.report

}

def create_fastq_channel(LinkedHashMap row) {

    // IndivID;SampleID;libraryID;rgID;rgPU;platform;platform_model;Center;Date;R1;R2

    def meta = [:]
    meta.patient_id = row.IndivID
    meta.sample_id = row.SampleID
    meta.library_id = row.libraryID
    meta.readgroup_id = row.rgID
    meta.center = row.Center
    meta.date = row.Date
    meta.platform_unit = row.rgPU

    def array = []
    array = [ meta, file(row.R1), file(row.R2) ]

    return array
}

def create_pacbio_channel(LinkedHashMap row) {

    // IndivID;SampleID;libraryID;rgID;rgPU;platform;platform_model;Center;Date;R1;R2

    def meta = [:]
    meta.patient_id = row.IndivID
    meta.sample_id = row.SampleID

    def array = []
    array = [ meta, file(row.R1) ]

    return array
}
