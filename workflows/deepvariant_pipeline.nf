nextflow.enable.dsl=2

fasta       = params.genomes[params.assembly].fasta
fai         = params.genomes[params.assembly].fai
dict        = params.genomes[params.assembly].dict

ch_fasta    = Channel.fromList( [ file(fasta , checkIfExists: true), file(fai, checkIfExits: true), file(dict, checkIfExists: true) ] ).collect()

bwa2_index  = params.genomes[params.assembly].bwa2_index
mmi_index   = params.genomes[params.assembly].mmi

ch_bed      = params.intervals ? Channel.fromPath(params.intervals).collect() : Channel.fromPath(params.genomes[params.assembly].bed).collect()

tools       = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase().replaceAll('-', '').replaceAll('_', '')} : []

// Make sure this assembly as a configured graph reference
if (params.short_read_aligner == "vg" ) {

    if (params.genomes[params.assembly].vg_path) {
        gbz = file(params.genomes["gbz"], checkIfExists: true)
        hapl = file(params.genomes["hapl"], checkIfExists: true)
        ch_vg_index = Channel.from([gbz,hapl])
        vg_path = params.genomes[params.assembly].vg_path
    } else {
        log.info "Requested VG aligner, but no graph index is configured!"
        System.exit(1)
    }

} else {
    ch_vg_index = Channel.empty()
    vg_path = null
}

// Data is from Pacbio long reads
if (params.pacbio) {

    params.dv_model = "PACBIO"

    Channel.fromPath(params.samples)
        .splitCsv(sep: ';', header: true)
        .map { create_pacbio_channel(it) }
        .set { reads }

    Channel.fromPath(params.tandem_repeats)
        .ifEmpty { exit 1; "Missing tandem repeats for pacbio SV calling..." }
        .set { tandem_repeats }

// Data is from regular paired-end short reads
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
include { BWA2_ALIGN }                                  from "../subworkflows/bwa2_align"
include { PICARD_WGS_METRICS  }                         from "../modules/picard/collect_wgs_metrics"
include { MULTIQC }                                     from "../modules/multiqc/main"
include { MOSDEPTH }                                    from "../modules/mosdepth/main"
include { BAM_SELECT_READS }                            from "../modules/helper/bam_select_reads"
include { MANTA }                                       from "../modules/manta/main"
include { BGZIP_INDEX }                                 from "../modules/htslib/bgzip_index"
include { PBSV_SIG; PBSV_CALL }                         from "../modules/sv/main"
include { CUSTOM_DUMPSOFTWAREVERSIONS }                 from "./../modules/custom/dumpsoftwareversions/main"
include { FASTP }				                        from "./../modules/fastp/main"
include { VG }                                          from "../subworkflows/vg"

// Set default channels
ch_vcf      = Channel.from([])
ch_bam      = Channel.from([])
ch_versions = Channel.from([])
ch_reports  = Channel.from([])

workflow DEEPVARIANT_PIPELINE {

    main:

    // Call variants from Pacbio long reads
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

    // Short-read workflow
    } else {

        // Trim reads
        FASTP(
            reads
        )

        ch_versions = ch_versions.mix(FASTP.out.versions)
        ch_reports = ch_reports.mix(FASTP.out.json)

        // Align with VG
        if (params.short_read_aligner == "vg") {
            VG(
                FASTP.out.reads,
                ch_vg_index,
                vg_path,
                ch_fasta
            )
            ch_versions = ch_versions.mix(VG.out.versions)
            ch_reports  = ch_reports.mix(VG.out.stats)
            ch_bam = ch_bam.mix(VG.out.bam)
        } else if (params.short_read_aligner == "bwa2" ) {
            BWA2_ALIGN(
                FASTP.out.reads,
                bwa2_index,
                ch_fasta
            )
            ch_reports  = ch_reports.mix(BWA2_ALIGN.out.stats)
            ch_versions = ch_versions.mix(BWA2_ALIGN.out.versions)
            ch_bam = ch_bam.mix(BWA2_ALIGN.out.bam)
        }

        // Create a BAM file that only covers the regions provided in the BED file (instead of the whole genome)
        if ('intersect' in tools) {
            BAM_SELECT_READS(
                ch_bam,
                ch_bed.collect(),
                ch_fasta
            )

            ch_versions = ch_versions.mix(BAM_SELECT_READS.out.versions)

        }

        // Call variants using DeepVariant
        if ('deepvariant' in tools) {
            DEEPVARIANT_SHORT_READS(
                ch_bam,
                ch_bed.collect(),
                ch_fasta
            )
            ch_vcf = ch_vcf.mix(DEEPVARIANT_SHORT_READS.out.vcf)

            ch_versions = ch_versions.mix(DEEPVARIANT_SHORT_READS.out.versions)
        }

        // Call structural variants using Manta
        if ('manta' in tools) {
            BGZIP_INDEX(
                ch_bed
            )
            
            ch_versions = ch_versions.mix(BGZIP_INDEX.out.versions)

            MANTA(
                ch_bam,
                BED_COMPRESS_AND_INDEX.out.bed.collect(),
                ch_fasta
            )

            ch_versions = ch_versions.mix(MANTA.out.versions)

        }
    }

    // Fast WGS Coverage with MosDepth
    MOSDEPTH(
        ch_bam,
        ch_fasta,
        ch_bed.collect()
    )
    
    // Picard WGS metrics
    PICARD_WGS_METRICS(
        ch_bam,
        ch_fasta,
        ch_bed.collect()
    )
    
    ch_reports = ch_reports.mix(PICARD_WGS_METRICS.out.stats)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    ch_reports = ch_reports.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml)
    
    MULTIQC(
        ch_reports.collect()
    )

    emit:
    qc = MULTIQC.out.report

}

def create_fastq_channel(LinkedHashMap row) {

    // patient;sample;library;rgid;rgpu;R1;R2

    def meta = [:]
    meta.patient_id = row.patient
    meta.sample_id = row.sample
    meta.library_id = row.library
    meta.readgroup_id = row.rgid
    meta.platform_unit = row.rgpu

    def array = []
    array = [ meta, file(row.R1,  checkIfExists: true), file(row.R2, checkIfExists: true) ]

    return array
}

def create_pacbio_channel(LinkedHashMap row) {

    // patient;sample;R1

    def meta = [:]
    meta.patient_id = row.patient
    meta.sample_id = row.sample

    def array = []
    array = [ meta, file(row.R1) ]

    return array
}
