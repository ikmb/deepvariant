nextflow.enable.dsl=2

if (!params.genome) {
	exit 1, "Must provide a genome assembly name!"
}
if (!params.samples) {
	exit 1, "Must provide a samples file in CSV format"
}

params.fasta = params.genomes[params.genome].fasta
params.bwa2_index = params.genomes[params.genome].bwa2_index
params.dict = params.genomes[params.genome].dict

params.mmi = params.genomes[params.genome].mmi

params.bed = params.intervals ?: params.genomes[params.genome].bed

params.fai = params.genomes[params.genome].fai
params.fastagz = params.genomes[params.genome].fastagz
params.gzfai = params.genomes[params.genome].gzfai
params.gzi = params.genomes[params.genome].gzi

params.dbsnp = params.genomes[params.genome].dbsnp

params.mitochondrion = params.genomes[params.genome].mitochondrion

params.tandem_repeats = params.genomes[params.genome].tandem_repeats

params.cnv_annotation = file(params.genomes[ params.genome ].cnv_annotation )
params.cnv_mappable = file(params.genomes[ params.genome ].cnv_mappable)
params.cnv_blacklist = file(params.genomes[ params.genome ].cnv_blacklist )
params.cnv_exclusion = file(params.genomes[ params.genome ].cnv_exclusion )

if (!params.fasta || !params.fai || !params.dict || !params.fastagz || !params.gzfai || !params.gzi) {
	exit 1, "Missing one or several mandatory options..."
}

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

bed = Channel.fromPath(params.bed)

// Import workflows
include { DEEPVARIANT_SHORT_READS } from "./workflows/deepvariant_illumina/main.nf"
include { DEEPVARIANT_PACBIO } from "./workflows/deepvariant_pacbio/main.nf"
include { CNVKIT } from "./workflows/cnvkit/main.nf" 
include { MULTIQC ; MOSDEPTH ; PICARD_WGS_METRICS  } from "./modules/qc/main.nf"
include { VCF_STATS } from "./modules/vcf/main.nf"
include { VEP } from "./modules/vep/main.nf"

// Initialize channels
Channel
	.fromPath(params.bed)
	.ifEmpty {exit 1; "Could not find a BED file"}
	.set { bed }

fai = Channel
    .fromPath(params.fai)
    .ifEmpty{exit 1, "Fai file not found: ${params.fai}"}

fastaGz = Channel
    .fromPath(params.fastagz)
    .ifEmpty{exit 1, "Fastagz file not found: ${params.fastagz}"}

gzFai = Channel
    .fromPath(params.gzfai)
    .ifEmpty{exit 1, "gzfai file not found: ${params.gzfai}"}

gzi = Channel
    .fromPath(params.gzi)
    .ifEmpty{exit 1, "gzi file not found: ${params.gzi}"}

// The main pipeline logic

workflow {

	main:

	if (params.pacbio) {
		DEEPVARIANT_PACBIO(reads,bed,tandem_repeats,fastaGz,gzFai,gzi,fai)
		bam = DEEPVARIANT_PACBIO.out.bam
		vcf = DEEPVARIANT_PACBIO.out.vcf
		gvcf = DEEPVARIANT_PACBIO.out.gvcf
		vcf_dv = DEEPVARIANT_PACBIO.out.vcf_dv
	} else {
		DEEPVARIANT_SHORT_READS(reads,bed,fastaGz,gzFai,gzi,fai)
		bam = DEEPVARIANT_SHORT_READS.out.bam
		vcf = DEEPVARIANT_SHORT_READS.out.vcf
		gvcf = DEEPVARIANT_SHORT_READS.out.gvcf
		vcf_dv = DEEPVARIANT_SHORT_READS.out.vcf_dv
	}

	// effect prediction
	if (params.vep) {
		VEP(vcf)
	}

	MOSDEPTH(bam,bed.collect())
	PICARD_WGS_METRICS(bam,bed.collect())
	VCF_STATS(vcf)

	MULTIQC(
		MOSDEPTH.out.mix(
			VCF_STATS.out.stats,PICARD_WGS_METRICS.out
		).collect()
	)

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
