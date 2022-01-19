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

params.bed = params.genomes[params.genome].bed

params.fai = params.genomes[params.genome].fai
params.fastagz = params.genomes[params.genome].fastagz
params.gzfai = params.genomes[params.genome].gzfai
params.gzi = params.genomes[params.genome].gzi

params.dbsnp = params.genomes[params.genome].dbsnp

params.mitochondrion = params.genomes[params.genome].mitochondrion

params.tandem_repeats = params.genomes[params.genome].tandem_repeats

if (!params.fasta || !params.fai || !params.dict || !params.fastagz || !params.gzfai || !params.gzi) {
	exit 1, "Missing one or several mandatory options..."
}

if (params.pacbio) {

	params.dv_model = "PACBIO"

	Channel.fromPath(params.samples)
		.splitCsv(sep: ';', header: true)
		.map{ row-> tuple(row.IndivID,row.SampleID,file(row.R1) ) }
		.set { reads }

	Channel.fromPath(params.tandem_repeats)
		.ifEmpty { exit 1; "Missing tandem repeats for pacbio SV calling..." }
		.set { tandem_repeats }

} else {

	params.dv_model = "WGS"

	Channel.fromPath(params.samples)
		.splitCsv(sep: ';', header: true)
		.map{ row-> tuple( row.IndivID,row.SampleID,row.libraryID,row.rgID,row.platform_unit,row.platform,row.platform_model,row.center,row.run_date,file(row.R1),file(row.R1)  ) }
		.set { reads }
}

// Import workflows
include { DEEPVARIANT_SHORT_READS } from "./workflows/deepvariant_illumina/main.nf" params(params)
include { DEEPVARIANT_PACBIO } from "./workflows/deepvariant_pacbio/main.nf" params(params)
include { multiqc ; wgs_coverage } from "./modules/qc/main.nf" params(params)
include { vcf_stats } from "./modules/vcf/main.nf" params(params)
include { vep } from "./modules/vep/main.nf" params(params)

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
	} else {
		DEEPVARIANT_ILLUMINA(reads,bed,fastaGz,gzFai,gzi,fai)
		bam = DEEPVARIANT_ILLUMINA.out.bam
		vcf = DEEPVARIANT_ILLUMINA.out.vcf
		gvcf = DEEPVARIANT_ILLUMINA.out.gvcf
	}

	if (params.vep) {
		vep(vcf)
	}

	if (params.joint_calling) {

	}

	wgs_coverage(bam,bed)
	vcf_stats(vcf)

	multiqc(wgs_coverage.out.concat(vcf_stats.out).collect())

}
