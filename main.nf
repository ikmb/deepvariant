

if (!params.genome) {
	exit 1, "Must provide a genome assembly name!"
}
if (!params.samples) {
	exit 1, "Must provide a samples file in CSV format"
}

params.fasta = params.genomes[params.genome].fasta
params.dict = params.genomes[params.genome].dict

params.bed = params.genomes[params.genome].bed

params.fai = params.genomes[params.genome].fai
params.fastagz = params.genomes[params.genome].fastagz
params.gzfai = params.genomes[params.genome].gzfai
params.gzi = params.genomes[params.genome].gzi

params.dbsnp = params.genomes[params.genome].dbsnp

if (!params.fasta || !params.fai || !params.dict || !params.fastagz || !params.gzfai || !params.gzi) {
	exit 1, "Missing one or several mandatory options..."
}

Channel
	.fromPath(params.samples)
	.splitCsv(sep: ';', header: true)
	.set { inputFastp }

Channel
	.fromPath(params.bed)
	.ifEmpty {exit 1; "Could not find a BED file"}
	.into { BedFile; BedToMerge; BedCoverage }

faiToExamples = Channel
    .fromPath(params.fai)
    .ifEmpty{exit 1, "Fai file not found: ${params.fai}"}

fastaGz = Channel
    .fromPath(params.fastagz)
    .ifEmpty{exit 1, "Fastagz file not found: ${params.fastagz}"}
    .into {fastaGzToExamples; fastaGzToVariants }

gzFai = Channel
    .fromPath(params.gzfai)
    .ifEmpty{exit 1, "gzfai file not found: ${params.gzfai}"}
    .into{gzFaiToExamples; gzFaiToVariants }

gzi = Channel
    .fromPath(params.gzi)
    .ifEmpty{exit 1, "gzi file not found: ${params.gzi}"}
    .into {gziToExamples; gziToVariants}

process runFastp {

	label 'fastp'

        scratch true

        input:
	set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, center, run_date, fastqR1, fastqR2 from inputFastp

        output:
	set val(indivID), val(sampleID), val(libraryID), val(rgID), val(platform_unit), val(platform), val(platform_model), val(center), val(run_date),file(left),file(right) into inputBwa
	set file(json),file(html) into outputReportTrimming

        script:
	left = file(fastqR1).getBaseName() + ".trimmed.fastq.gz"
	right = file(fastqR2).getBaseName() + ".trimmed.fastq.gz"
	json = indivID + "_" + sampleID + "_" + libraryID + "-" + rgID + ".fastp.json"
	html = indivID + "_" + sampleID + "_" + libraryID + "-" + rgID + ".fastp.html"
        """
		fastp --in1 $fastqR1 --in2 $fastqR2 --out1 $left --out2 $right --detect_adapter_for_pe -w ${task.cpus} -j $json -h $html
        """

}

process runBwa {

	label 'bwa'

        scratch true

        input:
	set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, center, run_date,file(fastqR1),file(fastqR2) from inputBwa

        output:
	set indivID, sampleID, file(outfile) into runBWAOutput
	val(sample_name) into SampleNames

        script:
	outfile = indivID + "_" + sampleID + "." + libraryID + "_" + rgID + ".aligned.bam"
	sample_name = "${indivID}_${sampleID}"
        """
		bwa mem -H $params.dict -M -R "@RG\\tID:${rgID}\\tPL:ILLUMINA\\tPU:${platform_unit}\\tSM:${indivID}_${sampleID}\\tLB:${libraryID}\\tDS:${params.fasta}\\tCN:${center}" \
		-t ${task.cpus} ${params.fasta} $fastqR1 $fastqR2 \
		| samtools fixmate -@ 4 -m - - \
                | samtools sort -@ 4 -O bam -o $outfile -
        """

}

process runMD {

	publishDir "${params.outdir}/${indivID}/${sampleID}/", mode: 'copy'

	label 'samtools'

	input:
	set val(indivID), val(sampleID), file(bam) from runBWAOutput

	output:
	set indivID, sampleID, file(bam_md),file(bam_index) into (BamMD,BamStats)
	set file(bam_md),file(bam_index) into BamMDCoverage

	script:
	bam_md = bam.getBaseName() + ".md.cram"
	bam_index = bam_md + ".crai"

	"""
		samtools markdup -O CRAM --reference $params.fasta -@ ${task.cpus} $bam $bam_md
		samtools index $bam_md
	"""

}


process runDeepvariant {

	publishDir "${params.outdir}/${indivID}/${sampleID}/DeepVariant", mode: 'copy'

	label 'deepvariant'

	input:
	set indivID, sampleID, file(bam),file(bai) from BamMD
	file(bed) from BedFile.collect()
	file fai from faiToExamples.collect()
	file fastagz from fastaGzToExamples.collect()
	file gzfai from gzFaiToExamples.collect()
	file gzi from gziToExamples.collect()

	output:
	set indivID,sampleID,file(gvcf) into DvVCF
	file(gvcf) into MergeGVCF
	file(vcf)

	script:
	gvcf = bam.getBaseName() + ".g.vcf.gz"
	vcf = bam.getBaseName() + ".vcf.gz"
	
	"""
		/opt/deepvariant/bin/run_deepvariant \
		--model_type=WGS \
		--ref=$fastagz \
		--reads $bam \
		--output_vcf=$vcf \
		--output_gvcf=$gvcf \
		--regions=$bed \
		--num_shards=${task.cpus}
	"""
}

if (params.joint_calling) {

	process runMergeGvcf {

		scratch true 
        	label 'glnexus'

	        input:
		file(gvcfs) from MergeGVCF.collect()
		file(bed) from BedToMerge.collect()

        	output:
		file(merged_vcf) into MergedVCF

        	script:
		merged_vcf = "deepvariant.merged.vcf.gz"

        	"""
			/usr/local/bin/glnexus_cli \
			--config DeepVariantWGS \
			--bed $bed \
			$gvcfs | bcftools view - | bgzip -c > $merged_vcf
 
	        """
	}

	process annotateIDs {

                publishDir "${params.outdir}/DeepVariant", mode: 'copy'

        	label 'glnexus'

	        input:
        	file (vcf) from MergedVCF

	        output:
		file(vcf_annotated) into VcfAnnotated
		file(vcf_annotated_index)

	        script:
        	vcf_annotated = vcf.getBaseName() + ".rsids.vcf.gz"
		vcf_annotated_index = vcf_annotated + ".tbi"
	        """
			tabix $vcf
                	bcftools annotate -c ID -a $params.dbsnp -O z -o $vcf_annotated $vcf
			tabix $vcf_annotated
	        """
	}

	process VcfGetSample {

                publishDir "${params.outdir}/DeepVariant", mode: 'copy'

		label 'glnexus'

		input:
		file(vcf) from VcfAnnotated
		val(sample_name) from SampleNames

		output:
		set  file(vcf_sample),file(vcf_sample_index) into VcfSample

		script:
		vcf_sample = sample_name + ".vcf.gz"
		vcf_sample_index = vcf_sample + ".tbi"

		"""
			bcftools view -o $vcf_sample -O z -t ${task.cpus} -a -s $sample_name $vcf
			tabix $vcf_sample
		"""

	}

	process VcfStats {

		label 'glnexus'

		input:
		set file(vcf),file(tbi) from VcfSample

		output:
		file(vcf_stats) into VcfInfo

		script:
		vcf_stats = vcf.getBaseName() + ".stats"

		"""
			bcftools stats $vcf > $vcf_stats
		"""

	}

} else {

	VcfInfo = Channel.empty()
}

process runWgsCoverage {

        label 'mosdepth'

        input:
        set file(bam),file(bai) from BamMDCoverage
	file(bed) from BedCoverage.collect()

        output:
        set file(genome_bed_coverage),file(genome_global_coverage) into Coverage

        script:
        base_name = bam.getBaseName()
        genome_bed_coverage = base_name + ".mosdepth.region.dist.txt"
	genome_global_coverage = base_name + ".mosdepth.global.dist.txt"
	

        """
                mosdepth -t ${task.cpus} -n -f $params.fasta -x -Q 10 -b $bed $base_name $bam
        """

}

process runMultiQC {

	label 'multiqc'

        publishDir "${params.outdir}/MultiQC", mode: 'copy'

	input:
	file('*') from outputReportTrimming.collect()
	file('*') from Coverage.collect()
	file('*') from VcfInfo.collect()

	output:
	file("multiqc_report.html")

	script:

	"""
		multiqc .
	"""

}
