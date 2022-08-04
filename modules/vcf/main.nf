process GLNEXUS {

	scratch true
	label 'glnexus'

	input:
	path(gvcfs)
	path(bed)

	output:
	path(merged_vcf), emit: vcf

	script:
	 merged_vcf = "deepvariant.merged.vcf.gz"

	"""
	/usr/local/bin/glnexus_cli \
		--config DeepVariantWGS \
		--bed $bed \
		$gvcfs | bcftools view - | bgzip -c > $merged_vcf

	"""
}

process VCF_ADD_DBSNP {

	publishDir "${params.outdir}/${indivID}/${sampleID}/Variants", mode: 'copy'

	input:
	tuple val(indivID),val(sampleID),path(vcf),path(tbi)

	output:
	tuple val(indivID),val(sampleID),path(vcf_annotated),path(vcf_annotated_index), emit: vcf

	script:
	vcf_annotated = vcf.getBaseName() + ".rsids.vcf.gz"
	vcf_annotated_index = vcf_annotated + ".tbi"

	"""
		bcftools annotate -c ID -a $params.dbsnp -O z -o $vcf_annotated $vcf
		tabix $vcf_annotated
	"""
}

process VCF_GET_SAMPLE {

	publishDir "${params.outdir}/${indivID}/${sampleID}/${folder}", mode: 'copy'

	label 'glnexus'

	input:
	tuple val(indivID),val(sampleID),path(vcf),path(tbi)
	val(folder)

	output:
	tuple val(indivID),val(sampleID),path(vcf_sample),path(vcf_sample_index), emit: vcf

	script:
	vcf_sample = vcf.getSimpleName() + "." + sampleID  + ".vcf.gz"
	vcf_sample_index = vcf_sample + ".tbi"

	"""
		bcftools view -o $vcf_sample -O z -a -s $sampleID $vcf
		tabix $vcf_sample
	"""

}

process VCF_STATS {

	label 'glnexus'

	input:
	tuple val(indivID),val(sampleID),path(vcf),path(tbi)

	output:
	path(vcf_stats), emit: stats

	script:
	vcf_stats = vcf.getBaseName() + ".stats"

	"""
		bcftools stats $vcf > $vcf_stats
	"""

}
	
process VCF_INDEX {
	
	input:
	tuple val(indivID),val(sampleID),path(vcf)

	output:
	tuple val(indivID),val(sampleID),path(vcf),path(tbi), emit: vcf

	script:
	tbi = vcf + ".tbi"

	"""
		tabix $vcf
	"""

}

process VCF_PASS {

	publishDir "${params.outdir}/${indivID}/${sampleID}/Variants", mode: 'copy'
	input:
	tuple val(indivID),val(sampleID),path(vcf),path(tbi)

	output:
	tuple val(indivID),val(sampleID),path(vcf_f),path(tbi_f), emit: vcf

	script:
	vcf_f = vcf.getSimpleName() + ".pass.vcf.gz"
	tbi_f = vcf_f + ".tbi"

	"""
		bcftools view -f "PASS" -O z -o $vcf_f $vcf
		tabix $vcf_f
	"""

}
	
