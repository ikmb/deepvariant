process glnexus {

	scratch true
	label 'glnexus'

	input:
	path(gvcfs)
	path(bed)

	output:
	path(merged_vcf)

	script:
	 merged_vcf = "deepvariant.merged.vcf.gz"

	"""
	/usr/local/bin/glnexus_cli \
		--config DeepVariantWGS \
		--bed $bed \
		$gvcfs | bcftools view - | bgzip -c > $merged_vcf

	"""
}

process vcf_add_dbsnp {

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

process vcf_get_sample {

	publishDir "${params.outdir}/${indivID}/${sampleID}/${folder}", mode: 'copy'

	label 'glnexus'

	input:
	tuple val(indivID),val(sampleID),path(vcf),path(tbi)
	val(folder)

	output:
	tuple val(indivID),val(sampleID),path(vcf_sample),path(vcf_sample_index)

	script:
	vcf_sample = vcf.getSimpleName() + "." + sampleID  + ".vcf.gz"
	vcf_sample_index = vcf_sample + ".tbi"

	"""
		bcftools view -o $vcf_sample -O z -a -s $sampleID $vcf
		tabix $vcf_sample
	"""

}

process vcf_stats {

	label 'glnexus'

	input:
	tuple val(indivID),val(sampleID),path(vcf),path(tbi)

	output:
	path(vcf_stats)

	script:
	vcf_stats = vcf.getBaseName() + ".stats"

	"""
		bcftools stats $vcf > $vcf_stats
	"""

}
	
process vcf_index {
	
	input:
	tuple val(indivID),val(sampleID),path(vcf)

	output:
	tuple val(indivID),val(sampleID),path(vcf),path(tbi)

	script:
	tbi = vcf + ".tbi"

	"""
		tabix $vcf
	"""

}

process vcf_pass {

	publishDir "${params.outdir}/${indivID}/${sampleID}/Variants", mode: 'copy'
	input:
	tuple val(indivID),val(sampleID),path(vcf),path(tbi)

	output:
	tuple val(indivID),val(sampleID),path(vcf_f),path(tbi_f)

	script:
	vcf_f = vcf.getSimpleName() + ".pass.vcf.gz"
	tbi_f = vcf_f + ".tbi"

	"""
		bcftools view -f "PASS" -O z -o $vcf_f $vcf
		tabix $vcf_f
	"""

}
	
