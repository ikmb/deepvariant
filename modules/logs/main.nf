process get_software_versions {

	executor 'local'

	publishDir "${params.outdir}/Summary/versions", mode: 'copy'

	output:
	file("v*.txt")
	file(yaml_file) into software_versions_yaml

	script:
	yaml_file = "software_versions_mqc.yaml"

	"""
		echo $workflow.manifest.version &> v_ikmb_deepvariant.txt
		echo $workflow.nextflow.version &> v_nextflow.txt
		fastp -v &> v_fastp.txt
		echo "Deepvariant 1.2.0" &> v_deepvariant.txt
		echo "GLNexus 1.3.1" &> v_glnexus.txt
		samtools --version &> v_samtools.txt
		multiqc --version &> v_multiqc.txt
		bwa-mem2 > v_bwa.txt 2>&1 || true
		parse_versions.pl >  $yaml_file
	"""
}


