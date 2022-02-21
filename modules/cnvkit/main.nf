// autobin to get reference coverages
process cnvkit_autobin {

	label 'cnvkit'

	publishDir "${params.outdir}/CnvKit/Ref", mode: 'copy'

	input:
	path(bams)

	output:
	path(targets)

	script:
	targets = baits.getBaseName() + ".target.bed"

	"""
		cp ${params.cnv_annotation} . 
		gunzip -c refFlat.txt.gz > refFlat.txt
		cnvkit.py access $params.fasta -o access.bed -x ${params.cnv_blacklist} -x ${params.cnv_exclusion}
		cnvkit.py autobin -m wgs -b 50000 -g access.bed --annotate refFlat.txt *.bam
	"""
}

// get per sample coverage for targets and antitargets
process cnvkit_coverage {

	label 'cnvkit'

	publishDir "${params.outdir}/${indivID}/${sampleID}/CnvKit/Processing", mode: 'copy'

	input:
	tuple val(indivID), val(sampleID), path(bam),path(bai)
	path(targets)

	output:
	path(cnn)
	tuple val(indivID),val(sampleID),path(cnn)

	script:
	cnn = bam.getBaseName() + ".targetcoverage.cnn"

	"""
		cnvkit.py coverage $bam $targets -o $cnn -p ${task.cpus}
	"""

}

process cnvkit_reference {

	label 'cnvkit'

	publishDir "${params.outdir}/CnvKit", mode: 'copy'

	input:
	path('*')

	output:
	path(cnn)

	script:
	cnn = "cnvkit_ref_" + params.run_name + ".cnn"

	"""

		cnvkit.py reference *targetcoverage.cnn --fasta $params.fasta -o $cnn
	"""
}

// correct biases and stuff
process cnvkit_process {

	label 'cnvkit'
	input:
	tuple val(indivID),val(sampleID),path(cnn),path(cnn_anti)
	path(ref)

	output:
	tuple val(indivID),val(sampleID),path(cnr),path(cns)

	script:
	cns = cnn.getBaseName() + ".cns"
	cnr = cnn.getBaseName() + ".cnr"
	
	"""
		cnvkit.py fix $cnn $cnn_anti $ref -o $cnr
		cnvkit.py segment -m hmm-germline $cnr -o $cns -p ${task.cpus}
	"""

}

// Segmetrics
process cnvkit_segmetrics {

	label 'cnvkit'

	publishDir "${params.outdir}/${indivID}/${sampleID}/CnvKit/Processing", mode: 'copy'

	input:
	tuple val(indivID),val(sampleID),path(cnr),path(cns)

	output:
	tuple val(indivID),val(sampleID),path(cnr),path(seg_cns)

	script:
	seg_cns = cns.getBaseName() + ".segmetrics.cns"

	"""
		cnvkit.py segmetrics -s $cns $cnr --ci
	"""
}

// Attach confidence interfals
process cnvkit_call {

	label 'cnvkit'

	publishDir "${params.outdir}/${indivID}/${sampleID}/CnvKit/Processing", mode: 'copy'

	input:
	tuple val(indivID),val(sampleID),path(cnr),path(cns)

	output:
	tuple val(indivID),val(sampleID),path(cnr),path(call_cns)

	script:
	call_cns = cns.getBaseName() + ".call.cns"

	"""
		cnvkit.py call $cns --filter ci
	"""

}

// Metrics per gene
process cnvkit_genemetrics {

	label 'cnvkit'

	publishDir "${params.outdir}/${indivID}/${sampleID}/CnvKit/Metrics", mode: 'copy'

	input:
	tuple val(indivID), val(sampleID), path(cnr),path(cns)

	output:
	path(metrics)

	script:

	metrics = cnr.getBaseName() + ".genemetrics.txt"

	"""
		cnvkit.py genemetrics -s $cns $cnr -t 0.2 > $metrics
	"""

}

// find putative breaks
process cnvkit_breaks {

	label 'cnvkit'

	publishDir "${params.outdir}/${indivID}/${sampleID}/CnvKit/Metrics", mode: 'copy'

	input:
	tuple val(indivID),val(sampleID),path(cnr),path(cns)

	output:
	path(breaks)

	script:

	breaks = cnr.getBaseName() + ".breaks.txt"

	"""
		cnvkit.py breaks $cns $cnr > $breaks
	"""

}

// make useful output formats
process cnvkit_export {

	label 'cnvkit'

	publishDir "${params.outdir}/${indivID}/${sampleID}/CnvKit", mode: 'copy'

	input:
	tuple val(indivID), val(sampleID),path(cnr),path(call_cns)

	output:
	tuple path(bed),path(vcf)

	script:
	bed = call_cns.getBaseName() + ".bed"
	vcf = call_cns.getBaseName() + ".vcf"

	"""
		cnvkit.py export bed $call_cns -o $bed
		cnvkit.py export vcf $call_cns -i $sampleID -o $vcf
	"""

}

process cnvkit_plots {

	label 'cnvkit'

	publishDir "${params.outdir}/${indivID}/${sampleID}/CnvKit/Plots", mode: 'copy'

	input:
	tuple val(indivID),val(sampleID),path(cnr),path(call_cns)

	output:
	path(scatter)
	path(diagram)

	script:
	scatter = call_cns.getBaseName() + ".scatter.pdf"
	diagram = call_cns.getBaseName() + ".diagram.pdf"

	"""
		cnvkit.py scatter --y-min -4 --y-max 4 -o $scatter -s $call_cns $cnr
		cnvkit.py diagram -o $diagram -s $call_cns $cnr
	"""
}

