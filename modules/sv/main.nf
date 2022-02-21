process pbsv_sig {

	label 'pbsv'

	input:
	tuple val(indivID),val(sampleID),path(bam),path(pbi)
	path(repeat_ref)

	output:
	tuple val(indivID),val(sampleID),path(sig)
	path(sig)

	script:
	sig = bam.getBaseName() + ".svsig.gz"

	"""
		pbsv discover --ccs --tandem-repeats $repeat_ref $bam $sig
	"""
}

process pbsv_call {

	publishDir "${params.outdir}/SVs", mode: 'copy'

	label 'pbsv'

	input:
	path(sigs)

	output:
	path(vcf)

	script:
	vcf = "SVs.vcf"

	"""
		pbsv call --ccs -j ${task.cpus} ${params.fasta} ${sigs} $vcf
	"""

}

process manta {

	label 'manta'

	publishDir "${params.outdir}/${indivID}/${sampleID}/SV", mode: 'copy'

	input:
	tuple val(indivID),val(sampleID),path(bam),path(bai)
	path(bed)

	output:
	tuple val(indivID),val(sampleID),path(sv),file(sv_tbi)
	tuple val(indivID),val(sampleID),path(indel),file(indel_tbi)
	tuple val(indivID),val(sampleID),path(sv_can),file(sv_can_tbi)

	script:
	
	sv = indivID + "_" + sampleID + ".diploidSV.vcf.gz"
        sv_tbi = sv + ".tbi"
        indel = indivID + "_" + sampleID + ".candidateSmallIndels.vcf.gz"
        indel_tbi = indel + ".tbi"
        sv_can = indivID + "_" + sampleID + ".candidateSV.vcf.gz"
        sv_can_tbi = sv_can + ".tbi"
        
	"""
        	configManta.py --bam $bam --referenceFasta ${FASTA} --runDir manta --callRegions $bed_gz --exome

                manta/runWorkflow.py -j ${task.cpus}

                cp manta/results/variants/diploidSV.vcf.gz $sv
                cp manta/results/variants/diploidSV.vcf.gz.tbi $sv_tbi
                cp manta/results/variants/candidateSmallIndels.vcf.gz $indel
                cp manta/results/variants/candidateSmallIndels.vcf.gz.tbi $indel_tbi
                cp manta/results/variants/candidateSV.vcf.gz $sv_can
                cp manta/results/variants/candidateSV.vcf.gz.tbi $sv_can_tbi

	"""

}

