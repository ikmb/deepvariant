process PBSV_SIG {

	tag "${meta.patient_id}|${meta.sample_id}"

	label 'pbsv'

	input:
	tuple val(meta),path(bam),path(pbi)
	path(repeat_ref)

	output:
	tuple val(meta),path(sig), emit: sig
	path(sig)

	script:
	sig = bam.getBaseName() + ".svsig.gz"

	"""
		pbsv discover --ccs --tandem-repeats $repeat_ref $bam $sig
	"""
}

process PBSV_CALL {

	publishDir "${params.outdir}/SVs", mode: 'copy'

	label 'pbsv'

	input:
	path(sigs)

	output:
	path(vcf), emit: vcf

	script:
	vcf = "SVs.vcf"

	"""
		pbsv call --ccs -j ${task.cpus} ${params.fasta} ${sigs} $vcf
	"""

}

process MANTA {

	tag "${meta.patient_id}|${meta.sample_id}"

	label 'manta'

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/SV", mode: 'copy'

	input:
	tuple val(meta),path(bam),path(bai)
	tuple path(bed_gz), path(bed_gz_tbi)

	output:
	tuple val(meta),path(sv),file(sv_tbi), emit: sv
	tuple val(meta),path(indel),file(indel_tbi), emit: indel
	tuple val(meta),path(sv_can),file(sv_can_tbi), emit: sv_can

	script:
	
	sv = meta.patient_id + "_" + meta.sample_id + ".diploidSV.vcf.gz"
        sv_tbi = sv + ".tbi"
        indel = meta.patient_id + "_" + meta.sample_id + ".candidateSmallIndels.vcf.gz"
        indel_tbi = indel + ".tbi"
        sv_can = meta.patient_id + "_" + meta.sample_id + ".candidateSV.vcf.gz"
        sv_can_tbi = sv_can + ".tbi"
        
	"""
        	configManta.py --bam $bam --referenceFasta ${params.fasta} --runDir manta --callRegions $bed_gz --exome

                manta/runWorkflow.py -j ${task.cpus}

                cp manta/results/variants/diploidSV.vcf.gz $sv
                cp manta/results/variants/diploidSV.vcf.gz.tbi $sv_tbi
                cp manta/results/variants/candidateSmallIndels.vcf.gz $indel
                cp manta/results/variants/candidateSmallIndels.vcf.gz.tbi $indel_tbi
                cp manta/results/variants/candidateSV.vcf.gz $sv_can
                cp manta/results/variants/candidateSV.vcf.gz.tbi $sv_can_tbi

	"""

}

