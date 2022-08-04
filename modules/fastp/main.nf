process FASTP {

	tag "${meta.patient_id}|${meta.sample_id}|${meta.readgroup_id}"

	label 'fastp'

        input:
        tuple val(meta), path(fastqR1), path(fastqR2)

        output:
        tuple val(meta),file(left),file(right), emit: reads
        tuple path(json),path(html), emit: qc

        script:
        left = file(fastqR1).getBaseName() + ".trimmed.fastq.gz"
        right = file(fastqR2).getBaseName() + ".trimmed.fastq.gz"
        json = meta.patient_id + "_" + meta.sample_id + "_" + meta.library_id + "-" + meta.readgroup_id + ".fastp.json"
	html = meta.patient_id + "_" + meta.sample_id + "_" + meta.library_id + "-" + meta.readgroup_id + ".fastp.html"

        """
                fastp --in1 $fastqR1 --in2 $fastqR2 --out1 $left --out2 $right --detect_adapter_for_pe -w ${task.cpus} -j $json -h $html
        """

}

