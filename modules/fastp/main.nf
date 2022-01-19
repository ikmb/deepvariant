process trim {

	label 'fastp'

        input:
        tuple val(indivID), val(sampleID), val(libraryID), val(rgID), val(platform_unit), val(platform), val(platform_model), val(center), val(run_date), path(fastqR1), path(fastqR2)

        output:
        tuple val(indivID), val(sampleID), val(libraryID), val(rgID), val(platform_unit), val(platform), val(platform_model), val(center), val(run_date),file(left),file(right)
        tuple path(json),path(html)

        script:
        left = file(fastqR1).getBaseName() + ".trimmed.fastq.gz"
        right = file(fastqR2).getBaseName() + ".trimmed.fastq.gz"
        json = indivID + "_" + sampleID + "_" + libraryID + "-" + rgID + ".fastp.json"
        html = indivID + "_" + sampleID + "_" + libraryID + "-" + rgID + ".fastp.html"

        """
                fastp --in1 $fastqR1 --in2 $fastqR2 --out1 $left --out2 $right --detect_adapter_for_pe -w ${task.cpus} -j $json -h $html
        """

}

