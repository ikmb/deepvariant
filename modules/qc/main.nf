process wgs_coverage {

	label 'mosdepth'

	input:
	tuple val(indivID),val(sampleID),path(bam),path(bai)
	path(bed)

        output:
        tuple path(genome_bed_coverage),path(genome_global_coverage)

        script:
        base_name = bam.getBaseName()
        genome_bed_coverage = base_name + ".mosdepth.region.dist.txt"
        genome_global_coverage = base_name + ".mosdepth.global.dist.txt"

        """
                mosdepth -t ${task.cpus} -n -f $params.fasta -x -Q 10 -b $bed $base_name $bam
        """
}

process multiqc {

        publishDir "${params.outdir}/MultiQC", mode: 'copy'

        input:
        path('*')

        output:
        path("multiqc_report.html")

        script:

        """
                multiqc .
        """
}

