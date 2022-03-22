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

process picard_wgs_metrics {

	label 'picard'

	input:
	tuple val(indivID),val(sampleID),path(bam),path(bai)
        path(bed)

	output:
	path(picard_stats)

	script:
	base_name = bam.getBaseName()
	picard_stats = base_name + "_wgs_metrics.txt"
	intervals = bed.getBaseName() + ".interval_list"

	"""
		picard BedToIntervalList I=$bed O=$intervals SD=${params.dict}
		picard CollectWgsMetrics I=bam REFERENCE_SEQUENCE=${params.fasta} O=$picard_stats INTERVALS=$intervals
	"""
}

process multiqc {

        label 'multiqc'

        publishDir "${params.outdir}/MultiQC", mode: 'copy'

        input:
        path('*')

        output:
        path("multiqc_report.html")

        script:

        """
		cp $baseDir/conf/multiqc_config.yaml multiqc_config.yaml
                multiqc .
        """
}

