process align {

        input:
        tuple val(indivID), val(sampleID), val(libraryID), val(rgID), val(platform_unit), val(platform), val(platform_model), val(center), val(run_date),path(fastqR1),path(fastqR2)

        output:
        tuple val(indivID), val(sampleID), file(bam)
        val(sample_name)

        script:
        bam = indivID + "_" + sampleID + "." + libraryID + "_" + rgID + ".aligned.cram"
        sample_name = "${indivID}_${sampleID}"

        """
                bwa-mem2 mem -K 1000000 -H $params.dict -M -R "@RG\\tID:${rgID}\\tPL:ILLUMINA\\tPU:${platform_unit}\\tSM:${indivID}_${sampleID}\\tLB:${libraryID}\\tDS:${params.fasta}\\tCN:${center}" \
                -t ${task.cpus} ${params.bwa2_index} $fastqR1 $fastqR2 \
                | samtools fixmate -@ ${task.cpus} -m - - \
                | samtools sort -m 4G --reference $params.fasta -O CRAM -@ ${task.cpus} -o $bam -
        """

}
