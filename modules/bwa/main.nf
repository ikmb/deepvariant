process BWA {

	tag "${meta.patient_id}|${meta.sample_id}|${meta.readgroup_id}"

        input:
        tuple val(meta),path(fastqR1),path(fastqR2)

        output:
        tuple val(meta), file(bam), emit: bam
        val(meta), emit: sname

        script:
        bam = meta.patient_id + "_" + meta.sample_id + "." + meta.library_id + "_" + meta.readgroup_id + ".aligned.cram"

        """
                bwa-mem2 mem -K 1000000 -H $params.dict -M -R "@RG\\tID:${meta.readgroup_id}\\tPL:ILLUMINA\\tPU:${meta.platform_unit}\\tSM:${meta.patient_id}_${meta.sample_id}\\tLB:${meta.library_id}\\tDS:${params.fasta}\\tCN:CCGA" \
                -t ${task.cpus} ${params.bwa2_index} $fastqR1 $fastqR2 \
                | samtools fixmate -@ ${task.cpus} -m - - \
                | samtools sort -m 4G --reference $params.fasta -O CRAM -@ ${task.cpus} -o $bam -
        """

}
