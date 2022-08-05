include { CNVKIT_AUTOBIN ;  CNVKIT_COVERAGE ; CNVKIT_REFERENCE ; CNVKIT_PROCESS ; CNVKIT_SEGMETRICS ; CNVKIT_CALL ; CNVKIT_GENEMETRICS ;  CNVKIT_BREAKS ; CNVKIT_EXPORT ; CNVKIT_PLOTS } from "./../../modules/cnvkit/main.nf"


workflow CNVKIT {

	take:
		bed
		bam

	main:
		CNVKIT_AUTOBIN(bed,bam.map { i,s,b,n -> [ b,n] }.collect())
		CNVKIT_COVERAGE(bam,cnvkit_autobin.out.collect())
		CNVKIT_REFERENCE(cnvkit_coverage.out[0].collect())	
		CNVKIT_PROCESS(cnvkit_coverage.out[1],cnvkit_reference.out.collect())
		CNVKIT_SEGMETRICS(cnvkit_process.out)
		CNVKIT_CALL(cnvkit_segmetrics.out)
		CNVKIT_GENEMETRICS(cnvkit_call.out)
		CNVKIT_BREAKS(cnvkit_call.out)
		CNVKIT_EXPORT(cnvkit_call.out)
		CNVKIT_PLOTS(cnvkit_call.out)

	emit:
		bed = CNVKIT_EXPORT.out[0]
		vcf = CNVKIT_EXPORT.out[1]				

}
