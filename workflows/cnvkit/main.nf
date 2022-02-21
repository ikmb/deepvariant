include { cnvkit_autobin ;  cnvkit_coverage ; cnvkit_reference ; cnvkit_process ; cnvkit_segmetrics ; cnvkit_call ; cnvkit_genemetrics ;  cnvkit_breaks ; cnvkit_export ; cnvkit_plots } from "./../../modules/cnvkit/main.nf" params(params)


workflow CNVKIT {

	take:
		bed
		bam

	main:
		cnvkit_autobin(bed,bam.map { i,s,b,n -> [ b,n] }.collect())
		cnvkit_coverage(bam,cnvkit_autobin.out.collect())
		cnvkit_reference(cnvkit_coverage.out[0].collect())	
		cnvkit_process(cnvkit_coverage.out[1],cnvkit_reference.out.collect())
		cnvkit_segmetrics(cnvkit_process.out)
		cnvkit_call(cnvkit_segmetrics.out)
		cnvkit_genemetrics(cnvkit_call.out)
		cnvkit_breaks(cnvkit_call.out)
		cnvkit_export(cnvkit_call.out)
		cnvkit_plots(cnvkit_call.out)

	emit:
		bed = cnvkit_export.out[0]
		vcf = cnvkit_export.out[1]				

}
