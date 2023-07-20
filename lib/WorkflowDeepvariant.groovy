//
// This file holds several functions specific to the workflow/esga.nf in the nf-core/esga pipeline
//

class WorkflowDeepvariant {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        
        genomeExistsError(params, log)

        if (!params.run_name) {
            log.info  "Must provide a run_name (--run_name)"
            System.exit(1)
        }
    
        if (!params.tools) {
            log.info "No analysis tools specified, performing only alignments and QC!"
        }

    }

    private static void genomeExistsError(params, log) {
        if (params.assembly && !params.genomes.containsKey(params.assembly)) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome assembly '${params.assembly}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available assemblies  are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }
    }

}
