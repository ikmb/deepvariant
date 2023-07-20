//
// This file holds several functions specific to the workflow/esga.nf in the nf-core/esga pipeline
//

class WorkflowMain {

    //
    // Check and validate parameters
    //
    //
    // Validate parameters and print summary to screen
    //
    public static void initialise(workflow, params, log) {
        // Print help to screen if required
        
        log.info header(workflow)

        if (params.help) {
            log.info help(workflow, params, log)
            System.exit(0)
        }

    }

    public static String header(workflow) {
        def headr = ''
        def info_line = "IKMB Research WGS pipeline | version ${workflow.manifest.version}"
        headr = """
    ===============================================================================
    ${info_line}
    ===============================================================================
    """
        return headr
    }

    public static String help(workflow, params, log) {
        def command = "nextflow run ${workflow.manifest.name} --samples Samples.csv --assembly GRCh38 --kit xGen_v2 -profile diagnostic"
        def help_string = ''
        // Help message
        help_string = """

           Usage: nextflow run ikmb/exome-seq --assembly GRCh38 --kit xGen_v2 --samples Samples.csv
           This example will perform an exome analysis against the ALT-free hg38 assembly, assuming that exome reads were generated with
           the IDT xGen v2 kit and using DeepVariant with GLNexus.

           Required parameters:
           --samples                      A sample list in CSV format (see website for formatting hints)
           --assembly                     Name of the reference assembly to use
           --tools                        Comma-separated list of tools to run. 
           Optional parameters:
           --run_name                     A descriptive name for this pipeline run
           Output:
           --outdir                       Local directory to which all output is written (default: results)
        """
        return help_string
    }

}
