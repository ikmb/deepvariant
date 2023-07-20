#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/**
===============================
WGS Pipeline
===============================

This Pipeline performs variant calling
using Google DeepVariant

### Homepage / git
git@github.com:ikmb/deepvariant.git
### Implementation
Re-Implemented in Q3 2023

Author: Marc P. Hoeppner, m.hoeppner@ikmb.uni-kiel.de

**/

def summary = [:]

WorkflowMain.initialise(workflow, params, log)
WorkflowDeepvariant.initialise( params, log)

//
// Summary of all options
//
summary['runName'] = params.run_name
summary['Samples'] = params.samples
summary['Current home'] = "$HOME"
summary['Current user'] = "$USER"
summary['Current path'] = "$PWD"
summary['Assembly'] = params.assembly
summary['CommandLine'] = workflow.commandLine

if (workflow.containerEngine) {
        summary['Container'] = "$workflow.containerEngine - $workflow.container"
}
summary['SessionID'] = workflow.sessionId

def multiqc_report = Channel.from([])

include { DEEPVARIANT_PIPELINE } from "./workflows/deepvariant_pipeline.nf"

workflow {

    DEEPVARIANT_PIPELINE()
    multiqc_report = multiqc_report.mix(DEEPVARIANT_PIPELINE.out.qc).toList()

}

workflow.onComplete {

        log.info "========================================="
        log.info "Duration:		$workflow.duration"
        log.info "========================================="

}

