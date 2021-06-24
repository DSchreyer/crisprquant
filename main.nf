#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/crisprquant
========================================================================================
 nf-core/crisprquant Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/crisprquant
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////

def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    // TODO nf-core: Update typical command used to run pipeline
    def command = "nextflow run nf-core/crisprquant --input samplesheet.csv --library -profile docker"
    log.info Schema.params_help(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --        GENOME PARAMETER VALUES           -- */
////////////////////////////////////////////////////

params.fasta = Checks.get_genome_attribute(params, 'fasta')

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////

def summary_params = Schema.params_summary_map(workflow, params, json_schema)
log.info Schema.params_summary_log(workflow, params, json_schema)

////////////////////////////////////////////////////
/* --          PARAMETER CHECKS                -- */
////////////////////////////////////////////////////

// Check that conda channels are set-up correctly
if (params.enable_conda) {
    Checks.check_conda_channels(log)
}

// Check AWS batch settings
Checks.aws_batch(workflow, params)

// Check the hostnames against configured profiles
Checks.hostname(workflow, params, log)

// Check genome key exists if provided
Checks.genome_exists(params, log)

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.library) { ch_library = file(params.library) } else { exit 1, 'Library file not specified!' }
//if (params.fasta) { ch_fasta = file(params.fasta) } else { exit 1, 'Genome fasta file not specified!' }

if (params.pinapl_config) { ch_pinaplconfig = file(params.pinapl_config) } else { exit 1, 'pinAPL-py config file does not exist' }
if (params.pinapl_datasheet) { ch_datasheet = file(params.pinapl_datasheet) } else { exit 1, 'pinAPL-py datasheet does not exist' }
////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

////////////////////////////////////////////////////
/* --       IMPORT MODULES / SUBWORKFLOWS      -- */
////////////////////////////////////////////////////

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? " --title \"$params.multiqc_title\"" : ''

// Local: Modules
include { GET_SOFTWARE_VERSIONS } from './modules/local/process/get_software_versions' addParams( options: [publish_files : ['csv':'']] )

// Local: Sub-workflows
include { INPUT_CHECK           } from './modules/local/subworkflow/input_check'       addParams( options: [:]                          )
include { GET_LIBRARY_FASTA     } from './modules/local/subworkflow/get_library_fasta'         addParams( options: [:]                          )
include { CONVERT_LIBRARY_FILE  } from './modules/local/subworkflow/convert_library_file'   addParams( options: [:]                             )
include { MAGECK_COUNT          } from './modules/local/subworkflow/mageck_count'   addParams( options: [:]                             )
include { PINAPLPY              } from './modules/local/pinaplpy'   addParams( options: [:]                             )

// nf-core/modules: Modules
include { FASTQC                } from './modules/nf-core/software/fastqc/main'        addParams( options: modules['fastqc']            )
include { MULTIQC               } from './modules/nf-core/software/multiqc/main'       addParams( options: multiqc_options              )
include { BOWTIE2_BUILD         } from './modules/nf-core/software/bowtie2/build/main' addParams( options: [:]                          )
include { BOWTIE2_ALIGN         } from './modules/nf-core/software/bowtie2/align/main' addParams( options: [:]                          )
include { CUTADAPT              } from './modules/nf-core/software/cutadapt/main'      addParams( options: [:]                          )
// include { TRIMGALORE            } from './modules/nf-core/software/trimgalore/main'      addParams( options: [:]                          )


////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

// Info required for completion email and summary
def multiqc_report = []

workflow {

    ch_software_versions = Channel.empty()

    /*
     * SUBWORKFLOW: Read in samplesheet, validate and stage input files
     */
    INPUT_CHECK (
        ch_input
    )
    CONVERT_LIBRARY_FILE (
        ch_library
    )

    GET_LIBRARY_FASTA (
        CONVERT_LIBRARY_FILE.out.mageck_library
    )

    // TODO: GUNZIP FASTQ FILES

    /*
     * MODULE: Run FastQC
     */
    //INPUT_CHECK.out.reads.view()

    // INPUT_CHECK.out.reads.buffer { it == "*.fastq.gz" }.view()

    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))

    BOWTIE2_BUILD (
        GET_LIBRARY_FASTA.out
    )

    CUTADAPT (
        INPUT_CHECK.out.reads
    )

    BOWTIE2_ALIGN (
        CUTADAPT.out.reads, BOWTIE2_BUILD.out.index
    )

    // generate channel for mageck 
//    BOWTIE2_ALIGN.out.bam.map{
//    meta, file ->
//    ["group", meta, file]}.groupTuple().map{
//        group, meta, file ->
//        [meta, file]
//    }.set{aligned_reads}

    CUTADAPT.out.reads.map {
    meta, file ->
    ["group", meta, file]}.groupTuple().map{
        group, meta, file ->
        [meta, file]
    }.set{aligned_reads}

    MAGECK_COUNT (
        aligned_reads, CONVERT_LIBRARY_FILE.out.mageck_library
    )

    PINAPLPY (
        ch_pinaplconfig, ch_datasheet, CONVERT_LIBRARY_FILE.out.pinapl_library, aligned_reads
    )

    /*
     * MODULE: Pipeline reporting
     */
    GET_SOFTWARE_VERSIONS ( 
        ch_software_versions.map { it }.collect()
    )

    /*
     * MultiQC
     */  
    if (!params.skip_multiqc) {
        workflow_summary    = Schema.params_summary_multiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
        ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
        
        MULTIQC (
            ch_multiqc_files.collect()
        )
        multiqc_report       = MULTIQC.out.report.toList()
        ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
    }
}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
    Completion.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    Completion.summary(workflow, params, log)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
