// Import generic module functions
include { saveFiles } from '../process/functions'

params.options = [:]

/*
 * Create Pinapl config file
 */
process GET_LIBRARY_FASTA {
    tag "$library"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'library', publish_id:'') }

    conda     (params.enable_conda ? "conda-forge::python=3.8.3 anaconda::yaml=0.25" : null)
    container "quay.io/biocontainers/python:3.8.3"

    input:
    path library

    output:
    path '*.fasta'


    script:  // This script is bundled with the pipeline, in nf-core/crisprquant/bin/
    """
    create_pinapl_config.py --LibFilename $pinapl_library --seq_5_end 
    """
}