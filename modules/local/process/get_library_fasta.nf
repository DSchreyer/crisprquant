// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

/*
 * Generate Fasta file from Library txt file
 */
process GET_LIBRARY_FASTA {
    tag "$library"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'library', publish_id:'') }

    conda     (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "quay.io/biocontainers/python:3.8.3"

    input:
    path library
    
    output:
    path '*.fasta'


    script:  // This script is bundled with the pipeline, in nf-core/crisprquant/bin/
    """
    get_library_fasta.py $library library.fasta
    """
}