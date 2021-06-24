// Import generic module functions
include { saveFiles } from '../process/functions'

params.options = [:]

/*
 * Convert library files to respect input requirements of MAGECK, PINAPL-PY, CRISPR-ANALYZER
 */
process CONVERT_LIBRARY_FILE {
    tag "$library"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'library', publish_id:'') }

    conda     (params.enable_conda ? "conda-forge::python=3.8.3 pandas=1.2.4" : null)
    container "quay.io/biocontainers/python:3.8.3"

    input:
    path library
    
    output:
    path 'pinapl_library.csv', emit:pinapl_library
    path 'mageck_library.csv', emit:mageck_library


    script:  // This script is bundled with the pipeline, in nf-core/crisprquant/bin/
    """
    file_conversions.py $library pinapl_library.csv mageck_library.csv
    """
}