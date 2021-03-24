// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

/*
 * Reformat design file and check validity
 */
process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', publish_id:'') }

    conda     (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "quay.io/biocontainers/python:3.8.3"

    input:
    path samplesheet
    
    output:
    path '*.csv'


    script:  // This script is bundled with the pipeline, in nf-core/crisprquant/bin/
    """
    check_samplesheet.py $samplesheet samplesheet.valid.csv
    """
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def get_samplesheet_paths(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample
    meta.single_end   = row.single_end.toBoolean()

    def array = []
    if (!file(row.fastq).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq}"
    }
    if (meta.single_end) {
        array = [ meta, file(row.fastq) ]
    }
    return array
}