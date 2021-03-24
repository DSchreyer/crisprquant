/*
 * Check input samplesheet and get read channels
 */

params.options = [:]

include {
    FASTA_REF
    } from '../process/get_library_fasta' addParams( options: params.options )

workflow INPUT_CHECK {
    take:
    library // file: /path/to/library.txt

    main:
    FASTA_REF ( library )
        .splitCsv ( header:true, sep:',' )
        .map { get_samplesheet_paths(it) }
        .set { reads }

    emit:
    fasta // channel: [ val(meta), [ reads ] ]
}
