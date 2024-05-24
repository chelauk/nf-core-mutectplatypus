//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        | splitCsv ( header:true )
        | map { row ->
                meta = row.subMap('patient','sample','id','status')
                [ meta, [ file(row.bam, checkIfExists:true), file(row.bai, checkIfExists: true)]]
                 }
        .set { bams }

    emit:
    bams                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}