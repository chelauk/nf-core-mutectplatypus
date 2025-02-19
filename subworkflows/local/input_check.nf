//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
        SAMPLESHEET_CHECK ( samplesheet )
        SAMPLESHEET_CHECK.out.csv
            | splitCsv(header: true)
            | map { row ->
                meta = row.subMap('patient', 'sample', 'id', 'status')
                if (row.containsKey('sequenza_control')) {
                    meta['sequenza_control'] = row['sequenza_control'] // Only add if present
                }
                meta.id = meta.id ?: "${meta.patient}_${meta.sample}"
                [meta, [file(row.bam, checkIfExists: true), file(row.bai, checkIfExists: true)]]
            }
            | set { bams }

    emit:
        bams                                     // channel: [ val(meta), [ reads ] ]
        versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}
