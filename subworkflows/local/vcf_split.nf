
include { VCF_SPLIT } from '../../modules/local/split_vcf/main.nf'

workflow SPLIT_VCF {

    take:
    combined_vcf

    main:
        VCF_SPLIT(combined_vcf)
        split_vcf  = VCF_SPLIT.out.vcf.flatMap { patient, vcfs ->
            vcfs.collect { vcf ->
            def parts = vcf.toString().tokenize('/')[-1].tokenize('_')
            def sample = parts[2]
           tuple(patient, sample, vcf)
           }
        }

    emit:
    split_vcf                                  // channel: [ val(patient), sample, vcf ]
}
