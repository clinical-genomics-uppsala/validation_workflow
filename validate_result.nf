
nextflow.enable.dsl=2
  
params.input = "test_input.tsv"
params.publish_dir = './results'


process validate_vcf {
    input:
      tuple val(input_file), val(checksum)

    output:
      tuple val(input_file), val(checksum), env(md5) 

    when: 
      input_file =~ /vcf$/
      
    """
    md5=\$(cat ${input_file} | 
           awk 'BEING{START_PRINT=0}{if(START_PRINT) print(\$0); if(/^#CHROM/) START_PRINT=1; }'|  
           awk '{if(\$8 ~/&/) {split(\$8, arr, "[&|]"); joined=""; asort(arr, arr_s);for (i=1; i <= length(arr_s); i++){ joined=joined"&"arr_s[i];} \$8=joined}; print(\$0)}' |  
           md5sum |
           awk '{print(\$1)}')
    """
}

process validate_vcf_gz {
    input:
      tuple val(input_file), val(checksum)

    output:
      tuple val(input_file), val(checksum), env(md5)

    when: 
      input_file =~ /vcf.gz$/
      
    """
    md5=\$(zcat ${input_file} | 
           awk 'BEING{START_PRINT=0}{if(START_PRINT) print(\$0); if(/^#CHROM/) START_PRINT=1; }'|  
           awk '{if(\$8 ~/&/) {split(\$8, arr, "[&|]"); joined=""; asort(arr, arr_s);for (i=1; i <= length(arr_s); i++){ joined=joined"&"arr_s[i];} \$8=joined}; print(\$0)}' |  
           md5sum |
           awk '{print(\$1)}')
    """
}

process validate_mutations_and_coverage {
    input:
      tuple val(input_file), val(checksum)

    output:
      tuple val(input_file), val(checksum), env(md5)

    when:
      input_file =~ /coverage_and_mutations\.tsv$/

    """
    md5=\$(cat ${input_file} |
           awk '{if(\$5 ~/&/) {split(\$5, arr, "&"); joined=""; asort(arr, arr_s);for (i=1; i <= length(arr_s); i++){ joined=joined"&"arr_s[i];} \$5=joined};' | 
           md5sum |
           awk '{print(\$1)}')
    """
}

process validate_metrics {
    input:
      tuple val(input_file), val(checksum)

    output:
      tuple val(input_file), val(checksum), env(md5)

    when:
      input_file =~ /insert_size_metrics\.txt$|HsMetrics\.txt$|alignment_summary_metrics\.txt$|duplication_metrics\.txt$/

    """
    md5=\$(cat ${input_file} |
           awk 'BEING{START_PRINT=0}{if(/## METRICS/) START_PRINT=1; if(START_PRINT) print(\$0);}' |
           md5sum | 
           awk '{print(\$1)}')
    """

}

process validate_multiqc {
    input:
      tuple val(input_file), val(checksum)

    output:
      tuple val(input_file), val(checksum), env(md5)

    when:
        input_file =~ /multiqc_.*\.html$/


    """
    md5=\$(cat ${input_file} |
        sed 's/generated on [0-9:, -]*//' | sed 's/mqc_analysis_path.*code/mqc_analysis_pathcode/g' | sed 's/able[_ ][A-Za-z]*/able_/g' |
        md5sum |
        awk '{print(\$1)}')
    """
}

process validate_samtool_stats {
    input:
      tuple val(input_file), val(checksum)

    output:
      tuple val(input_file), val(checksum), env(md5)

    when:
      input_file =~ /samtools-stats\.txt$/

    """
    md5=\$(cat ${input_file} |
           awk 'BEING{START_PRINT=0}{if(START_PRINT) print(\$0); if(/command line was/) START_PRINT=1; }' |
           md5sum |
           awk '{print(\$1)}')
    """

}

process validate_collection_of_files {
    input:
      tuple val(input_file), val(checksum)

    output:
      tuple val(input_file), val(checksum), env(md5)

    when:
        input_file =~ /msisensor_pro\.score\.tsv$|cnvkit\.loh\.cns$|cnvkit\.scatter\.png$|gatk_cnv\.seg$|cnv_report\.tsv$|\.gene_fuse_report\.tsv$|hrd_score\.txt$|TMB\.txt$|hrd.*\.txt$|table$|cnv\.html$/
 
    """
    md5=\$(cat ${input_file} |
           md5sum |
           awk '{print(\$1)}')
    """
}

process validate_genefuse {
    input:
      tuple val(input_file), val(checksum)

    output:
      tuple val(input_file), val(checksum), env(md5)

    when:
        input_file =~ /gene_fuse_fusions\.txt$/
 
    """
    md5=\$(cat ${input_file} | 
         awk '{if(/^# genefuse/) exit(0); print(\$0)}' |
         md5sum |
         awk '{print(\$1)}')
    """
}

process validate_bam {
    input:
      tuple val(input_file), val(checksum)

    output:
       tuple val(input_file), val(checksum), env(md5)

    when:
        input_file =~ /\.bam$/
 
    """
    md5=\$(samtools view ${input_file} |
           md5sum | awk '{print(\$1)}'); 
    echo \$md5 | awk 'BEGIN{ORS=""}{if(\$1 == "${checksum}") {print("Validated: ${input_file}"); exit(0)} else { print("Failed validation: ${input_file}: ${checksum} != "\$1); exit(1)}}'
    """
}

process checksum_bam {
    input:
      tuple val(input_file), val(checksum)

    output:
       tuple val(input_file), val(checksum), env(md5)

    when:
        input_file =~ /\.bam$/
 
    """
    md5=\$(samtools view ${input_file} |
           md5sum |
           awk '{print(\$1)}')
    """
}

process validate_checksum {
    input:
        tuple val(input_file), val(checksum), val(calculated_checksum)
   
    output:
        stdout

    """
    echo ${calculated_checksum} |
       awk 'BEGIN{ORS=""}{if(\$1 == "${checksum}") {print("Validated: ${input_file}"); exit(0)} else { print("Failed validation: ${input_file}: ${checksum} != "\$1); exit(1)}}' 
    """
}

process create_checksum_file {
    input:
        tuple val(input_file), val(old_md5), val(new_md5)
    
    output:
       path 'new_test_validation.tsv'

    """
    echo "${input_file}\t${new_md5}" > new_test_validation.tsv 
    """
}



workflow validate {
      ch = Channel.fromPath(params.input).splitCsv( header:true, sep:"\t" ).map(row -> [file(row.file), row.checksum])
      ch | (
            validate_vcf & 
            validate_vcf_gz &
            validate_metrics &
            validate_multiqc &
            validate_mutations_and_coverage &
            validate_samtool_stats &
            validate_collection_of_files &
            validate_genefuse &
            validate_bam ) | mix | validate_checksum
}

workflow create_validation_data {
    ch = Channel.fromPath(params.input).splitCsv( header:true, sep:"\t" ).map(row -> [file(row.file), row.checksum])
    ch | (
            validate_vcf &  
            validate_vcf_gz &
            validate_metrics &
            validate_multiqc &
            validate_mutations_and_coverage &
            validate_samtool_stats &
            validate_collection_of_files &
            validate_genefuse &
            checksum_bam ) | mix |
            create_checksum_file | collectFile(name: "test.txt", newLine: false, storeDir: "result") | view
}

workflow {
   main:
      validate()
      //create_validation_data()

}

