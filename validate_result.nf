
nextflow.enable.dsl=2
  
params.input = "test_input.tsv"
params.publish_dir = './results'


process validate_vcf {
    input:
      tuple val(input_file), val(checksum)

    output:
      tuple val(input_file), val(checksum), env(md5) 

    when: 
      input_file =~ /vcf$/ && !(input_file =~ /svdb_query.vcf$/) && !(input_file =~ /id_snps.vcf$/) && !(input_file =~ /jumble.vcf$/)
      
    """
    md5=\$(cat ${input_file} | 
           awk 'BEING{START_PRINT=0}{if(START_PRINT) print(\$0); if(/^#CHROM/) START_PRINT=1; }'|  
           awk '{if(\$8 ~/&/) {split(\$8, arr, "[&|]"); joined=""; asort(arr, arr_s);for (i=1; i <= length(arr_s); i++){ joined=joined"&"arr_s[i];} \$8=joined}; print(\$0)}' |  
           md5sum |
           awk '{print(\$1)}')
    """
}

process validate_svdb_query_vcf {
    input:
      tuple val(input_file), val(checksum)

    output:
      tuple val(input_file), val(checksum), env(md5)

    when: 
      input_file =~ /svdb_query.vcf$/
      
    """
    md5=\$(cat ${input_file} | 
           awk 'BEING{START_PRINT=0}{if(START_PRINT) print(\$0); if(/^#CHROM/) START_PRINT=1; }'|  
           awk '{print(\$1,\$2,\$5)}' |  
           md5sum |
           awk '{print(\$1)}')
    """
}

process validate_id_snps_vcf {
    input:
      tuple val(input_file), val(checksum)

    output:
      tuple val(input_file), val(checksum), env(md5)

    when: 
      input_file =~ /id_snps.vcf$/
      
    """
    md5=\$(cat ${input_file} | 
           awk '\$0 !~ /^#/ && NF' |  
           awk '{if(\$10 ~/:/) {split(\$10,arr,":")}; new=arr[1]; \$10=new; print(\$1,\$2,\$4,\$5,\$10)}' |  
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
      input_file =~ /vcf.gz$/ && !(input_file =~ /germline.vcf.gz$/)
      
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
      input_file =~ /coverage_and_mutations\./

    """
    md5=\$(cat ${input_file} |
           awk '{if(\$5 ~/&/) {split(\$5, arr, "&"); joined=""; asort(arr, arr_s);for (i=1; i <= length(arr_s); i++){ joined=joined"&"arr_s[i];} \$5=joined}; print(\$0)}' | 
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
        sed 's/generated on [0-9:, -]*//' | sed 's/mqc_analysis_path.*code/mqc_analysis_pathcode/g' | sed 's/able[A-Za-zw_ ]*/able/g' |
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
        input_file =~ /tc\.txt$|deletions\.tsv$|report\.[ct]sv$|cnvkit\.loh\.cns$|cnvkit\.scatter\.png$|gatk_cnv\.seg$|cnv_report\.tsv$|\.gene_fuse_report\.tsv$|TMB\.txt$|\.table$|\.purity.txt$|amplifications\.tsv$|fuseq_wes\.report\.csv$|fuseq_wes\.unfiltered\.results\.csv$|exon_skipping\.tsv$|fusion_report\.tsv$|ouse_keeping_gene_coverage\.tsv$|arriba.fusions\.tsv$|fusions\.txt$|fusion_predictions\.txt$|purecn_purity_ploidity\.csv$|sample_mixup_check\.tsv$/
 
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
         sed 's/^>[0-9]*,//g'  |
         sort |
         md5sum |
         awk '{print(\$1)}')
    """
}


process validate_dna_bam {
    input:
      tuple val(input_file), val(checksum)

    output:
       tuple val(input_file), val(checksum), env(md5)

    when:
        input_file =~ /[TN]\.bam$/
 
    """
    md5=\$(samtools view ${input_file} |
           md5sum | awk '{print(\$1)}'); 
    echo \$md5 | awk 'BEGIN{ORS=""}{if(\$1 == "${checksum}") {print("Validated: ${input_file}"); exit(0)} else { print("Failed validation: ${input_file}: ${checksum} != "\$1); exit(1)}}'
    """
}

process checksum_dna_bam {
    input:
      tuple val(input_file), val(checksum)

    output:
       tuple val(input_file), val(checksum), env(md5)

    when:
        input_file =~ /[TN]\.bam$/
 
    """
    md5=\$(samtools view ${input_file} |
           md5sum |
           awk '{print(\$1)}')
    """
}


process validate_rna_bam {
    input:
      tuple val(input_file), val(checksum)

    output:
       tuple val(input_file), val(checksum), env(md5)

    when:
        input_file =~ /fusion\.bam$/

    """
    md5=\$(samtools view ${input_file} | sort |
           md5sum | awk '{print(\$1)}');
    echo \$md5 | awk 'BEGIN{ORS=""}{if(\$1 == "${checksum}") {print("Validated: ${input_file}"); exit(0)} else { print("Failed validation: ${input_file}: ${checksum} != "\$1); exit(1)}}'
    """
}


process checksum_rna_bam {
    input:
      tuple val(input_file), val(checksum)

    output:
       tuple val(input_file), val(checksum), env(md5)

    when:
        input_file =~ /fusion\.bam$/
 
    """
    md5=\$(samtools view ${input_file} | sort |
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
            validate_svdb_query_jumble_vcf &
            validate_id_snps_vcf &
            validate_vcf_gz &
            validate_metrics &
            validate_multiqc &
            validate_mutations_and_coverage &
            validate_samtool_stats &
            validate_collection_of_files &
            validate_genefuse &
            validate_rna_bam &
            validate_dna_bam ) | mix | validate_checksum
}

workflow create_validation_data {
    ch = Channel.fromPath(params.input).splitCsv( header:true, sep:"\t" ).map(row -> [file(row.file), row.checksum])
    ch | (
            validate_vcf &
            validate_svdb_query_jumble_vcf &
            validate_id_snps_vcf &
            validate_vcf_gz &
            validate_metrics &
            validate_multiqc &
            validate_mutations_and_coverage &
            validate_samtool_stats &
            validate_collection_of_files &
            validate_genefuse &
            checksum_dna_bam &
            checksum_rna_bam ) | mix |
            create_checksum_file | collectFile(name: "test.txt", newLine: false, storeDir: "result") | view
}

workflow {
   main:
      validate()
      //create_validation_data()

}

