#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ---------------------------------------------------------
// 1. Pipeline Parameters & Inputs
// ---------------------------------------------------------
params.sequence = null
params.outdir = "results"
params.allowed_HD_error = 15
params.extension_dir = "C_output_dir"
params.diag_tol = 2000
params.sequence_id = null
params.chromosome_id = null
params.prefix = "repeats"
params.start = null
params.end = null
params.mapping = null	   
params.org = null          
params.db_path = "genome_clusters.db"

input_sequence = Channel.fromPath(params.sequence, checkIfExists: true)

// ---------------------------------------------------------
// 2. Script Channels 
// ---------------------------------------------------------
script_extract_seq  = file("${projectDir}/scripts/0_extract_chr.py")
script_full_chr      = file("${projectDir}/scripts/1_Exact_search.py")
script_analyze_kmers = file("${projectDir}/scripts/2_Region_creation.py")
script_final_nodes   = file("${projectDir}/scripts/4_Final_kmers.py")
script_clusters      = file("${projectDir}/scripts/5_Clustering.py")

// C Files
script_c_main        = file("${projectDir}/scripts/3_Extension.c")
dir_c_hash           = file("${projectDir}/scripts/zhash-c")

// ---------------------------------------------------------
// 3. Workflow Definition
// ---------------------------------------------------------
workflow {
    // 0. Extract Sequence of Interest
    EXTRACT_SEQ(script_extract_seq, input_sequence)

    // A. Compile the C program first
    COMPILE_C(script_c_main, dir_c_hash)

    // B. Run Python Step 1
    FULL_CHR(EXTRACT_SEQ.out, script_full_chr)

    // C. Run Python Step 2 (Takes output of Step 1)
    ANALYZE_KMERS(FULL_CHR.out, EXTRACT_SEQ.out, script_analyze_kmers)

    // D. Run C Step (Takes Output of Step 2 + Compiled Binary + Sequence)
    EXTENSION(ANALYZE_KMERS.out.main_out, COMPILE_C.out, EXTRACT_SEQ.out)

    // E. Run Final Nodes (Takes Output of Extension)
    FINAL_NODES(EXTENSION.out, script_final_nodes)

    // F. Clusters
    CLUSTERS(FINAL_NODES.out, EXTENSION.out, EXTRACT_SEQ.out, script_clusters)
}

// ---------------------------------------------------------
// PROCESS: Step 0 (Extract sequence)
// ---------------------------------------------------------
process EXTRACT_SEQ {
    publishDir "${params.outdir}", mode: 'copy'
    input:
    path script_file
    path seq

    output:
    path "${params.chromosome_id}_seq.fasta"

    script:
    def start_flag = params.start ? "--start ${params.start}" : ""
    def end_flag   = params.end   ? "--end ${params.end}"   : ""
    def output_filename = "${params.chromosome_id}_seq.fasta"
    """
    python ${script_file} ${seq} ${output_filename} ${params.sequence_id} ${start_flag} ${end_flag}
    """
}

// ---------------------------------------------------------
// PROCESS: Compile C
// ---------------------------------------------------------
process COMPILE_C {
    cache true 

    input:
    path c_main
    path zhash_lib_dir
    //path h_dep
    
    output:
    path "3_Extension"

    script:
    """
    # Compile linking both C files. 
    gcc -O3 -o 3_Extension ${c_main} ${zhash_lib_dir}/src/zhash.c -lm
    """
}

// ---------------------------------------------------------
// PROCESS: Step 1 (Exact Search for Kmers)
// ---------------------------------------------------------
process FULL_CHR {
    publishDir "${params.outdir}", mode: 'copy'
    input:
    path seq
    path script_file

    output:
    path "1_Exact_search_kmers.tsv"

    script:
    """
    python ${script_file} ${seq} 1_Exact_search_kmers.tsv
    """
}

// ---------------------------------------------------------
// PROCESS: Step 2 (Region Creation)
// ---------------------------------------------------------
process ANALYZE_KMERS {
    publishDir "${params.outdir}", mode: 'copy'
    input:
    path regions_tsv
    path seq
    path script_file

    output:
    path "2_Regions_of_interest.tsv", emit: main_out
    path "*.png", emit: plots, optional: true

    script:
    """
    python ${script_file} ${regions_tsv} 2_Regions_of_interest ${seq}
    """
}

// ---------------------------------------------------------
// PROCESS: Step 3 (Extension - C Program)
// ---------------------------------------------------------
process EXTENSION {
    publishDir "${params.outdir}", mode: 'copy'
    input:
    path analyzed_tsv
    path exe           // The compiled binary from COMPILE_C
    path seq
    
    output:
    path "${params.extension_dir}" 

    script:
    """
    # Ensure binary is executable
    chmod +x ${exe}
    
    # Make the output directory
    mkdir -p ${params.extension_dir}
    
    # Run the C program
    ./${exe} ${analyzed_tsv} ${seq} ${params.extension_dir} ${params.allowed_HD_error}
    """
}

// ---------------------------------------------------------
// PROCESS: Step 4 (Final Kmer Extraction)
// ---------------------------------------------------------
process FINAL_NODES {
    publishDir "${params.outdir}", mode: 'copy'
    input:
    path extension_dir // This is a directory
    path script_file

    output:
    path "4_Final_kmers.txt"

    script:
    """
    python ${script_file} ${extension_dir} 4_Final_kmers.txt --workers ${task.cpus}
    """
}

// ---------------------------------------------------------
// PROCESS: Step 5 (Clustering)
// ---------------------------------------------------------
process CLUSTERS {
    publishDir "${params.outdir}", mode: 'copy' // Save this final file to your actual computer

    input:
    path final_nodes_file
    path extension_dir
    path seq
    path script_file

    output:
    path "5_Clusters.bed", emit: main_out
    path "5_Clusters_details.tsv", optional: true

    script:
    """
    python ${script_file} ${final_nodes_file} ${extension_dir} ${seq} ${params.sequence_id} \
        --mapping ${params.mapping} \
        --org "${params.org}" \
        --db_path ${params.db_path} \
        --output 5_Clusters.bed \
        --workers ${task.cpus}
    """
}