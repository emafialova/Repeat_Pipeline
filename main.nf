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
//params.gff_file = null
params.prefix = "repeats"
params.start = null
params.end = null
params.mapping = null	   
params.org = null          
params.db_path = "genome_clusters.db"

input_sequence = Channel.fromPath(params.sequence, checkIfExists: true)
//gff_file = Channel.fromPath(params.gff_file, checkIfExists: true)

// ---------------------------------------------------------
// 2. Script Channels 
// ---------------------------------------------------------
script_extract_seq  = file("${projectDir}/scripts/A1_0_extract_chr.py")
script_full_chr      = file("${projectDir}/scripts/Approach_1_full_chr.py")
script_analyze_kmers = file("${projectDir}/scripts/Approach_1_analyze_kmers.py")
script_final_nodes   = file("${projectDir}/scripts/A1_2_final_nodes_regions.py")
script_clusters      = file("${projectDir}/scripts/A1_4_trial.py")
//script_yass          = file("${projectDir}/scripts/YASS_ground_truth.py")
//script_sens_spec     = file("${projectDir}/scripts/Sens_spec_calculation.py")
//script_CDS           = file("${projectDir}/scripts/A1_5_CDS_overlap.py")

// C Files
script_c_main        = file("${projectDir}/scripts/A1_2_sec.c")
//script_c_zhash       = file("${projectDir}/zhash-c/src/zhash.c")
dir_c_hash           = file("${projectDir}/scripts/zhash-c")
//script_h_zhash       = file("${projectDir}/zhash-c/src/zhash.h")

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
    //CLUSTERS(FINAL_NODES.out, EXTENSION.out, input_sequence, script_clusters)

    /*
    // G. YASS (Parallel independent step)
    YASS(EXTRACT_SEQ.out, script_yass)

    // H. Stats
    SENS_SPEC(CLUSTERS.out.main_out, YASS.out, script_sens_spec, EXTRACT_SEQ.out)
    */
    // I. CDS Overlap
    //CDS_OVERLAP(script_CDS, CLUSTERS.out.main_out, gff_file)
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
    path "A1_2_sec"

    script:
    """
    # Compile linking both C files. 
    # Because 'h_dep' is in inputs, Nextflow stages it in this dir, 
    # so we don't need complex -I flags.
    gcc -O3 -o A1_2_sec ${c_main} ${zhash_lib_dir}/src/zhash.c -lm
    """
}

// ---------------------------------------------------------
// PROCESS: Step 1 (Full Chr)
// ---------------------------------------------------------
process FULL_CHR {
    publishDir "${params.outdir}", mode: 'copy'
    input:
    path seq
    path script_file

    output:
    path "A1_window_kmers_10k_5k_final_10_10_regions.tsv"

    script:
    """
    python ${script_file} ${seq} A1_window_kmers_10k_5k_final_10_10_regions.tsv
    """
}

// ---------------------------------------------------------
// PROCESS: Step 2 (Analyze Kmers)
// ---------------------------------------------------------
process ANALYZE_KMERS {
    publishDir "${params.outdir}", mode: 'copy'
    input:
    path regions_tsv
    path seq
    path script_file

    output:
    path "output_windows_10k_5k_final_10_10_regions.tsv", emit: main_out
    path "*.png", emit: plots, optional: true

    script:
    """
    python ${script_file} ${regions_tsv} output_windows_10k_5k_final_10_10 ${seq}
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
// PROCESS: Step 4 (Final Nodes)
// ---------------------------------------------------------
process FINAL_NODES {
    publishDir "${params.outdir}", mode: 'copy'
    input:
    path extension_dir // This is a directory
    path script_file

    output:
    path "output_final_nodes_C7_new.txt"

    script:
    """
    python ${script_file} ${extension_dir} output_final_nodes_C7_new.txt --workers ${task.cpus}
    """
}

// ---------------------------------------------------------
// PROCESS: Step 5 (Clusters)
// ---------------------------------------------------------
process CLUSTERS {
    publishDir "${params.outdir}", mode: 'copy' // Save this final file to your actual computer

    input:
    path final_nodes_file
    path extension_dir
    path seq
    path script_file

    output:
    path "cluster_islands_C7_new.tsv", emit: main_out
    path "cluster_islands_C7_new_details.tsv", optional: true

    script:
    """
    python ${script_file} ${final_nodes_file} ${extension_dir} ${seq} ${params.sequence_id} \
        --mapping ${params.mapping} \
        --org "${params.org}" \
        --db_path ${params.db_path} \
        --output cluster_islands_C7_new.tsv \
        --workers ${task.cpus}
    """
}

// ---------------------------------------------------------
// PROCESS: YASS
// ---------------------------------------------------------
process YASS {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path seq
    path script_file

    output:
    path "CPT1B_yass_trial1.tsv"

    script:
    """
    python ${script_file} ${seq} CPT1B_yass_trial1.tsv --diag_tol ${params.diag_tol}
    """
}

// ---------------------------------------------------------
// PROCESS: Sens/Spec
// ---------------------------------------------------------
process SENS_SPEC {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path clusters
    path yass
    path script_file
    path seq

    output:
    //stdout // Capture print output to screen or log
    path "confusion_matrix.txt" 

    script:
    """
    python ${script_file} ${yass} ${clusters} ${seq} "confusion_matrix.txt" 
    """
}

// ---------------------------------------------------------
// PROCESS: CDS
// ---------------------------------------------------------
process CDS_OVERLAP {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path script_file
    path regions_bed_file
    path gff_file

    output:
    //stdout // Capture print output to screen or log
    path "repeats_non_transcript.bed" 
    path "repeats_in_transcript.bed"

    script:
    """
    python ${script_file} ${regions_bed_file} ${gff_file} --prefix ${params.prefix} --chromosome ${params.chromosome_id}
    """
}
