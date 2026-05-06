import argparse
import sys
import os
import glob
import matplotlib.pyplot as plt
from Bio import SeqIO
from matplotlib.lines import Line2D
import statistics

def get_genome_data(fasta_path):
    """
    Reads a FASTA file and returns a dictionary {chr_name: {'length': bp, 'gc': %}}.
    """
    print(f"Reading data from {fasta_path}...")
    data = {}
    try:
        for record in SeqIO.parse(fasta_path, "fasta"):
            seq = record.seq
            gc_content = (seq.count('G') + seq.count('C') + seq.count('g') + seq.count('c')) / len(seq) * 100
            data[record.id] = {'length': len(record), 'gc': gc_content}
    except FileNotFoundError:
        print(f"Error: Fasta file {fasta_path} not found.")
        sys.exit(1)
    return data

def process_directory(base_dir, genome_data_dict, species_label, target_filename="FINAL_clusters.tsv"):
    """
    Scans a directory for chromosome folders, counts clusters, and returns detailed stats.
    """
    print(f"Processing {species_label} results in {base_dir}...")
    
    subfolders = glob.glob(os.path.join(base_dir, "*"))
    results = [] 

    for folder in subfolders:
        if not os.path.isdir(folder):
            continue
            
        chr_name = os.path.basename(folder)
        
        if chr_name not in genome_data_dict:
            print(f"  Skipping folder {chr_name}: Sequence not found in FASTA.")
            continue
            
        length_bp = genome_data_dict[chr_name]['length']
        gc_val = genome_data_dict[chr_name]['gc']
        tsv_path = os.path.join(folder, target_filename)
        
        count_clusters = 0
        if os.path.exists(tsv_path):
            with open(tsv_path, 'r') as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 4:
                        count_clusters += 1
        
        if length_bp > 0:
            length_mbp = length_bp / 1_000_000
            density = count_clusters / length_mbp
            results.append({
                'chr': chr_name, 
                'density': density, 
                'count': count_clusters, 
                'gc': gc_val
            })
            
    return results

def calculate_group_stats(data_list, g1, g2):
    """Calculates averages for the three groups."""
    stats = {
        "Macro": {"gc": [], "count": [], "dens": []},
        "Micro": {"gc": [], "count": [], "dens": []},
        "Dot":   {"gc": [], "count": [], "dens": []}
    }
    
    for item in data_list:
        if item['chr'] in g1:
            target = stats["Macro"]
        elif item['chr'] in g2:
            target = stats["Micro"]
        else:
            target = stats["Dot"]
            
        target["gc"].append(item['gc'])
        target["count"].append(item['count'])
        target["dens"].append(item['density'])
    
    print("\n" + "="*60)
    print(f"{'Group':<10} | {'Avg GC%':<10} | {'Avg Clusters':<12} | {'Avg Dens (Mbp)':<15}")
    print("-" * 60)
    for group, vals in stats.items():
        if not vals["gc"]: continue
        avg_gc = sum(vals["gc"]) / len(vals["gc"])
        avg_count = sum(vals["count"]) / len(vals["count"])
        avg_dens = sum(vals["dens"]) / len(vals["dens"])

        # Calculate standard deviations (default to 0.0 if only 1 item exists)
        sd_gc = statistics.stdev(vals["gc"]) if len(vals["gc"]) > 1 else 0.0
        sd_count = statistics.stdev(vals["count"]) if len(vals["count"]) > 1 else 0.0
        sd_dens = statistics.stdev(vals["dens"]) if len(vals["dens"]) > 1 else 0.0
        
        # Format strings for the table
        gc_str = f"{avg_gc:.2f} ± {sd_gc:.2f}"
        count_str = f"{avg_count:.1f} ± {sd_count:.1f}"
        dens_str = f"{avg_dens:.2f} ± {sd_dens:.2f}"
        
        print(f"{group:<10} | {gc_str:<18} | {count_str:<22} | {dens_str:<22}")
    print("="*85 + "\n")
    print("="*60 + "\n")

def main():
    parser = argparse.ArgumentParser(description="Compare Cluster Density: Chicken vs Human")
    parser.add_argument("--chicken_dir", required=True, help="Path to Chicken results folder")
    parser.add_argument("--chicken_fasta", required=True, help="Path to Chicken genome FASTA")
    parser.add_argument("--second_dir", required=True, help="Path to Human results folder")
    parser.add_argument("--second_fasta", required=True, help="Path to Human genome FASTA")

    args = parser.parse_args()

    # Load Genome Data (Length + GC)
    chicken_info = get_genome_data(args.chicken_fasta)
    second_info = get_genome_data(args.second_fasta)

    # Process Data
    chicken_data = process_directory(args.chicken_dir, chicken_info, "Chicken")
    second_data = process_directory(args.second_dir, second_info, "Human")

    # Define Chromosomal Groups for selected genomes
    chicken_group_1 = ['CP100555.1', 'CP100556.1', 'CP100557.1', 'CP100558.1', 'CP100559.1', 'CP100560.1', 'CP100561.1', 'CP100562.1', 'CP100563.1', 'CP100594.1'] 
    chicken_group_2 = ['CP100564.1', 'CP100565.1', 'CP100566.1', 'CP100567.1', 'CP100568.1', 'CP100569.1', 'CP100571.1', 'CP100572.1', 'CP100573.1', 'CP100574.1', 'CP100575.1', 'CP100576.1', 'CP100577.1', 'CP100578.1', 'CP100579.1', 'CP100580.1', 'CP100581.1', 'CP100582.1', 'CP100587.1'] 
    
    c_color_1, c_color_2, c_color_3 = "#c6dbef", "#6baed6", "#1f77b4"

    other_group_macro_dn = ['NC_088098.1', 'NC_088099.1', 'NC_088100.1', 'NC_088101.1', 'NC_088102.1', 'NC_088103.1', 'NC_088104.1', 'NC_088105.1', 'NC_088106.1', 'NC_088107.1', 'NC_088130.1', 'NC_088132.1'] 
    other_group_micro_dn = ['NC_133029.1', 'NC_133034.1', 'NC_133035.1', 'NC_133036.1', 'NC_133037.1', 'NC_133038.1', 'NC_133039.1', 'NC_133040.1', 'NC_133041.1', 'NC_133042.1', 'NC_133043.1', 'NC_133044.1', 'NC_133045.1', 'NC_133046.1', 'NC_133047.1', 'NC_133048.1', 'NC_133049.1', 'NC_133050.1', 'NC_133051.1']
    #other_group_macro = ['NC_133024.1', 'NC_133025.1', 'NC_133026.1', 'NC_133027.1', 'NC_133028.1', 'NC_133030.1', 'NC_133031.1', 'NC_133032.1', 'NC_133033.1', 'NC_133063.1', 'NC_133064.1']
    other_group_micro = ['NC_133029.1', 'NC_133034.1', 'NC_133035.1', 'NC_133036.1', 'NC_133037.1', 'NC_133038.1', 'NC_133039.1', 'NC_133040.1', 'NC_133041.1', 'NC_133042.1', 'NC_133043.1', 'NC_133044.1', 'NC_133045.1', 'NC_133046.1', 'NC_133047.1', 'NC_133048.1', 'NC_133049.1', 'NC_133050.1', 'NC_133051.1']
    other_group_macro = ['NC_057849.1', 'NC_057850.1', 'NC_057851.1','NC_057852.1','NC_051245.2','NC_051246.2','NC_057853.1','NC_057854.1','NC_057855.1','NC_051250.2','NC_051251.2']
    cm_rest = ['NC_051252.2','NC_051253.2','NC_051254.2','NC_057856.1','NC_057857.1','NC_051257.2','NC_051258.2','NC_051259.2','NC_051260.2','NC_051261.2','NC_051262.2','NC_051263.2','NC_051264.2','NC_057858.1','NC_057859.1','NC_057860.1','NC_051268.2']

    mapping = {key:i+1 for i,key in enumerate(other_group_macro+cm_rest)}
    print(mapping)
    
    # Calculate and Print Stats Table
    calculate_group_stats(chicken_data, chicken_group_1, chicken_group_2)
    #calculate_group_stats(human_data, other_group_macro, other_group_micro)
    calculate_group_stats(second_data, other_group_macro, cm_rest)

    # Colors for selected genomes 
    o_color_1, o_color_2, o_color_3 = "#9dff97", "#47af2b", "#2EBE48"
    o_color_1, o_color_2, o_color_3 = "#f38763", "#e74629", "#941305"
    h_color = "#f69ae7"
    am_color = "#f88a0d"
    dp_color = "#D9FF00"
    cm_color_1 = "#A37EEA"
    cm_color_2 = "#7D3EF0"

    macro_densities = [item['density'] for item in chicken_data if item['chr'] in chicken_group_1]
    micro_densities = [item['density'] for item in chicken_data if item['chr'] in chicken_group_2]
    dot_densities = [item['density'] for item in chicken_data if item['chr'] not in chicken_group_1 + chicken_group_2]

    avg_macro = sum(macro_densities) / len(macro_densities) if macro_densities else 0
    avg_micro = sum(micro_densities) / len(micro_densities) if micro_densities else 0
    avg_dot = sum(dot_densities) / len(dot_densities) if dot_densities else 0
    
    # Plotting Prep
    chicken_data.sort(key=lambda x: x['chr']) 
    #sort human data based on mapping
    second_data.sort(key=lambda x: mapping.get(x['chr'], float('inf')))
    all_chromosomes, all_densities, bar_colors = [], [], []

    all_chromosomes = ["GG_Avg_Macro", "GG_Avg_Micro", "GG_Avg_Dot"]
    all_densities = [avg_macro, avg_micro, avg_dot]
    bar_colors = [c_color_1, c_color_2, c_color_3]
    
    # Chromosome lists for selected genomes
    gg_chromosome_numbers = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,"Z","W"]
    hs_chromosome_numbers = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y"]
    dn_chromosome_numbers = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,33,30,31,"W",32,"Z",34,35,36,37,38,39]
    tg_chromosome_numbers = [1,"1A",2,3,4,"4A",5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,"Z","W"]
    am_chromosome_numbers = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
    dp_chromosome_numbers = [1,2,3,4,5,6,7,8,9,10,11,12,13,14]
    cm_chromosome_numbers = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28]
    
    # For Full Chromosomal Chicken Data
    #for item in chicken_data:
    #    all_chromosomes.append(item['chr'])
    #    all_densities.append(item['density'])
    #    if item['chr'] in chicken_group_1:
    #        bar_colors.append(c_color_1)
    #    elif item['chr'] in chicken_group_2:
    #        bar_colors.append(c_color_2)
    #    else:
    #        bar_colors.append(c_color_3)
    
    for i, item in enumerate(second_data):
        # Try to get the friendly name from dn_chromosome_numbers list
        if i < len(cm_chromosome_numbers):
            friendly_name = str(cm_chromosome_numbers[i])
        else:
            friendly_name = item['chr']
           
        all_chromosomes.append(friendly_name)
        all_densities.append(item['density'])
        if item['chr'] in other_group_macro:
            bar_colors.append(cm_color_1)
        elif item['chr'] in other_group_micro:
            bar_colors.append(o_color_2)
        else:
            bar_colors.append(cm_color_2)

    plt.figure(figsize=(max(10, len(all_chromosomes)*0.3), 8)) 
    plt.bar(all_chromosomes, all_densities, color=bar_colors, edgecolor='black')
        
    # Get current axes and modify the first 3 labels
    ax = plt.gca()
    for i, label in enumerate(ax.get_xticklabels()):
        if i < 3:
            label.set_rotation(90)
            label.set_verticalalignment('top') # Ensures they don't overlap the axis
        else:
            label.set_rotation(0)

    legend_elements = [
        Line2D([0], [0], color=c_color_1, lw=6, label='GG - Avg. Macrochromosomes'),
        Line2D([0], [0], color=c_color_2, lw=6, label='GG - Avg. Microchromosomes'),
        Line2D([0], [0], color=c_color_3, lw=6, label='GG - Avg. Dot chromosomes'),
        #Line2D([0], [0], color=o_color_1, lw=6, label='TG - Macrochromosomes'),
        #Line2D([0], [0], color=o_color_2, lw=6, label='TG - Microchromosomes'),
        #Line2D([0], [0], color=o_color_3, lw=6, label='TG - Dot crochromosomes')
        Line2D([0], [0], color=cm_color_1, lw=6, label='CM - Macrochromosomes'),
        Line2D([0], [0], color=cm_color_2, lw=6, label='CM - Microchromosomes'),
    ]
    plt.legend(handles=legend_elements)
    plt.xlabel("Chromosome", fontsize=16); plt.ylabel("Clusters per Mbp", fontsize=16)
    plt.title("Cluster Density: Chelonia midas", fontsize=20); 
    #plt.xticks(rotation=0, fontsize=10)
    plt.tight_layout()
    plt.savefig("comparison_cluster_density_colored_CM.png")
    print(f"Plot saved.")

if __name__ == '__main__':
    main()