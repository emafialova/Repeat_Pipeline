#!/usr/bin/env python3
import argparse
import os

def bulk_rename_files(target_dir, old_name, new_name, dry_run=False):
    print(f"Scanning directory: {os.path.abspath(target_dir)}")
    print(f"Target: Renaming '{old_name}' -> '{new_name}'\n")
    
    rename_count = 0
    
    # os.walk systematically steps through every subfolder
    for root, dirs, files in os.walk(target_dir):
        if old_name in files:
            old_file_path = os.path.join(root, old_name)
            new_file_path = os.path.join(root, new_name)
            
            # Extract the chromosome folder name for clean logging
            chrom_folder = os.path.basename(root)
            
            if dry_run:
                print(f"[DRY RUN] Would rename in {chrom_folder}: '{old_name}' -> '{new_name}'")
                rename_count += 1
            else:
                try:
                    os.rename(old_file_path, new_file_path)
                    print(f"SUCCESS [{chrom_folder}]: Renamed file.")
                    rename_count += 1
                except Exception as e:
                    print(f"ERROR [{chrom_folder}]: Could not rename: {e}")

    print(f"\nFinished! Total files processed: {rename_count}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Bulk rename specific execution files across chromosomal directories.")
    parser.add_argument("-d", "--dir", required=True, help="Path to the root organism directory (e.g., ../Results/HS)")
    parser.add_argument("--dry-run", action="store_true", help="Preview the changes without actually renaming files")
    
    args = parser.parse_args()
    
    # Define your exact mapping rules here
    
    OLD_FILENAME_1 = "A1_window_kmers_10k_5k_final_10_10_regions.tsv"
    NEW_FILENAME_1 = "1_Exact_search_kmers.tsv"
    OLD_FILENAME_2 = "output_windows_10k_5k_final_10_10_regions.tsv"
    NEW_FILENAME_2 = "2_Regions_of_interest.tsv"
    OLD_FILENAME_3 = "C_trial_7_limit.tar.gz"
    NEW_FILENAME_3 = "3_Extension_dir.tar.gz"
    OLD_FILENAME_4 = "output_final_nodes_C7.txt"
    NEW_FILENAME_4 = "4_Final_kmers.tsv"
    OLD_FILENAME_5 = "cluster_islands_C7.tsv"
    NEW_FILENAME_5 = "5_Clusters.bed"
    OLD_FILENAME_6 = "cluster_islands_C7_details.tsv"
    NEW_FILENAME_6 = "5_Clusters_details.tsv"
    
    mapping_dict = {
        OLD_FILENAME_1: NEW_FILENAME_1,
        OLD_FILENAME_2: NEW_FILENAME_2,
        OLD_FILENAME_3: NEW_FILENAME_3,
        OLD_FILENAME_4: NEW_FILENAME_4,
        OLD_FILENAME_5: NEW_FILENAME_5,
        OLD_FILENAME_6: NEW_FILENAME_6
    }
    
    for old_name, new_name in mapping_dict.items():
        bulk_rename_files(args.dir, old_name, new_name, dry_run=args.dry_run)