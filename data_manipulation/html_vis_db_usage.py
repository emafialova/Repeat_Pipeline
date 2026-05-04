#!/usr/bin/env python3
"""
DB-Driven HTML visualizer for kmer occurrences.
Uses relative database positions mapped to absolute genomic coordinates.
Fixed: Implemented clash detection and a 12-color high-contrast palette.
Added: Cluster metrics summary table.
Added: In-text highlighting for Core K-mers.
Added: Phylogenetic Tree visualization (Untruncated) + MSA (Single Block).
"""

import argparse
import html
import hashlib
import sqlite3
import json
import io
import tempfile
import subprocess
import sys
from Bio import SeqIO, Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# HElper function
def wrap_text(text, width):
    """Splits a string into chunks of a specific width joined by <br>."""
    if not text or len(text) <= width:
        return text
    return "<br>".join([text[i:i+width] for i in range(0, len(text), width)])

# ---------- DATABASE DATA FETCHING ----------

def get_cluster_data_from_db(db_path, human_id):
    """Retrieves all necessary data for one cluster from the SQLite DB."""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    cursor.execute("SELECT internal_id, chrom, start_pos, end_pos FROM instances WHERE human_id = ?", (human_id,))
    meta = cursor.fetchone()
    if not meta:
        print(f"❌ Error: ID {human_id} not found in database.")
        sys.exit(1)
    
    internal_id, chrom, start, end = meta
    
    # --- FETCH METRICS ---
    cursor.execute("""
        SELECT rep_sequence, gc, cons, med_len, avg_raw_lev, avg_norm_lev, 
               cov, region_size, kmer_num, avg_occ, med_occ 
        FROM cluster_details WHERE internal_id = ?""", (internal_id,))
    details = cursor.fetchone()
    
    # --- FETCH FINAL KMERS ---
    cursor.execute("SELECT kmer_seq, positions, is_rep FROM final_kmer_info WHERE internal_id = ?", (internal_id,))
    rows = cursor.fetchall()
    positions_info = {row[0]: json.loads(row[1]) for row in rows}
    rep_seq = next((row[0] for row in rows if row[2] == 1), None)

    # --- FETCH CORE KMERS ---
    cursor.execute("SELECT core_seq, positions FROM core_kmer_info WHERE internal_id = ?", (internal_id,))
    core_rows = cursor.fetchall()
    core_info = {row[0]: json.loads(row[1]) for row in core_rows}

    # --- FETCH DOMINANT CORE ---
    cursor.execute("""
        SELECT core_kmer, COUNT(*) as c 
        FROM final_kmer_info 
        WHERE internal_id = ? 
        GROUP BY core_kmer 
        ORDER BY c DESC LIMIT 1
    """, (internal_id,))
    dom_row = cursor.fetchone()
    dominant_core = dom_row[0] if dom_row else None
    
    conn.close()
            
    return {
        "chrom": chrom, "start": start, "end": end,
        "rep_seq": rep_seq, "final_seqs": list(positions_info.keys()),
        "positions_info": positions_info,
        "core_info": core_info,
        "metrics": details,
        "dominant_core": dominant_core
    }

# ---------- TREE CONSTRUCTION HELPERS ----------

def run_muscle_msa(records):
    """Perform MSA using MUSCLE via subprocess."""
    with tempfile.TemporaryDirectory() as tmpdir:
        input_fasta = f"{tmpdir}/input.fasta"
        output_fasta = f"{tmpdir}/aligned.fasta"

        with open(input_fasta, "w") as f:
            for rec in records:
                f.write(f">{rec.id}\n{str(rec.seq)}\n")

        # Suppress MUSCLE output to keep console clean
        subprocess.run(["muscle", "-align", input_fasta, "-output", output_fasta], 
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        alignment = AlignIO.read(output_fasta, "fasta")
        return alignment

def generate_tree_and_msa(final_kmers, core_kmer, kmer_to_name):
    """Builds a tree rooted at core_kmer and returns ASCII string + Alignment string."""
    records = []
    if core_kmer:
        name = kmer_to_name.get(core_kmer, "CORE")
        records.append(SeqRecord(Seq(core_kmer), id=name))
    
    for k in final_kmers:
        if k != core_kmer: 
            name = kmer_to_name.get(k, k[:10])
            records.append(SeqRecord(Seq(k), id=name))

    if len(records) < 2:
        return "Not enough sequences to build a tree.", ""

    try:
        # 1. Align
        alignment = run_muscle_msa(records)
        
        # 2. Format MSA (Single Block, Correctly Aligned)
        # We prepare the labels first to calculate the correct padding width
        msa_rows = []
        for rec in alignment:
            label = rec.id
            # Swap placeholder ID for full name BEFORE calculating padding
            if label == "CORE" and core_kmer:
                label = f"{core_kmer}_(Core)"
            msa_rows.append({'label': label, 'seq': str(rec.seq)})
            
        # Calculate max length based on the FINAL labels
        max_id_len = max(len(r['label']) for r in msa_rows)
        
        msa_lines = []
        for row in msa_rows:
            padded_id = row['label'].ljust(max_id_len + 4) # +4 for extra visual spacing
            msa_lines.append(f"{padded_id} {row['seq']}")
        msa_str = "\n".join(msa_lines)

        # 3. Build Tree
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(alignment)
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)

        # 4. Root Tree
        if core_kmer:
            root_clade = next((c for c in tree.find_clades() if c.name == "Core_Motif"), None)
            if root_clade:
                tree.root_with_outgroup(root_clade)
                tree.ladderize()

        # 5. Draw to String
        f = io.StringIO()
        # Increased to 1000 to strictly prevent truncation of long sequences
        Phylo.draw_ascii(tree, file=f, column_width=150) 
        tree_str = f.getvalue()
        
        if core_kmer:
            tree_str = tree_str.replace("CORE", f"{core_kmer} (Core)")
            # Note: We do NOT replace in msa_str here, because we handled it above
            
        return tree_str, msa_str
    except Exception as e:
        return f"Error building tree: {e}. (Is MUSCLE installed?)", ""

# ---------- THE VISUAL CORE ----------

def make_span_blocks_for_lane(chunk_local_start, chunk_local_end, lane_intervals):
    """Creates a single HTML row for a specific lane (rectangles above sequence)."""
    chunk_len = chunk_local_end - chunk_local_start
    mask = [None] * chunk_len 
    
    for (s, e, color) in lane_intervals:
        overlap_s = max(chunk_local_start, s)
        overlap_e = min(chunk_local_end, e)
        if overlap_s < overlap_e:
            for i in range(overlap_s - chunk_local_start, overlap_e - chunk_local_start):
                mask[i] = color

    html_line = []
    i = 0
    while i < chunk_len:
        current_color = mask[i]
        j = i
        while j < chunk_len and mask[j] == current_color:
            j += 1
        width = j - i
        if current_color is None:
            html_line.append("&nbsp;" * width)
        else:
            html_line.append(f'<span style="background:{current_color}; color:{current_color};">{"█" * width}</span>')
        i = j
    return "".join(html_line)

def highlight_sequence_text(chunk_seq, chunk_local_start, chunk_local_end, core_occurrences):
    """Wraps nucleotides in HTML spans to highlight core kmers."""
    seq_len = len(chunk_seq)
    char_colors = [None] * seq_len
    
    for occ in core_occurrences:
        s, e, color = occ['start'], occ['end'], occ['color']
        overlap_s = max(chunk_local_start, s)
        overlap_e = min(chunk_local_end, e)
        if overlap_s < overlap_e:
            for i in range(overlap_s - chunk_local_start, overlap_e - chunk_local_start):
                if i < seq_len: char_colors[i] = color

    html_out = []
    i = 0
    while i < seq_len:
        current_color = char_colors[i]
        j = i
        while j < seq_len and char_colors[j] == current_color:
            j += 1
        segment = html.escape(chunk_seq[i:j])
        if current_color:
            html_out.append(f"<span style='background-color:{current_color}; font-weight:bold;'>{segment}</span>")
        else:
            html_out.append(segment)
        i = j
    return "".join(html_out)

def get_core_colors(core_seqs):
    """Assigns colors to core kmers."""
    if len(core_seqs) == 1:
        return {core_seqs[0]: "#F90707"} 
    MULTI_PALETTE = ["#F90303", "#3805F2", "#01FB01", "#F57B00", "#7D03F7", "#F8F802"]
    return {k: MULTI_PALETTE[i % len(MULTI_PALETTE)] for i, k in enumerate(sorted(core_seqs))}

# ---------- MAIN HTML BUILD ----------

def build_html(fasta_path, human_id, org_name, db_data, out_html, chunk_size=100):
    chrom, start, stop = db_data['chrom'], db_data['start'], db_data['end']+1
    final_seqs = sorted(db_data['final_seqs'])
    positions_info = db_data['positions_info']
    core_info = db_data['core_info']
    metrics = db_data['metrics']
    rep_seq = db_data['rep_seq'] 

    kmer_to_name = {}
    other_kmers = [k for k in final_seqs if k != rep_seq]
    kmer_to_name[rep_seq] = "Representative_kmer"
    for idx, k in enumerate(other_kmers, 1):
        kmer_to_name[k] = f"Final_kmer_{idx}"
    if db_data['dominant_core'] not in kmer_to_name:
        kmer_to_name[db_data['dominant_core']] = "Core_Motif"
    
    # Colors
    #PASTEL_PALETTE = ["#FFB3BA", "#F1BAFF", "#FFFFBA", "#BAFFC9", "#BAE1FF", "#D4A5A5", "#FFC8A2", "#E2F0CB", "#B5EAD7", "#C7CEEA", "#F3D1DC", "#A8E6CF"]
    ROCKET_CATEGORICAL = [
        "#fde293", "#e26d8d", "#fba973", "#c882a4", 
        "#ffeda0", "#fa7b67", "#fbb4b9", "#fd8d3c", 
        "#df65b0", "#ffc5a1", "#a4739d", "#f4858e"
    ]
    #kmer_colors = {k: ROCKET_CATEGORICAL[i % len(ROCKET_CATEGORICAL)] for i, k in enumerate(final_seqs)}
    kmer_colors = {}
    for i, k in enumerate(final_seqs):
        if i < len(ROCKET_CATEGORICAL):
            # 1. Use your predefined Rocket colors for the first 12
            kmer_colors[k] = ROCKET_CATEGORICAL[i]
        else:
            # 2. For 13+, constrain colors to the "Rocket" zone (Purple to Yellow)
            # The safe window is from 270 degrees up to 420 degrees (which wraps to 60)
            # We use the Golden Ratio fraction (0.618) to pick widely separated hues inside this window
            fraction = (i * 0.6180339887) % 1.0
            hue = int(270 + (fraction * 150)) % 360
            
            # Alternate the lightness slightly (80% and 90%) to keep adjacent colors distinct
            lightness = 80 if i % 2 == 0 else 90
            
            kmer_colors[k] = f"hsl({hue}, 80%, {lightness}%)"
    core_colors = get_core_colors(list(core_info.keys()))

    record_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    if chrom not in record_dict:
        print(f"❌ FASTA Error: Could not find '{chrom}'.")
        sys.exit(1)
    full_seq = str(record_dict[chrom].seq).upper()
    
    html_parts = [
        "<!doctype html><html><head><meta charset='utf-8'>",
        f"<title>Cluster {human_id}</title>",
        "<style>",
        "body { font-family: 'Courier New', monospace; padding: 40px; background: #fff; }",
        ".header-box { background: #f8f9fa; border: 1px solid #dee2e6; padding: 20px; margin-bottom: 30px; border-radius: 8px; }",
        ".chunk-container { margin-bottom: 40px; border-bottom: 1px solid #f0f0f0; padding-bottom: 15px; }",
        ".kmer-row { line-height: 1.0; font-size: 14px; white-space: pre; height: 14px; margin-bottom: 2px; }",
        ".seq-row { line-height: 1.2; font-size: 14px; white-space: pre; background: #fdfdfd; }", 
        ".coord { color: #aaa; font-size: 11px; display:inline-block; width:10ch; margin-right: 15px; text-align: right; }",
        ".legend { border: 1px solid #ddd; padding: 15px; margin-top: 20px; background: #fafafa; }",
        ".legend-item { display:inline-block; margin-right:20px; margin-bottom: 8px; padding: 4px 10px; border-radius: 4px; border: 1px solid #ccc; font-size: 13px; max-width: 90vw; word-break: break-all; vertical-align: top; }",
        ".metrics-table { border-collapse: collapse; width: 100%; margin-top: 50px; font-family: monospace; }",
        ".metrics-table td { border: 1px solid #ddd; padding: 12px; text-align: left; font-family: monospace; white-space: pre-wrap; word-break: normal;}",
        ".metrics-table th, .metrics-table td { border: 1px solid #ddd; padding: 12px; text-align: left; }",
        ".metrics-table th { background-color: #83CBEB; color: black; }",
        ".metrics-table tr:nth-child(even) { background-color: #f2f2f2; }",
        ".tree-box { background: #fdfdfd; border: 1px solid #eee; padding: 20px; margin-top: 20px; overflow-x: auto; white-space: pre; font-family: monospace; }",
        "</style></head><body>"
    ]

    html_parts.append(f"""<div class='header-box'>
        <h1 style='margin:0 0 10px 0; color:#2c3e50;'>Cluster Visualization: {human_id}</h1>
        <table style='width:100%; border-collapse: collapse; font-family: sans-serif;'>
            <tr><td><strong>Organism:</strong> {org_name}</td><td><strong>Chromosome:</strong> {chrom}</td></tr>
            <tr><td><strong>Genomic Range:</strong> {chrom}:{start:}-{stop:}</td><td><strong>Total Span:</strong> {stop-start:,} bp</td></tr>
        </table></div>""")

    # Legend
    html_parts.append("<div class='legend'><strong>Final Motif Legend (Rectangles):</strong><br>")
    for k in final_seqs:
        name = kmer_to_name[k]
        wrapped_k = wrap_text(k, chunk_size)
        label = f"<b>{name}</b><br><b>{wrapped_k}</b>" if k == db_data['rep_seq'] else f"{name}<br>{wrapped_k}"
        html_parts.append(f"<span class='legend-item' style='background:{kmer_colors[k]};'>█ {label}</span>")
    html_parts.append("<br><br><strong>Core Motif Legend (Text Highlight):</strong><br>")
    for c_seq, c_col in core_colors.items():
        wrapped_c = wrap_text(c_seq, chunk_size)
        html_parts.append(f"<span class='legend-item' style='background:{c_col};'>{wrapped_c}</span>")
    html_parts.append("</div><hr style='border:none; border-top: 2px double #eee; margin: 40px 0;'>")

    # Data Prep
    all_final_occs = []
    for k, positions in positions_info.items():
        for p in positions: all_final_occs.append({'start': p, 'end': p + len(k), 'color': kmer_colors[k]})
    all_final_occs.sort(key=lambda x: x['start'])

    all_core_occs = []
    for c_seq, positions in core_info.items():
        for p in positions: all_core_occs.append({'start': p, 'end': p + len(c_seq), 'color': core_colors[c_seq]})
    all_core_occs.sort(key=lambda x: x['start'])

    # Visualization Loop
    total_span = stop - start
    for i in range(0, total_span, chunk_size):
        c_local_start = i
        c_local_end = min(total_span, i + chunk_size)
        c_abs_start = start + c_local_start
        chunk_str_raw = full_seq[c_abs_start : start + c_local_end]

        chunk_final_occs = [o for o in all_final_occs if o['start'] < c_local_end and o['end'] > c_local_start]
        lanes = []; lane_data = [] 
        for occ in chunk_final_occs:
            assigned = False
            for idx, end_pos in enumerate(lanes):
                if occ['start'] >= end_pos:
                    lanes[idx] = occ['end']; lane_data[idx].append((occ['start'], occ['end'], occ['color'])); assigned = True; break
            if not assigned: lanes.append(occ['end']); lane_data.append([(occ['start'], occ['end'], occ['color'])])

        chunk_core_occs = [o for o in all_core_occs if o['start'] < c_local_end and o['end'] > c_local_start]
        highlighted_seq_html = highlight_sequence_text(chunk_str_raw, c_local_start, c_local_end, chunk_core_occs)

        html_parts.append("<div class='chunk-container'>")
        for lane_intervals in lane_data:
            row_html = make_span_blocks_for_lane(c_local_start, c_local_end, lane_intervals)
            html_parts.append(f"<div class='kmer-row'><span class='coord'>{c_abs_start}</span>{row_html}</div>")
        html_parts.append(f"<div class='seq-row'><span class='coord'>{c_abs_start}</span>{highlighted_seq_html}</div>")
        html_parts.append("</div>")

    # Metrics Table
    if metrics:
        m_labels = ["Representative", "GC %", "Conservation", "Med Length", "Avg Raw Lev", "Avg Norm Lev", "Coverage", "Region Size", "Kmer Num", "Avg Occ", "Med Occ"]
        html_parts.append("<h2>Cluster Metrics Summary</h2>")
        html_parts.append("<table class='metrics-table'><tr>")
        for label in m_labels: html_parts.append(f"<th>{label}</th>")
        html_parts.append("</tr><tr>")
        #html_parts.append(f"<td>{metrics[0]}</td>") 
        wrapped_rep = wrap_text(str(metrics[0]), 70)
        html_parts.append(f"<td>{wrapped_rep}</td>")
        for val in metrics[1:]:
            display_val = f"{val:.2f}" if isinstance(val, float) else val
            html_parts.append(f"<td>{display_val}</td>")
        html_parts.append("</tr></table>")

    # --- PHYLOGENETICS SECTION ---
    if final_seqs:
        print("🌳 Building tree and MSA...")
        tree_ascii, msa_str = generate_tree_and_msa(final_seqs, db_data['dominant_core'], kmer_to_name)
        
        # Combine all unique sequences we want to colorize
        all_seqs = list(final_seqs)
        if db_data['dominant_core'] and db_data['dominant_core'] not in all_seqs:
            all_seqs.append(db_data['dominant_core'])
            
        # Sort by name length descending to avoid partial string replacement bugs
        all_seqs.sort(key=lambda x: len(kmer_to_name[x]), reverse=True)
        
        for seq in all_seqs:
            name = kmer_to_name[seq]
            # Fetch color, falling back to core_colors if it's the core motif
            color = kmer_colors.get(seq, core_colors.get(seq, "#cccccc"))
            
            # Wrap the name in a colored span
            styled_name = f"<span style='background-color:{color}; border-radius:3px; padding:0 3px;'>{name}</span>"
            
            if seq == rep_seq:
                styled_name = f"<b>{styled_name}</b>"
                
            # Perform the replacement in both text outputs
            tree_ascii = tree_ascii.replace(name, styled_name)
            msa_str = msa_str.replace(name, styled_name)

        html_parts.append("<h2>Phylogenetic Tree (Neighbor-Joining)</h2>")
        html_parts.append(f"<div class='tree-box'>{tree_ascii}</div>")
        
        html_parts.append("<h2>Multiple Sequence Alignment (MUSCLE)</h2>")
        html_parts.append(f"<div class='tree-box'>{msa_str}</div>")

    html_parts.append("</body></html>")

    with open(out_html, "w", encoding="utf8") as f:
        f.write("".join(html_parts))
    print(f"✅ Visualization generated: {out_html}")

def main():
    parser = argparse.ArgumentParser(description="DB-Driven HTML Cluster Reports.")
    parser.add_argument("--db", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--id", required=True)
    parser.add_argument("--org", required=True)
    parser.add_argument("--out", default="report.html")
    parser.add_argument("--chunk", type=int, default=100)
    args = parser.parse_args()

    data = get_cluster_data_from_db(args.db, args.id)
    build_html(args.fasta, args.id, args.org, data, args.out, args.chunk)

if __name__ == "__main__":
    main()