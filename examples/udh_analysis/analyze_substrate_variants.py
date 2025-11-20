#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
노이즈 클러스터(변이체) 상세 분석
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from collections import Counter
from datetime import datetime
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist, squareform

def analyze_variants(tsv_file):
    """변이체 분석"""

    df = pd.read_csv(tsv_file, sep='\t', header=None,
                     names=['sequence_id', 'substrate_site', 'cluster'])

    print("=" * 70)
    print("UDH Substrate Binding Site Variant Analysis")
    print("=" * 70)
    print()

    # Position definitions
    SUBSTRATE_DIRECT_PDB = [75, 111, 112, 113, 136, 165, 174, 258]
    SUBSTRATE_PROXIMAL_PDB = [76, 163, 164, 175]
    all_positions = sorted(SUBSTRATE_DIRECT_PDB + SUBSTRATE_PROXIMAL_PDB)

    # Standard cluster
    standard = df[df['cluster'] == 0]['substrate_site'].values
    variants = df[df['cluster'] == -1]['substrate_site'].values

    print(f"Standard cluster (C0): {len(standard)} sequences")
    print(f"Variant cluster (Noise): {len(variants)} sequences")
    print()

    # Standard consensus
    standard_consensus = []
    for pos in range(12):
        residues = [seq[pos] for seq in standard if pos < len(seq) and seq[pos] != 'X']
        if residues:
            most_common = Counter(residues).most_common(1)[0][0]
            standard_consensus.append(most_common)
    standard_consensus_str = ''.join(standard_consensus)

    print(f"Standard consensus: {standard_consensus_str}")
    print()

    # Analyze variant patterns
    print("=" * 70)
    print("Variant Pattern Analysis")
    print("=" * 70)
    print()

    # Count unique variant patterns
    variant_patterns = Counter([str(seq) for seq in variants])
    print(f"Total unique variant patterns: {len(variant_patterns)}")
    print(f"Most common variants (top 10):")
    for i, (pattern, count) in enumerate(variant_patterns.most_common(10), 1):
        print(f"  {i}. {pattern} ({count} sequences)")
    print()

    # Position-wise variation analysis
    print("=" * 70)
    print("Position-wise Variation in Variants")
    print("=" * 70)
    print()

    variant_aa_by_pos = []
    for pos in range(12):
        residues = [seq[pos] for seq in variants if pos < len(seq) and seq[pos] != 'X']
        counter = Counter(residues)
        variant_aa_by_pos.append(counter)

        pdb_pos = all_positions[pos]
        pos_type = "DIRECT" if pdb_pos in SUBSTRATE_DIRECT_PDB else "PROXIMAL"
        standard_aa = standard_consensus[pos]

        print(f"Position {pos+1} (PDB {pdb_pos}, {pos_type}):")
        print(f"  Standard: {standard_aa}")
        print(f"  Variants: {dict(counter.most_common(5))}")

        # Calculate diversity
        if residues:
            most_common_freq = counter.most_common(1)[0][1]
            diversity = 1 - (most_common_freq / len(residues))
            print(f"  Diversity: {diversity:.3f} (1.0 = maximum)")
        print()

    # Hierarchical clustering of variants
    print("=" * 70)
    print("Hierarchical Clustering of Variants")
    print("=" * 70)
    print()

    if len(variants) > 3:
        # Calculate distance matrix
        def hamming_distance(seq1, seq2):
            return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

        n = len(variants)
        dist_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(i+1, n):
                dist = hamming_distance(variants[i], variants[j])
                dist_matrix[i, j] = dist
                dist_matrix[j, i] = dist

        # Hierarchical clustering
        condensed_dist = pdist(dist_matrix)
        linkage_matrix = linkage(condensed_dist, method='average')

        # Form sub-clusters (cut at distance threshold)
        sub_clusters = fcluster(linkage_matrix, t=3, criterion='distance')

        sub_cluster_counts = Counter(sub_clusters)
        print(f"Number of sub-clusters: {len(sub_cluster_counts)}")
        print(f"Sub-cluster sizes: {dict(sub_cluster_counts.most_common(10))}")
        print()

        # Analyze each sub-cluster
        print("Sub-cluster consensus sequences:")
        for sc_id in sorted(set(sub_clusters)):
            sc_indices = [i for i, c in enumerate(sub_clusters) if c == sc_id]
            sc_seqs = [variants[i] for i in sc_indices]

            if len(sc_seqs) >= 3:  # Only show sub-clusters with 3+ members
                # Calculate consensus
                consensus = []
                for pos in range(12):
                    residues = [seq[pos] for seq in sc_seqs if pos < len(seq) and seq[pos] != 'X']
                    if residues:
                        most_common = Counter(residues).most_common(1)[0][0]
                        consensus.append(most_common)
                    else:
                        consensus.append('X')

                consensus_str = ''.join(consensus)

                # Compare to standard
                differences = sum(c1 != c2 for c1, c2 in zip(standard_consensus_str, consensus_str))

                print(f"  Sub-cluster {sc_id} ({len(sc_seqs)} sequences):")
                print(f"    Consensus: {consensus_str}")
                print(f"    Standard:  {standard_consensus_str}")
                print(f"    Diff:      {''.join(['|' if c1 != c2 else ' ' for c1, c2 in zip(standard_consensus_str, consensus_str)])}")
                print(f"    Hamming distance: {differences}")
                print()

        return df, variants, variant_aa_by_pos, sub_clusters, linkage_matrix, all_positions, SUBSTRATE_DIRECT_PDB, standard_consensus_str

    return df, variants, variant_aa_by_pos, None, None, all_positions, SUBSTRATE_DIRECT_PDB, standard_consensus_str

def create_improved_visualizations(df, variants, variant_aa_by_pos, sub_clusters, linkage_matrix,
                                   all_positions, direct_positions, standard_consensus, output_file):
    """개선된 시각화"""

    fig = plt.figure(figsize=(20, 14))

    # Color scheme
    standard_color = '#2ecc71'  # Green
    variant_color = '#e74c3c'   # Red

    # 1. Cluster distribution
    ax1 = plt.subplot(3, 3, 1)
    standard = df[df['cluster'] == 0]['substrate_site'].values

    sizes = [len(standard), len(variants)]
    labels = [f'Standard\n(C0)\n{len(standard)}', f'Variants\n(Noise)\n{len(variants)}']
    colors = [standard_color, variant_color]

    wedges, texts, autotexts = ax1.pie(sizes, labels=labels, colors=colors,
                                         autopct='%1.1f%%', startangle=90,
                                         textprops={'fontsize': 11, 'fontweight': 'bold'},
                                         explode=[0, 0.1])
    ax1.set_title('Cluster Distribution', fontsize=14, fontweight='bold')

    # 2. Position-wise diversity in variants
    ax2 = plt.subplot(3, 3, 2)

    diversity_scores = []
    for pos in range(12):
        residues = [seq[pos] for seq in variants if pos < len(seq) and seq[pos] != 'X']
        if residues:
            counter = Counter(residues)
            most_common_freq = counter.most_common(1)[0][1]
            diversity = 1 - (most_common_freq / len(residues))
            diversity_scores.append(diversity)
        else:
            diversity_scores.append(0)

    colors_div = [variant_color if all_positions[i] in direct_positions else '#f39c12' for i in range(12)]
    bars = ax2.bar(range(12), diversity_scores, color=colors_div, edgecolor='black', alpha=0.7)
    ax2.set_xticks(range(12))
    ax2.set_xticklabels(all_positions, fontsize=9)
    ax2.set_xlabel('PDB Position', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Diversity Score', fontsize=11, fontweight='bold')
    ax2.set_title('Position-wise Diversity in Variants\n(Red=Direct, Orange=Proximal)',
                  fontsize=12, fontweight='bold')
    ax2.grid(axis='y', alpha=0.3, linestyle='--')
    ax2.set_ylim([0, 1])

    # 3. Most common amino acids at each position (Variants only)
    ax3 = plt.subplot(3, 3, 3)
    ax3.axis('off')
    ax3.set_xlim(0, 10)
    ax3.set_ylim(0, 15)

    ax3.text(5, 14, 'Most Common Residues in Variants', ha='center',
            fontsize=13, fontweight='bold')

    y_pos = 12.5
    for pos in range(12):
        counter = variant_aa_by_pos[pos]
        if counter:
            top3 = counter.most_common(3)
            pdb_pos = all_positions[pos]
            pos_type = "D" if pdb_pos in direct_positions else "P"
            standard_aa = standard_consensus[pos]

            top_str = ', '.join([f"{aa}({cnt})" for aa, cnt in top3])
            ax3.text(0.5, y_pos, f"Pos {pos+1} ({pdb_pos},{pos_type}):", ha='left',
                    fontsize=9, family='monospace')
            ax3.text(3.5, y_pos, f"Std={standard_aa}", ha='left', fontsize=9,
                    family='monospace', color=standard_color, fontweight='bold')
            ax3.text(5, y_pos, f"Var: {top_str}", ha='left', fontsize=9, family='monospace')
            y_pos -= 1

    # 4. Dendrogram of variants (if available)
    ax4 = plt.subplot(3, 3, 4)
    if linkage_matrix is not None and len(variants) > 3:
        dendrogram(linkage_matrix, ax=ax4, no_labels=True, color_threshold=3)
        ax4.set_title(f'Hierarchical Clustering of Variants\n({len(variants)} sequences)',
                     fontsize=12, fontweight='bold')
        ax4.set_xlabel('Variant Index', fontsize=10)
        ax4.set_ylabel('Hamming Distance', fontsize=10)
    else:
        ax4.text(0.5, 0.5, 'Not enough variants\nfor clustering',
                ha='center', va='center', fontsize=12, transform=ax4.transAxes)
        ax4.set_title('Hierarchical Clustering', fontsize=12, fontweight='bold')

    # 5. Sub-cluster sizes (if available)
    ax5 = plt.subplot(3, 3, 5)
    if sub_clusters is not None:
        sub_cluster_counts = Counter(sub_clusters)
        sc_ids = [sc for sc, cnt in sub_cluster_counts.most_common(15)]
        sc_sizes = [sub_cluster_counts[sc] for sc in sc_ids]

        colors_sc = plt.cm.Set3(np.linspace(0, 1, len(sc_ids)))
        bars = ax5.barh(range(len(sc_ids)), sc_sizes, color=colors_sc, edgecolor='black')
        ax5.set_yticks(range(len(sc_ids)))
        ax5.set_yticklabels([f'SC-{sc}' for sc in sc_ids], fontsize=9)
        ax5.set_xlabel('Number of Sequences', fontsize=10, fontweight='bold')
        ax5.set_title('Variant Sub-clusters\n(Top 15)', fontsize=12, fontweight='bold')
        ax5.invert_yaxis()
        ax5.grid(axis='x', alpha=0.3, linestyle='--')

        # Add count labels
        for i, (bar, size) in enumerate(zip(bars, sc_sizes)):
            ax5.text(bar.get_width() + max(sc_sizes)*0.02, bar.get_y() + bar.get_height()/2,
                    str(size), va='center', fontsize=8)
    else:
        ax5.text(0.5, 0.5, 'Sub-clustering not performed',
                ha='center', va='center', fontsize=12, transform=ax5.transAxes)
        ax5.set_title('Variant Sub-clusters', fontsize=12, fontweight='bold')

    # 6. Heatmap: Standard cluster (discrete colors)
    ax6 = plt.subplot(3, 3, 6)

    aa_list = list('ACDEFGHIKLMNPQRSTVWY')
    freq_matrix_std = np.zeros((len(aa_list), 12))

    for pos in range(12):
        residues = [seq[pos] for seq in standard if pos < len(seq) and seq[pos] != 'X']
        counter = Counter(residues)
        total = len(residues)

        for i, aa in enumerate(aa_list):
            freq_matrix_std[i, pos] = counter.get(aa, 0) / total if total > 0 else 0

    # Use discrete colormap
    im = ax6.imshow(freq_matrix_std, cmap='Greens', aspect='auto', vmin=0, vmax=1)
    ax6.set_yticks(range(len(aa_list)))
    ax6.set_yticklabels(aa_list, fontsize=9)
    ax6.set_xticks(range(12))
    ax6.set_xticklabels(all_positions, fontsize=9, rotation=45)
    ax6.set_xlabel('PDB Position', fontsize=10, fontweight='bold')
    ax6.set_ylabel('Amino Acid', fontsize=10, fontweight='bold')
    ax6.set_title(f'Standard Cluster AA Frequency\n(n={len(standard)})',
                 fontsize=12, fontweight='bold')
    cbar = plt.colorbar(im, ax=ax6, fraction=0.046)
    cbar.set_label('Frequency', fontsize=9)

    # 7. Heatmap: Variant cluster (discrete colors)
    ax7 = plt.subplot(3, 3, 7)

    freq_matrix_var = np.zeros((len(aa_list), 12))

    for pos in range(12):
        residues = [seq[pos] for seq in variants if pos < len(seq) and seq[pos] != 'X']
        counter = Counter(residues)
        total = len(residues)

        for i, aa in enumerate(aa_list):
            freq_matrix_var[i, pos] = counter.get(aa, 0) / total if total > 0 else 0

    # Use discrete colormap
    im = ax7.imshow(freq_matrix_var, cmap='Reds', aspect='auto', vmin=0, vmax=1)
    ax7.set_yticks(range(len(aa_list)))
    ax7.set_yticklabels(aa_list, fontsize=9)
    ax7.set_xticks(range(12))
    ax7.set_xticklabels(all_positions, fontsize=9, rotation=45)
    ax7.set_xlabel('PDB Position', fontsize=10, fontweight='bold')
    ax7.set_ylabel('Amino Acid', fontsize=10, fontweight='bold')
    ax7.set_title(f'Variant Cluster AA Frequency\n(n={len(variants)})',
                 fontsize=12, fontweight='bold')
    cbar = plt.colorbar(im, ax=ax7, fraction=0.046)
    cbar.set_label('Frequency', fontsize=9)

    # 8. Conservation comparison
    ax8 = plt.subplot(3, 3, 8)

    std_conservation = []
    var_conservation = []

    for pos in range(12):
        # Standard
        residues = [seq[pos] for seq in standard if pos < len(seq) and seq[pos] != 'X']
        if residues:
            most_common_freq = Counter(residues).most_common(1)[0][1]
            std_conservation.append(most_common_freq / len(residues))
        else:
            std_conservation.append(0)

        # Variants
        residues = [seq[pos] for seq in variants if pos < len(seq) and seq[pos] != 'X']
        if residues:
            most_common_freq = Counter(residues).most_common(1)[0][1]
            var_conservation.append(most_common_freq / len(residues))
        else:
            var_conservation.append(0)

    x = np.arange(12)
    width = 0.35

    bars1 = ax8.bar(x - width/2, std_conservation, width, label='Standard',
                   color=standard_color, alpha=0.8, edgecolor='black')
    bars2 = ax8.bar(x + width/2, var_conservation, width, label='Variants',
                   color=variant_color, alpha=0.8, edgecolor='black')

    ax8.set_xticks(x)
    ax8.set_xticklabels(all_positions, fontsize=9, rotation=45)
    ax8.set_xlabel('PDB Position', fontsize=10, fontweight='bold')
    ax8.set_ylabel('Conservation Score', fontsize=10, fontweight='bold')
    ax8.set_title('Conservation Comparison', fontsize=12, fontweight='bold')
    ax8.legend(fontsize=10)
    ax8.grid(axis='y', alpha=0.3, linestyle='--')
    ax8.set_ylim([0, 1.05])

    # 9. Summary statistics
    ax9 = plt.subplot(3, 3, 9)
    ax9.axis('off')
    ax9.set_xlim(0, 10)
    ax9.set_ylim(0, 10)

    ax9.text(5, 9, 'Summary Statistics', ha='center', fontsize=14, fontweight='bold')

    stats_text = [
        f"Total sequences: {len(df)}",
        f"",
        f"Standard cluster (C0): {len(standard)} ({len(standard)/len(df)*100:.1f}%)",
        f"  Consensus: {standard_consensus}",
        f"  Avg conservation: {np.mean(std_conservation)*100:.1f}%",
        f"",
        f"Variant cluster: {len(variants)} ({len(variants)/len(df)*100:.1f}%)",
        f"  Unique patterns: {len(set([str(seq) for seq in variants]))}",
        f"  Avg conservation: {np.mean(var_conservation)*100:.1f}%",
        f"  Avg diversity: {np.mean(diversity_scores):.3f}",
    ]

    if sub_clusters is not None:
        sub_cluster_counts = Counter(sub_clusters)
        stats_text.append(f"  Sub-clusters: {len(sub_cluster_counts)}")

    y_pos = 7.5
    for line in stats_text:
        if line == "":
            y_pos -= 0.3
        else:
            ax9.text(0.5, y_pos, line, ha='left', fontsize=10, family='monospace')
            y_pos -= 0.6

    plt.suptitle('UDH Substrate Binding Site: Standard vs Variants - Detailed Analysis',
                 fontsize=18, fontweight='bold', y=0.998)
    plt.tight_layout()

    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nImproved visualization saved: {output_file}")

    return output_file

def main():
    tsv_file = Path("udh_substrate_asmc_eps025/groups_0.25_min_2.tsv")

    if not tsv_file.exists():
        print(f"Error: {tsv_file} not found")
        return

    # Analyze variants
    result = analyze_variants(tsv_file)
    if result is None or len(result) < 8:
        print("Error in analysis")
        return

    df, variants, variant_aa_by_pos, sub_clusters, linkage_matrix, all_positions, direct_positions, standard_consensus = result

    # Create improved visualizations
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_file = f"udh_variants_analysis_{timestamp}.png"
    viz_file = create_improved_visualizations(df, variants, variant_aa_by_pos, sub_clusters,
                                              linkage_matrix, all_positions, direct_positions,
                                              standard_consensus, output_file)

    print()
    print("=" * 70)
    print("Variant Analysis Complete!")
    print("=" * 70)
    print(f"Visualization: {viz_file}")
    print("=" * 70)

if __name__ == "__main__":
    main()
