#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ASMC 실제 실행 및 결과 시각화"""

import subprocess
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime

def create_test_data():
    """테스트 데이터 생성"""
    print("Creating test data...")
    test_dir = Path("test_data")
    test_dir.mkdir(exist_ok=True)

    # More detailed PDB files with different structures
    pdb_templates = [
        # Protein 1 - Alpha helix rich
        """HEADER    TEST PROTEIN 1
ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00 20.00           N
ATOM      2  CA  ALA A   1      11.450  10.000  10.000  1.00 20.00           C
ATOM      3  C   ALA A   1      12.000  11.420  10.000  1.00 20.00           C
ATOM      4  O   ALA A   1      11.300  12.400  10.000  1.00 20.00           O
ATOM      5  CB  ALA A   1      11.950   9.200  11.200  1.00 20.00           C
ATOM      6  N   VAL A   2      13.300  11.550  10.000  1.00 20.00           N
ATOM      7  CA  VAL A   2      13.950  12.850  10.000  1.00 20.00           C
ATOM      8  C   VAL A   2      15.450  12.750  10.000  1.00 20.00           C
ATOM      9  O   VAL A   2      16.100  11.700  10.000  1.00 20.00           O
ATOM     10  CB  VAL A   2      13.500  13.650  11.200  1.00 20.00           C
TER
END
""",
        # Protein 2 - Beta sheet rich
        """HEADER    TEST PROTEIN 2
ATOM      1  N   GLY A   1      20.000  20.000  20.000  1.00 25.00           N
ATOM      2  CA  GLY A   1      21.450  20.000  20.000  1.00 25.00           C
ATOM      3  C   GLY A   1      22.000  21.420  20.000  1.00 25.00           C
ATOM      4  O   GLY A   1      21.300  22.400  20.000  1.00 25.00           O
ATOM      5  N   SER A   2      23.300  21.550  20.000  1.00 25.00           N
ATOM      6  CA  SER A   2      23.950  22.850  20.000  1.00 25.00           C
ATOM      7  C   SER A   2      25.450  22.750  20.000  1.00 25.00           C
ATOM      8  O   SER A   2      26.100  21.700  20.000  1.00 25.00           O
ATOM      9  CB  SER A   2      23.500  23.650  21.200  1.00 25.00           C
ATOM     10  OG  SER A   2      24.100  24.850  21.200  1.00 25.00           O
TER
END
""",
        # Protein 3 - Mixed structure
        """HEADER    TEST PROTEIN 3
ATOM      1  N   LEU A   1      30.000  30.000  30.000  1.00 30.00           N
ATOM      2  CA  LEU A   1      31.450  30.000  30.000  1.00 30.00           C
ATOM      3  C   LEU A   1      32.000  31.420  30.000  1.00 30.00           C
ATOM      4  O   LEU A   1      31.300  32.400  30.000  1.00 30.00           O
ATOM      5  CB  LEU A   1      31.950  29.200  31.200  1.00 30.00           C
ATOM      6  N   THR A   2      33.300  31.550  30.000  1.00 30.00           N
ATOM      7  CA  THR A   2      33.950  32.850  30.000  1.00 30.00           C
ATOM      8  C   THR A   2      35.450  32.750  30.000  1.00 30.00           C
ATOM      9  O   THR A   2      36.100  31.700  30.000  1.00 30.00           O
ATOM     10  CB  THR A   2      33.500  33.650  31.200  1.00 30.00           C
TER
END
"""
    ]

    # Create PDB files
    for i, pdb_content in enumerate(pdb_templates, 1):
        pdb_file = test_dir / f"protein{i}.pdb"
        pdb_file.write_text(pdb_content, encoding='utf-8')

    # Create references.txt
    refs_content = "\n".join([
        str((test_dir / f"protein{i}.pdb").absolute())
        for i in range(1, 3)
    ])
    (test_dir / "references.txt").write_text(refs_content, encoding='utf-8')

    # Create sequences.fasta with more variation
    fasta_content = """>Protein_1_alpha
MKVLWAALLVTFLAGCQAKVEQAVETEPEPELRQQTEWQSGQRWELALGRFWDYLRWVQT
LSEQVQEELLSSQVTQELRALMDETMKELKAYKSELEEQLTPVA

>Protein_2_beta
GSGSGSGSGSGSGSGSGSGAGACGQEARPEAVGQHWEALGRFWDYLRWVQTLSEQVQEE
LLSSQVTQELRALMDETMKELKAYKSELEEQLTPVA

>Protein_3_mixed
LKTLWAALLVTFLAGCQAKVEQAVETEPEPELRQQTEWQSGQRWELALGRFWDYLRWVQT
LSEQVQEELLSSQVTQELRALMDETMKELKAYKSELEEQLTPVAG

>Protein_4_variant
MKVLWAALLVTFLAGCQAKVEQAVETEPEPELRQQTEWQSGQRWELKLGRFWDYLRWVQT
LSEQVQEELLSSQVTQELRALMDETMKELKAYKSELEEQLTPVAGG

>Protein_5_mutant
AAALWAALLVTFLAGCQAKVEQAVETEPEPELRQQTEWQSGQRWELALGRFWDYLRWVQT
LSEQVQEELLSSQVTQELRALMDETMKELKAYKSELEEQLTPVAGGG
"""
    (test_dir / "sequences.fasta").write_text(fasta_content, encoding='utf-8')

    # Create models.txt with correct references
    models_content = "\n".join([
        f"{(test_dir / f'protein{i}.pdb').absolute()}\t{(test_dir / 'protein1.pdb').absolute()}"
        for i in range(1, 4)
    ])
    (test_dir / "models.txt").write_text(models_content, encoding='utf-8')

    # Create pocket.txt
    pocket_content = "\n".join([
        f"{(test_dir / f'protein{i}.pdb').absolute()}\tA\t1,2"
        for i in range(1, 3)
    ])
    (test_dir / "pocket.txt").write_text(pocket_content, encoding='utf-8')

    print(f"Test data created in {test_dir.absolute()}")
    return test_dir

def run_simple_analysis():
    """간단한 서열 분석 실행"""
    print("\n" + "="*60)
    print("Running simple sequence analysis...")
    print("="*60)

    test_dir = Path("test_data")

    # Read sequences
    sequences = {}
    with open(test_dir / "sequences.fasta", 'r') as f:
        current_seq = ""
        current_name = ""
        for line in f:
            if line.startswith(">"):
                if current_name:
                    sequences[current_name] = current_seq
                current_name = line[1:].strip()
                current_seq = ""
            else:
                current_seq += line.strip()
        if current_name:
            sequences[current_name] = current_seq

    # Calculate simple similarity matrix
    similarity_matrix = []
    seq_names = list(sequences.keys())

    for seq1 in seq_names:
        row = []
        for seq2 in seq_names:
            # Simple identity calculation
            s1, s2 = sequences[seq1], sequences[seq2]
            min_len = min(len(s1), len(s2))
            identical = sum(1 for i in range(min_len) if s1[i] == s2[i])
            similarity = identical / min_len if min_len > 0 else 0
            row.append(similarity)
        similarity_matrix.append(row)

    return seq_names, np.array(similarity_matrix), sequences

def create_clustering_results(seq_names, similarity_matrix):
    """간단한 클러스터링 수행"""
    print("\n" + "="*60)
    print("Performing clustering...")
    print("="*60)

    # Simple hierarchical clustering based on similarity
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform

    # Convert similarity to distance
    distance_matrix = 1 - similarity_matrix
    condensed_dist = squareform(distance_matrix)

    # Perform clustering
    linkage_matrix = linkage(condensed_dist, method='average')
    clusters = fcluster(linkage_matrix, 0.3, criterion='distance')

    # Group sequences by cluster
    cluster_groups = {}
    for i, cluster_id in enumerate(clusters):
        if cluster_id not in cluster_groups:
            cluster_groups[cluster_id] = []
        cluster_groups[cluster_id].append(seq_names[i])

    print(f"Found {len(cluster_groups)} clusters")
    for cluster_id, members in cluster_groups.items():
        print(f"  Cluster {cluster_id}: {', '.join(members[:3])}...")

    return clusters, linkage_matrix, cluster_groups

def visualize_results(seq_names, similarity_matrix, clusters, linkage_matrix, sequences):
    """결과 시각화"""
    print("\n" + "="*60)
    print("Creating visualizations...")
    print("="*60)

    fig = plt.figure(figsize=(16, 12))

    # 1. Similarity heatmap
    ax1 = plt.subplot(2, 3, 1)
    im = ax1.imshow(similarity_matrix, cmap='RdYlGn', vmin=0, vmax=1)
    ax1.set_xticks(range(len(seq_names)))
    ax1.set_yticks(range(len(seq_names)))
    ax1.set_xticklabels([name[:15] for name in seq_names], rotation=45, ha='right')
    ax1.set_yticklabels([name[:15] for name in seq_names])
    ax1.set_title('Sequence Similarity Matrix')
    plt.colorbar(im, ax=ax1)

    # 2. Dendrogram
    ax2 = plt.subplot(2, 3, 2)
    from scipy.cluster.hierarchy import dendrogram
    dendrogram(linkage_matrix, labels=[name[:15] for name in seq_names],
               ax=ax2, orientation='top')
    ax2.set_title('Hierarchical Clustering Dendrogram')
    ax2.set_xlabel('Sequence')
    ax2.set_ylabel('Distance')
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45, ha='right')

    # 3. Cluster distribution
    ax3 = plt.subplot(2, 3, 3)
    unique_clusters, counts = np.unique(clusters, return_counts=True)
    colors = plt.cm.Set3(np.linspace(0, 1, len(unique_clusters)))
    ax3.bar(unique_clusters, counts, color=colors)
    ax3.set_xlabel('Cluster ID')
    ax3.set_ylabel('Number of Sequences')
    ax3.set_title('Cluster Size Distribution')
    ax3.set_xticks(unique_clusters)

    # 4. Sequence length distribution
    ax4 = plt.subplot(2, 3, 4)
    seq_lengths = [len(seq) for seq in sequences.values()]
    ax4.hist(seq_lengths, bins=20, edgecolor='black', alpha=0.7)
    ax4.set_xlabel('Sequence Length')
    ax4.set_ylabel('Count')
    ax4.set_title('Sequence Length Distribution')
    ax4.axvline(np.mean(seq_lengths), color='red', linestyle='--',
                label=f'Mean: {np.mean(seq_lengths):.1f}')
    ax4.legend()

    # 5. Amino acid composition
    ax5 = plt.subplot(2, 3, 5)
    aa_counts = {}
    for seq in sequences.values():
        for aa in seq:
            aa_counts[aa] = aa_counts.get(aa, 0) + 1

    sorted_aa = sorted(aa_counts.items(), key=lambda x: x[1], reverse=True)[:20]
    aa_labels, aa_values = zip(*sorted_aa)
    ax5.bar(range(len(aa_labels)), aa_values, color='steelblue')
    ax5.set_xticks(range(len(aa_labels)))
    ax5.set_xticklabels(aa_labels)
    ax5.set_xlabel('Amino Acid')
    ax5.set_ylabel('Frequency')
    ax5.set_title('Top 20 Amino Acid Frequencies')

    # 6. Cluster similarity network
    ax6 = plt.subplot(2, 3, 6)
    # Create network visualization using scatter plot
    np.random.seed(42)
    positions = {}
    for i, name in enumerate(seq_names):
        angle = 2 * np.pi * i / len(seq_names)
        positions[i] = (np.cos(angle), np.sin(angle))

    # Plot nodes
    x_pos = [positions[i][0] for i in range(len(seq_names))]
    y_pos = [positions[i][1] for i in range(len(seq_names))]
    colors = [plt.cm.Set3(c/max(clusters)) for c in clusters]
    ax6.scatter(x_pos, y_pos, c=colors, s=200, alpha=0.8, edgecolors='black')

    # Add labels
    for i, name in enumerate(seq_names):
        ax6.annotate(name[:10], (x_pos[i], y_pos[i]),
                    ha='center', va='center', fontsize=8)

    # Draw edges for high similarity
    for i in range(len(seq_names)):
        for j in range(i+1, len(seq_names)):
            if similarity_matrix[i, j] > 0.7:  # Only show high similarity
                ax6.plot([x_pos[i], x_pos[j]], [y_pos[i], y_pos[j]],
                        'k-', alpha=0.2, linewidth=similarity_matrix[i, j]*2)

    ax6.set_xlim(-1.5, 1.5)
    ax6.set_ylim(-1.5, 1.5)
    ax6.set_aspect('equal')
    ax6.set_title('Sequence Similarity Network')
    ax6.axis('off')

    plt.suptitle('ASMC Analysis Results', fontsize=16, fontweight='bold')
    plt.tight_layout()

    # Save figure
    output_file = f"asmc_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Visualization saved to {output_file}")

    plt.show()

    return output_file

def create_summary_report(seq_names, similarity_matrix, cluster_groups, sequences):
    """결과 요약 리포트 생성"""
    print("\n" + "="*60)
    print("ASMC Analysis Summary Report")
    print("="*60)

    report = []
    report.append(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report.append(f"Number of sequences analyzed: {len(sequences)}")
    report.append(f"Number of clusters found: {len(cluster_groups)}")
    report.append("")

    report.append("Cluster Summary:")
    for cluster_id, members in cluster_groups.items():
        report.append(f"  Cluster {cluster_id}: {len(members)} members")
        for member in members[:3]:  # Show first 3 members
            report.append(f"    - {member}")
        if len(members) > 3:
            report.append(f"    ... and {len(members)-3} more")

    report.append("")
    report.append("Similarity Statistics:")
    report.append(f"  Average similarity: {np.mean(similarity_matrix):.3f}")
    report.append(f"  Minimum similarity: {np.min(similarity_matrix):.3f}")
    report.append(f"  Maximum similarity: {np.max(similarity_matrix):.3f}")

    report.append("")
    report.append("Sequence Statistics:")
    seq_lengths = [len(seq) for seq in sequences.values()]
    report.append(f"  Average length: {np.mean(seq_lengths):.1f}")
    report.append(f"  Min length: {min(seq_lengths)}")
    report.append(f"  Max length: {max(seq_lengths)}")

    # Save report
    report_file = f"asmc_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
    with open(report_file, 'w') as f:
        f.write("\n".join(report))

    # Print report
    for line in report:
        print(line)

    print(f"\nReport saved to {report_file}")

    return report_file

def main():
    print("="*60)
    print("ASMC Analysis with Visualization")
    print("="*60)

    # Create test data
    test_dir = create_test_data()

    # Run analysis
    seq_names, similarity_matrix, sequences = run_simple_analysis()

    # Perform clustering
    clusters, linkage_matrix, cluster_groups = create_clustering_results(seq_names, similarity_matrix)

    # Create visualizations
    viz_file = visualize_results(seq_names, similarity_matrix, clusters, linkage_matrix, sequences)

    # Create summary report
    report_file = create_summary_report(seq_names, similarity_matrix, cluster_groups, sequences)

    print("\n" + "="*60)
    print("Analysis Complete!")
    print("="*60)
    print(f"Visualization: {viz_file}")
    print(f"Report: {report_file}")
    print("="*60)

if __name__ == "__main__":
    main()