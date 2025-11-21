#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
UDH (Uronate Dehydrogenase) 클러스터링 분석
1181개의 UDH 서열에 대한 ASMC 분석
"""

import subprocess
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
import seaborn as sns

def read_fasta(fasta_file):
    """FASTA 파일 읽기"""
    sequences = {}
    current_name = ""
    current_seq = ""

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_name:
                    sequences[current_name] = current_seq
                # Extract organism and protein name
                parts = line[1:].split('|')
                if len(parts) >= 2:
                    current_name = f"{parts[1]}_{parts[0]}"
                else:
                    current_name = line[1:50]  # First 50 chars
                current_seq = ""
            else:
                current_seq += line
        if current_name:
            sequences[current_name] = current_seq

    return sequences

def calculate_identity_matrix(sequences, sample_size=None):
    """서열 유사도 매트릭스 계산"""
    seq_names = list(sequences.keys())

    # 샘플링 (옵션)
    if sample_size and len(seq_names) > sample_size:
        print(f"Sampling {sample_size} sequences from {len(seq_names)}...")
        np.random.seed(42)
        indices = np.random.choice(len(seq_names), sample_size, replace=False)
        seq_names = [seq_names[i] for i in sorted(indices)]

    n = len(seq_names)
    print(f"Calculating identity matrix for {n} sequences...")

    identity_matrix = np.zeros((n, n))

    for i in range(n):
        if i % 50 == 0:
            print(f"  Progress: {i}/{n}")
        for j in range(i, n):
            seq1 = sequences[seq_names[i]]
            seq2 = sequences[seq_names[j]]

            # Pairwise identity
            min_len = min(len(seq1), len(seq2))
            if min_len > 0:
                identical = sum(1 for k in range(min_len) if seq1[k] == seq2[k])
                identity = identical / min_len
            else:
                identity = 0

            identity_matrix[i, j] = identity
            identity_matrix[j, i] = identity

    print(f"Identity matrix calculated!")
    return seq_names, identity_matrix

def perform_clustering(seq_names, identity_matrix, method='average', threshold=0.3):
    """계층적 클러스터링 수행"""
    print(f"\nPerforming hierarchical clustering...")
    print(f"  Method: {method}")
    print(f"  Threshold: {threshold}")

    # Convert identity to distance
    distance_matrix = 1 - identity_matrix

    # Condensed distance matrix for linkage
    condensed_dist = squareform(distance_matrix)

    # Hierarchical clustering
    linkage_matrix = linkage(condensed_dist, method=method)

    # Form clusters
    clusters = fcluster(linkage_matrix, threshold, criterion='distance')

    # Group by cluster
    cluster_groups = {}
    for i, cluster_id in enumerate(clusters):
        if cluster_id not in cluster_groups:
            cluster_groups[cluster_id] = []
        cluster_groups[cluster_id].append(seq_names[i])

    print(f"\nClustering complete!")
    print(f"  Number of clusters: {len(cluster_groups)}")
    print(f"  Cluster sizes: min={min(len(v) for v in cluster_groups.values())}, "
          f"max={max(len(v) for v in cluster_groups.values())}, "
          f"mean={np.mean([len(v) for v in cluster_groups.values()]):.1f}")

    return clusters, linkage_matrix, cluster_groups

def visualize_results(seq_names, identity_matrix, clusters, linkage_matrix, cluster_groups, sequences):
    """결과 시각화"""
    print("\nCreating visualizations...")

    # 너무 많은 서열은 시각화가 어려우므로 샘플링
    n_seqs = len(seq_names)
    if n_seqs > 100:
        print(f"  Sampling 100 sequences for heatmap from {n_seqs}...")
        sample_indices = np.random.choice(n_seqs, min(100, n_seqs), replace=False)
        sample_names = [seq_names[i] for i in sorted(sample_indices)]
        sample_matrix = identity_matrix[np.ix_(sample_indices, sample_indices)]
    else:
        sample_names = seq_names
        sample_matrix = identity_matrix

    fig = plt.figure(figsize=(20, 14))

    # 1. Similarity heatmap (sampled)
    ax1 = plt.subplot(2, 3, 1)
    im = ax1.imshow(sample_matrix, cmap='RdYlGn', vmin=0, vmax=1, aspect='auto')
    ax1.set_title(f'Sequence Identity Heatmap\n(Sample of {len(sample_names)} sequences)', fontsize=12)
    ax1.set_xlabel('Sequence Index')
    ax1.set_ylabel('Sequence Index')
    plt.colorbar(im, ax=ax1, label='Identity')

    # 2. Dendrogram (full or sampled)
    ax2 = plt.subplot(2, 3, 2)
    if n_seqs <= 50:
        dendrogram(linkage_matrix, ax=ax2, labels=[name[:20] for name in seq_names],
                   leaf_rotation=90, leaf_font_size=6)
        ax2.set_title('Hierarchical Clustering Dendrogram', fontsize=12)
    else:
        dendrogram(linkage_matrix, ax=ax2, no_labels=True)
        ax2.set_title(f'Hierarchical Clustering Dendrogram\n({n_seqs} sequences)', fontsize=12)
    ax2.set_xlabel('Sequence')
    ax2.set_ylabel('Distance')

    # 3. Cluster size distribution
    ax3 = plt.subplot(2, 3, 3)
    cluster_sizes = [len(members) for members in cluster_groups.values()]
    cluster_ids = list(cluster_groups.keys())

    if len(cluster_ids) <= 50:
        colors = plt.cm.tab20(np.linspace(0, 1, len(cluster_ids)))
        ax3.bar(cluster_ids, cluster_sizes, color=colors, edgecolor='black')
        ax3.set_xlabel('Cluster ID')
    else:
        ax3.hist(cluster_sizes, bins=30, edgecolor='black', alpha=0.7)
        ax3.set_xlabel('Cluster Size')
    ax3.set_ylabel('Count' if len(cluster_ids) > 50 else 'Number of Sequences')
    ax3.set_title(f'Cluster Size Distribution\n({len(cluster_groups)} clusters)', fontsize=12)

    # 4. Sequence length distribution
    ax4 = plt.subplot(2, 3, 4)
    seq_lengths = [len(seq) for seq in sequences.values()]
    ax4.hist(seq_lengths, bins=50, edgecolor='black', alpha=0.7, color='steelblue')
    ax4.set_xlabel('Sequence Length (aa)')
    ax4.set_ylabel('Frequency')
    ax4.set_title(f'Sequence Length Distribution\n(mean: {np.mean(seq_lengths):.1f} aa)', fontsize=12)
    ax4.axvline(np.mean(seq_lengths), color='red', linestyle='--',
                label=f'Mean: {np.mean(seq_lengths):.1f}')
    ax4.axvline(np.median(seq_lengths), color='orange', linestyle='--',
                label=f'Median: {np.median(seq_lengths):.1f}')
    ax4.legend()

    # 5. Identity distribution
    ax5 = plt.subplot(2, 3, 5)
    # Upper triangle only (excluding diagonal)
    upper_tri_indices = np.triu_indices_from(identity_matrix, k=1)
    identity_values = identity_matrix[upper_tri_indices]
    ax5.hist(identity_values, bins=50, edgecolor='black', alpha=0.7, color='green')
    ax5.set_xlabel('Pairwise Identity')
    ax5.set_ylabel('Frequency')
    ax5.set_title(f'Pairwise Identity Distribution\n(mean: {np.mean(identity_values):.3f})', fontsize=12)
    ax5.axvline(np.mean(identity_values), color='red', linestyle='--',
                label=f'Mean: {np.mean(identity_values):.3f}')
    ax5.legend()

    # 6. Top clusters composition
    ax6 = plt.subplot(2, 3, 6)
    # Show top 10 largest clusters
    sorted_clusters = sorted(cluster_groups.items(), key=lambda x: len(x[1]), reverse=True)
    top_n = min(15, len(sorted_clusters))
    top_cluster_ids = [f"C{c_id}" for c_id, _ in sorted_clusters[:top_n]]
    top_cluster_sizes = [len(members) for _, members in sorted_clusters[:top_n]]

    colors_top = plt.cm.Set3(np.linspace(0, 1, top_n))
    bars = ax6.barh(range(top_n), top_cluster_sizes, color=colors_top, edgecolor='black')
    ax6.set_yticks(range(top_n))
    ax6.set_yticklabels(top_cluster_ids)
    ax6.set_xlabel('Number of Sequences')
    ax6.set_title(f'Top {top_n} Largest Clusters', fontsize=12)
    ax6.invert_yaxis()

    # Add values on bars
    for i, (bar, size) in enumerate(zip(bars, top_cluster_sizes)):
        ax6.text(bar.get_width() + max(top_cluster_sizes)*0.01, bar.get_y() + bar.get_height()/2,
                f'{size}', va='center', fontsize=9)

    plt.suptitle('UDH (Uronate Dehydrogenase) Clustering Analysis',
                 fontsize=16, fontweight='bold', y=0.995)
    plt.tight_layout()

    # Save figure
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_file = f"udh_clustering_{timestamp}.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Visualization saved: {output_file}")

    return output_file

def create_report(seq_names, identity_matrix, cluster_groups, sequences):
    """분석 리포트 생성"""
    print("\nGenerating analysis report...")

    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    report_file = f"udh_analysis_report_{timestamp}.txt"

    report = []
    report.append("=" * 70)
    report.append("UDH (Uronate Dehydrogenase) Clustering Analysis Report")
    report.append("=" * 70)
    report.append(f"\nAnalysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report.append(f"Total sequences analyzed: {len(sequences)}")
    report.append(f"Sequences used in clustering: {len(seq_names)}")
    report.append(f"Number of clusters: {len(cluster_groups)}")

    # Sequence statistics
    seq_lengths = [len(seq) for seq in sequences.values()]
    report.append("\n" + "-" * 70)
    report.append("Sequence Statistics:")
    report.append(f"  Length range: {min(seq_lengths)} - {max(seq_lengths)} aa")
    report.append(f"  Mean length: {np.mean(seq_lengths):.1f} aa")
    report.append(f"  Median length: {np.median(seq_lengths):.1f} aa")

    # Identity statistics
    upper_tri = np.triu_indices_from(identity_matrix, k=1)
    identity_values = identity_matrix[upper_tri]
    report.append("\n" + "-" * 70)
    report.append("Pairwise Identity Statistics:")
    report.append(f"  Mean identity: {np.mean(identity_values):.3f}")
    report.append(f"  Median identity: {np.median(identity_values):.3f}")
    report.append(f"  Min identity: {np.min(identity_values):.3f}")
    report.append(f"  Max identity: {np.max(identity_values):.3f}")

    # Cluster statistics
    cluster_sizes = [len(members) for members in cluster_groups.values()]
    report.append("\n" + "-" * 70)
    report.append("Cluster Statistics:")
    report.append(f"  Number of clusters: {len(cluster_groups)}")
    report.append(f"  Cluster size range: {min(cluster_sizes)} - {max(cluster_sizes)}")
    report.append(f"  Mean cluster size: {np.mean(cluster_sizes):.1f}")
    report.append(f"  Median cluster size: {np.median(cluster_sizes):.1f}")

    # Top 10 largest clusters
    sorted_clusters = sorted(cluster_groups.items(), key=lambda x: len(x[1]), reverse=True)
    report.append("\n" + "-" * 70)
    report.append("Top 10 Largest Clusters:")
    for i, (cluster_id, members) in enumerate(sorted_clusters[:10], 1):
        report.append(f"\n  Cluster {cluster_id}: {len(members)} members")
        # Show first 5 members
        for j, member in enumerate(members[:5], 1):
            report.append(f"    {j}. {member[:60]}")
        if len(members) > 5:
            report.append(f"    ... and {len(members) - 5} more")

    # Save report
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("\n".join(report))

    # Print report
    for line in report:
        print(line)

    print(f"\nReport saved: {report_file}")
    return report_file

def save_cluster_results(cluster_groups, sequences, seq_names):
    """클러스터 결과를 파일로 저장"""
    print("\nSaving cluster results...")

    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

    # Save cluster assignments
    cluster_file = f"udh_clusters_{timestamp}.tsv"
    with open(cluster_file, 'w') as f:
        f.write("Sequence\tCluster_ID\tSequence_Length\n")
        for cluster_id, members in cluster_groups.items():
            for member in members:
                seq_len = len(sequences.get(member, ""))
                f.write(f"{member}\t{cluster_id}\t{seq_len}\n")
    print(f"  Cluster assignments saved: {cluster_file}")

    # Save cluster FASTAs (top 5 largest clusters)
    sorted_clusters = sorted(cluster_groups.items(), key=lambda x: len(x[1]), reverse=True)
    for i, (cluster_id, members) in enumerate(sorted_clusters[:5], 1):
        fasta_file = f"udh_cluster_{cluster_id}_{timestamp}.fasta"
        with open(fasta_file, 'w') as f:
            for member in members:
                if member in sequences:
                    f.write(f">{member}\n{sequences[member]}\n")
        print(f"  Cluster {cluster_id} FASTA saved: {fasta_file}")

    return cluster_file

def main():
    """Main analysis function"""
    print("=" * 70)
    print("UDH (Uronate Dehydrogenase) Clustering Analysis")
    print("=" * 70)

    # Paths
    fasta_file = Path("test_data/UDHs_filtered_std2.5.fasta")
    pdb_file = Path("test_data/AtUdh_pdb3rfv_chainA.pdb")

    # Check files exist
    if not fasta_file.exists():
        print(f"Error: FASTA file not found: {fasta_file}")
        return
    if not pdb_file.exists():
        print(f"Warning: PDB file not found: {pdb_file}")

    # Read sequences
    print(f"\nReading sequences from {fasta_file}...")
    sequences = read_fasta(fasta_file)
    print(f"  Total sequences: {len(sequences)}")

    # Calculate identity matrix (with optional sampling for large datasets)
    # For 1181 sequences, we'll sample to make it manageable
    sample_size = 500 if len(sequences) > 500 else None
    seq_names, identity_matrix = calculate_identity_matrix(sequences, sample_size=sample_size)

    # Perform clustering
    clusters, linkage_matrix, cluster_groups = perform_clustering(
        seq_names, identity_matrix, method='average', threshold=0.5
    )

    # Create visualizations
    viz_file = visualize_results(seq_names, identity_matrix, clusters,
                                 linkage_matrix, cluster_groups, sequences)

    # Generate report
    report_file = create_report(seq_names, identity_matrix, cluster_groups, sequences)

    # Save cluster results
    cluster_file = save_cluster_results(cluster_groups, sequences, seq_names)

    print("\n" + "=" * 70)
    print("Analysis Complete!")
    print("=" * 70)
    print(f"Visualization: {viz_file}")
    print(f"Report: {report_file}")
    print(f"Cluster data: {cluster_file}")
    print("=" * 70)

if __name__ == "__main__":
    main()