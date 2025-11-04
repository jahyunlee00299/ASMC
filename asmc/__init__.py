from importlib.metadata import distribution

__version__ = distribution("asmc").version

# Import all functions and exceptions from asmc.py
from .asmc import (
    build_ds,
    extract_pocket,
    build_pocket_text,
    renumber_residues,
    extract_aligned_pos,
    build_multiple_alignment,
    search_active_site_in_msa,
    read_alignment,
    read_matrix,
    pairwise_score,
    dissimilarity,
    dbscan_clustering,
    formatting_output,
    build_fasta,
    build_logo,
    merge_logo,
    RenumberResiduesError
)

# Export all functions
__all__ = [
    "build_ds",
    "extract_pocket",
    "build_pocket_text",
    "renumber_residues",
    "extract_aligned_pos",
    "build_multiple_alignment",
    "search_active_site_in_msa",
    "read_alignment",
    "read_matrix",
    "pairwise_score",
    "dissimilarity",
    "dbscan_clustering",
    "formatting_output",
    "build_fasta",
    "build_logo",
    "merge_logo",
    "RenumberResiduesError"
]