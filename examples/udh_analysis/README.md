# UDH (Uronate Dehydrogenase) Analysis Examples

This directory contains scripts for analyzing uronate dehydrogenase (UDH) enzyme family active sites and substrate specificity.

## Scripts

### Core Analysis
- **run_udh_asmc.py** - Run ASMC on UDH enzyme family
- **analyze_udh.py** - Analyze UDH clustering and diversity
- **extract_udh_active_sites.py** - Extract active site regions from sequences
- **extract_udh_substrate_sites.py** - Extract substrate binding site regions

### Result Analysis
- **analyze_asmc_results.py** - General ASMC result analysis
- **analyze_substrate_results.py** - Substrate-specific clustering analysis
- **analyze_substrate_variants.py** - Analyze substrate specificity variants

## Quick Start

### 1. Basic UDH Analysis

```bash
cd examples/udh_analysis

# Run ASMC on UDH sequences with known active sites
python run_udh_asmc.py \
    --references ../../data/udh_references.txt \
    --sequences ../../data/udh_active_sites.fasta \
    --pocket ../../data/udh_pocket.txt \
    --output ../../output/udh_results
```

### 2. Extract Active Site Sequences

```bash
# Extract active site regions from full-length sequences
python extract_udh_active_sites.py \
    --input ../../data/sequences.fasta \
    --pocket ../../data/udh_pocket.txt \
    --output ../../output/udh_active_sites.fasta
```

### 3. Analyze Results

```bash
# Analyze clustering results
python analyze_udh.py \
    --input ../../output/udh_results/ \
    --output ../../output/udh_analysis_report.txt
```

## Data Files

Required data files (in `../../data/`):
- `udh_references.txt` - Reference PDB structures
- `udh_pocket.txt` - Active site residue definitions
- `udh_active_sites.fasta` - UDH sequences with active site regions
- `udh_substrate_sites.fasta` - UDH sequences with substrate binding sites

## Reference Structure

The reference structure for UDH analysis is:
- **PDB: 3RFV** - Arabidopsis thaliana UDP-glucose dehydrogenase (AtUdh)
- **Chain A** - 266 residues
- **Active sites**: 22 residues involved in NAD+ binding and substrate coordination

See `../../data/pdb_structures/atudh_active_site.txt` for details.

## Workflow

```
1. Prepare Data
   ├── Reference PDB (AtUdh 3RFV)
   ├── Active site definition
   └── UDH sequence collection

2. Extract Active Sites
   └── extract_udh_active_sites.py

3. Run ASMC
   └── run_udh_asmc.py
       ├── Pocket detection (P2RANK or predefined)
       ├── Structural alignment (USalign)
       └── Clustering (DBSCAN)

4. Analyze Results
   ├── analyze_udh.py
   ├── analyze_substrate_results.py
   └── analyze_substrate_variants.py

5. Visualize
   └── Clustering plots
   └── Diversity analysis
```

## Expected Outputs

- Cluster assignments for UDH sequences
- Active site diversity plots
- Substrate specificity groupings
- FASTA files for each cluster

## Notes

- UDH family shows high structural diversity (RMSD up to 18Å between homologs)
- Substrate binding sites are more conserved than overall structure
- NAD+ binding (Rossmann fold) is the most conserved region

## References

- Campbell et al. (1997) Plant Cell - AtUdh characterization
- P2RANK for pocket detection: https://github.com/rdk/p2rank
- ASMC workflow: https://github.com/labgem/ASMC
