# AlphaFold Workflow for ASMC

Tools for integrating AlphaFold-predicted structures with ASMC active site analysis.

## Scripts

- **run_asmc_with_alphafold.py** - Run ASMC using AlphaFold models as input
- **create_models_txt_for_alphafold.py** - Generate models.txt from AlphaFold output directory

## Overview

AlphaFold provides high-quality protein structure predictions. This workflow allows you to:
1. Use AlphaFold models as reference structures
2. Skip homology modeling step (ASMC uses AlphaFold models directly)
3. Perform pocket detection and clustering on AI-predicted structures

## Quick Start

### Step 1: Prepare AlphaFold Models

```bash
# Option A: You already have AlphaFold models
ls alphafold_models/
# AF-P12345-F1-model_v4.pdb
# AF-Q67890-F1-model_v4.pdb
# ...

# Option B: Download from AlphaFold DB
# Visit: https://alphafold.ebi.ac.uk/
```

### Step 2: Create models.txt

```bash
cd examples/alphafold_workflow

python create_models_txt_for_alphafold.py \
    --alphafold-dir /path/to/alphafold_models/ \
    --output ../../data/alphafold_models.txt
```

Output format:
```
/path/to/alphafold_models/AF-P12345-F1-model_v4.pdb
/path/to/alphafold_models/AF-Q67890-F1-model_v4.pdb
```

### Step 3: Run ASMC with AlphaFold Models

```bash
python run_asmc_with_alphafold.py \
    --models ../../data/alphafold_models.txt \
    --references ../../data/references_alphafold.txt \
    --pocket ../../data/pocket_alphafold.txt \
    --output ../../output/alphafold_asmc_results/
```

## Workflow Details

### Standard ASMC Workflow
```
Sequences → MODELLER → Structures → P2RANK → Alignment → Clustering
```

### AlphaFold Workflow (Skip MODELLER)
```
AlphaFold Models → P2RANK → Alignment → Clustering
```

## Advantages

1. **Skip homology modeling** - AlphaFold provides high-quality structures
2. **Faster** - No need to run MODELLER
3. **High confidence** - AlphaFold pLDDT scores indicate reliability
4. **Scalable** - Process thousands of AlphaFold predictions

## Important Notes

### Chain Selection
AlphaFold models typically contain a single chain. ASMC defaults to chain A.

```bash
# Specify chain if needed
asmc run --models models.txt --chain A --references refs.txt
```

### Pocket Detection
P2RANK works well with AlphaFold structures:
- High pLDDT regions (>70) are reliable for pocket prediction
- Low pLDDT regions (<50) may have unreliable pockets

### Reference Structure
Use experimental structures (PDB) as references when available:
```
# references_alphafold.txt
/path/to/experimental/3RFV.pdb
```

## Example: UDH with AlphaFold

```bash
# 1. Get AlphaFold models for UDH homologs
# (Download from AlphaFold DB or run AlphaFold locally)

# 2. Create models list
python create_models_txt_for_alphafold.py \
    --alphafold-dir udh_alphafold_models/ \
    --output ../../data/alphafold_models.txt

# 3. Use AtUdh (3RFV) as reference
echo "../../data/pdb_structures/AtUdh_3rfv.pdb" > ../../data/references_alphafold.txt

# 4. Define active sites (or use P2RANK)
cp ../../data/pdb_structures/atudh_active_site.txt ../../data/pocket_alphafold.txt

# 5. Run ASMC
python run_asmc_with_alphafold.py \
    --models ../../data/alphafold_models.txt \
    --references ../../data/references_alphafold.txt \
    --pocket ../../data/pocket_alphafold.txt \
    --output ../../output/udh_alphafold_results/
```

## Output

Same as standard ASMC:
- Cluster assignments
- Active site alignments
- Diversity analysis
- Visualizations

## Troubleshooting

### Issue: Low clustering quality
- **Solution**: Filter models by pLDDT score (>70 recommended)

### Issue: P2RANK fails on AlphaFold models
- **Solution**: Use pre-defined pocket file from experimental structure

### Issue: Large RMSD values
- **Solution**: AlphaFold models may have different domain orientations; use local alignment

## References

- AlphaFold Database: https://alphafold.ebi.ac.uk/
- AlphaFold paper: Jumper et al. (2021) Nature 596:583-589
- ASMC: https://github.com/labgem/ASMC
