# P2RANK Integration for ASMC

This directory contains tools and documentation for integrating P2RANK pocket prediction with ASMC workflows.

## Files

- **INSTALL_P2RANK.md** - Complete installation guide for P2RANK on Windows
- **USE_PRANKWEB.md** - Guide for using PrankWeb (online P2RANK service)
- **compare_p2rank_results.py** - Script to compare P2RANK predictions with known active sites

## Quick Start

### Option 1: Use PrankWeb (Easiest)

1. Visit https://prankweb.cz/
2. Upload your PDB file
3. Download results
4. Compare with known sites:
   ```bash
   python compare_p2rank_results.py results.csv
   ```

### Option 2: Install P2RANK Locally

See `INSTALL_P2RANK.md` for detailed instructions.

## P2RANK Output Format

P2RANK generates a CSV file with predicted pockets:

```csv
name,rank,score,probability,sas_points,surf_atoms,center_x,center_y,center_z,residue_ids,surf_atom_ids
pocket1,1,18.3,0.95,156,89,12.5,34.2,56.7,A_10 A_12 A_34,1 5 12 ...
```

Column 9 (`residue_ids`) contains: `{chain}_{residue_number}` format

## Integration with ASMC

ASMC can automatically run P2RANK or use pre-computed pocket files:

```bash
# Auto-run P2RANK (requires prank in PATH)
asmc run -r references.txt -s sequences.fasta -o output/

# Use pre-computed pocket file
asmc run -r references.txt -s sequences.fasta -p pocket.txt -o output/
```

### Pocket File Format

```
pdb_file[TAB]chain[TAB]residue_numbers
AtUdh_3rfv.pdb	A	10,12,34,36,58,80,102,104,106,135,137,139
```

## Workflow: P2RANK → ASMC

1. **Run P2RANK on reference structure**
   ```bash
   prank predict reference.pdb -o p2rank_output/
   ```

2. **Extract pocket residues**
   ```python
   # ASMC does this automatically via asmc.extract_pocket()
   # Or manually parse predictions.csv
   ```

3. **Create pocket.txt**
   ```bash
   # Format: pdb[TAB]chain[TAB]residues
   echo "reference.pdb	A	10,12,34,36,58" > pocket.txt
   ```

4. **Run ASMC with pocket file**
   ```bash
   asmc run -r references.txt -s sequences.fasta -p pocket.txt -o output/
   ```

## Example: UDH Active Site Detection

See `../examples/udh_analysis/` for complete workflow using P2RANK to identify UDH active sites.

## References

- P2RANK GitHub: https://github.com/rdk/p2rank
- P2RANK Paper: Krivák & Hoksza (2018) J Cheminform 10:39
- PrankWeb: https://prankweb.cz/
