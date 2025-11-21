# ASMC Data Directory

This directory contains input data files for ASMC analysis.

## Structure

```
data/
├── pdb_structures/          # PDB structure files and active site definitions
├── *.fasta                  # Sequence files
├── *_models.txt            # Model file lists
├── *_references.txt        # Reference structure lists
└── *_pocket.txt            # Active site pocket definitions
```

## File Types

### PDB Structures (`pdb_structures/`)
- Experimental structures from PDB
- AlphaFold predicted models
- Active site definition files (`.txt`)

### Sequences (`.fasta`)
- `sequences.fasta` - General test sequences
- `udh_active_sites.fasta` - UDH active site regions
- `udh_substrate_sites.fasta` - UDH substrate binding sites

### Reference Lists (`*_references.txt`)
Format: One PDB file path per line
```
/path/to/reference1.pdb
/path/to/reference2.pdb
```

### Pocket Files (`*_pocket.txt`)
Format: `pdb_file[TAB]chain[TAB]residue_numbers`
```
AtUdh_3rfv.pdb	A	10,12,34,36,58,80,102
```

## Usage with ASMC

```bash
# Run ASMC with data from this directory
asmc run \
    -r data/references.txt \
    -s data/sequences.fasta \
    -p data/pocket.txt \
    -o output/
```

## Adding New Data

1. Place PDB files in `pdb_structures/`
2. Create reference list: `ls pdb_structures/*.pdb > new_references.txt`
3. Define pockets (manually or use P2RANK)
4. Add sequences in FASTA format
