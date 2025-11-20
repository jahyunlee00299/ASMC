# PDB Structures

This directory contains protein structure files used for ASMC analysis.

## Contents

- **AtUdh_3rfv.pdb** - Arabidopsis thaliana UDP-glucose dehydrogenase (reference)
- **pdb3rfv_chainA.pdb** - Chain A only (for pocket prediction)  
- **atudh_active_site.txt** - Active site residue definitions for AtUdh
- **results/** - Active site detection results from structure alignment

## Reference Structure: AtUdh (3RFV)

- **Organism**: Arabidopsis thaliana
- **Enzyme**: UDP-glucose dehydrogenase  
- **Function**: Oxidizes UDP-glucose to UDP-glucuronic acid
- **Chain A**: 266 residues
- **Active Sites**: 22 residues (NAD+ binding + substrate binding)

### Active Site Residues

From `atudh_active_site.txt`:
```
10, 12, 34, 36, 58, 80, 102, 104, 106, 135, 137, 139,
161, 163, 165, 187, 189, 211, 233, 235, 257, 259
```

Regions:
- NAD+ binding (Rossmann fold)
- UDP-glucose substrate binding  
- Catalytic residues

## Using These Structures

### With ASMC
```bash
# Create reference file
echo "$(pwd)/AtUdh_3rfv.pdb" > ../references.txt

# Use pocket definition
cp atudh_active_site.txt ../udh_pocket.txt

# Run ASMC
asmc run -r ../references.txt -s ../sequences.fasta -p ../udh_pocket.txt
```

### With Active Site Detection
```bash
cd ../../examples/active_site_detection

python find_udh_active_sites.py \
    --reference ../../data/pdb_structures/AtUdh_3rfv.pdb \
    --reference-sites ../../data/pdb_structures/atudh_active_site.txt \
    --target your_protein.pdb
```

## Adding New Structures

1. Download PDB file
2. Place in this directory
3. Add to reference list if needed
4. Define active sites (manually or with P2RANK)
