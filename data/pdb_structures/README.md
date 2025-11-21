# PDB Structures

This directory contains protein structure files used for ASMC analysis.

## Contents

- **AtUdh_3rfv.pdb** - Arabidopsis thaliana UDP-glucose dehydrogenase (reference)
- **pdb3rfv_chainA.pdb** - Chain A only (for pocket prediction)

### Active Site Definition Files
- **atudh_active_site.txt** - NAD+ binding site (Direct+Proximal, 31 residues)
- **atudh_nad_site.txt** - NAD+ binding site only (31 residues)
- **atudh_substrate_site.txt** - Substrate (galacturonate) binding site (12 residues)
- **atudh_combined_3_5.5A.txt** - Combined NAD+Substrate (3.0-5.5A cutoff, 40 residues)

- **results/** - Active site detection results from structure alignment

## Reference Structure: AtUdh (3RFV)

- **Organism**: Arabidopsis thaliana
- **Enzyme**: UDP-glucose dehydrogenase
- **Function**: Oxidizes UDP-glucose to UDP-glucuronic acid
- **Chain A**: 265 residues (PDB residue 2-266)
- **Ligands**: NAI (NADH), 15L (galacturonate analog)

## PDB Numbering Offset

**PDB 3RFV는 residue 2번부터 시작합니다.**

```
PDB_num = alignment_pos + 2
alignment_pos = PDB_num - 2
```

**예시:**
- Y136 (PDB) = Y134 (alignment) - 핵심 촉매 잔기
- NAI와 2.54A, 15L과 2.53A 거리

### AlphaFold 구조와의 차이

AlphaFold 구조는 residue 1번부터 시작하므로:
- 3RFV Y136 → AlphaFold ~Y134-140 (구조에 따라 다름)
- 서열 정렬 vs 구조적 거리로 대응 잔기가 다를 수 있음
- 예: A0A0N8HGG8에서 서열정렬=Y139, 구조거리=Y140

## Active Site Definitions (2025-11-21 Updated)

### 1. NAD+ Binding Site (atudh_nad_site.txt)
리간드 거리 기준 (NAI ligand):
- **Direct (≤3.5A)**: 21 residues
- **Proximal (3.5-5.0A)**: 10 residues
- **Total**: 31 residues

### 2. Substrate Binding Site (atudh_substrate_site.txt)
리간드 거리 기준 (15L galacturonate):
- **Direct (≤3.5A)**: S75, S111, N112, H113, **Y136**, S165, R174, F258
- **Proximal (3.5-5.0A)**: V76, I163, G164, M175
- **Total**: 12 residues

### 3. Combined Site (atudh_combined_3_5.5A.txt)
NAI + 15L, 3.0-5.5A cutoff:
- **Total**: 40 residues

## Using These Structures

### With ASMC
```bash
# Create reference file
echo "$(pwd)/AtUdh_3rfv.pdb" > ../references.txt

# Use pocket definition
cp atudh_substrate_site.txt ../udh_pocket.txt

# Run ASMC
asmc run -r ../references.txt -s ../sequences.fasta -p ../udh_pocket.txt
```

### With Active Site Detection
```bash
cd ../../examples/active_site_detection

python find_udh_active_sites.py \
    --reference ../../data/pdb_structures/AtUdh_3rfv.pdb \
    --reference-sites ../../data/pdb_structures/atudh_substrate_site.txt \
    --target your_protein.pdb
```

## Adding New Structures

1. Download PDB file
2. Place in this directory
3. Add to reference list if needed
4. Define active sites (manually or with P2RANK)
5. **Check residue numbering offset!**
