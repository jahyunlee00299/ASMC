# AtUdh (3RFV Chain A) Active Site - Detailed Analysis

## Overview
- **Protein**: UDP-glucose dehydrogenase from Arabidopsis thaliana
- **PDB ID**: 3RFV
- **Chain**: A
- **Total Active Site Residues**: 23

## Complete Active Site Residue List

| Position | Residue | Type | CA Coordinates (Å) | Atoms | Notes |
|----------|---------|------|-------------------|-------|-------|
| 10 | GLY | Small/Flexible | (67.75, 40.09, 12.40) | 4 | NAD+ binding region |
| 12 | ALA | Small/Hydrophobic | (66.05, 37.71, 17.12) | 5 | NAD+ binding region |
| 34 | ASP | Acidic | (66.74, 35.53, 11.55) | 8 | Substrate coordination |
| 36 | SER | Polar | (68.50, 32.18, 15.00) | 6 | NAD+ binding |
| 58 | VAL | Hydrophobic | (66.53, 35.03, 1.10) | 7 | Structural support |
| 80 | PHE | Aromatic | (86.34, 38.10, 8.03) | 11 | Hydrophobic pocket |
| 102 | GLY | Small/Flexible | (58.99, 44.53, -6.67) | 4 | Loop flexibility |
| 104 | PRO | Rigid | (60.28, 46.89, -1.62) | 7 | Structural rigidity |
| 106 | ILE | Hydrophobic | (64.66, 49.15, 2.42) | 8 | Hydrophobic pocket |
| 135 | LEU | Hydrophobic | (86.17, 44.72, 8.91) | 8 | Hydrophobic pocket |
| 137 | GLY | Small/Flexible | (82.71, 49.00, 9.31) | 4 | Loop flexibility |
| 139 | SER | Polar | (81.21, 45.63, 5.39) | 6 | Potential H-bonding |
| 161 | VAL | Hydrophobic | (70.94, 54.60, 8.41) | 7 | Structural support |
| 163 | ILE | Hydrophobic | (73.14, 54.03, 13.59) | 8 | Hydrophobic pocket |
| 165 | SER | Polar | (74.05, 53.75, 19.76) | 6 | Potential H-bonding |
| 187 | SER | Polar | (59.86, 57.58, 13.94) | 6 | NAD+ binding |
| 189 | ILE | Hydrophobic | (61.47, 54.25, 9.99) | 8 | Hydrophobic pocket |
| 211 | ALA | Small/Hydrophobic | (76.49, 66.89, 21.82) | 5 | Structural support |
| 233 | PHE | Aromatic | (71.43, 61.97, 30.79) | 11 | Substrate binding |
| 235 | ARG | Basic | (74.47, 63.70, 35.16) | 11 | Substrate coordination |
| 257 | THR | Polar | (87.01, 51.45, 21.07) | 7 | NAD+ binding |
| 259 | VAL | Hydrophobic | (89.50, 51.57, 16.56) | 7 | Structural support |
| 281 | - | - | - | - | **Missing from structure** |

## Residue Type Distribution

| Residue Type | Count | Percentage |
|--------------|-------|------------|
| Serine (SER) | 4 | 17.4% |
| Glycine (GLY) | 3 | 13.0% |
| Isoleucine (ILE) | 3 | 13.0% |
| Valine (VAL) | 3 | 13.0% |
| Phenylalanine (PHE) | 2 | 8.7% |
| Alanine (ALA) | 2 | 8.7% |
| Aspartate (ASP) | 1 | 4.3% |
| Proline (PRO) | 1 | 4.3% |
| Leucine (LEU) | 1 | 4.3% |
| Arginine (ARG) | 1 | 4.3% |
| Threonine (THR) | 1 | 4.3% |

## Functional Classification

### 1. NAD+ Binding Site (Rossmann Fold)
**Critical for cofactor binding and catalysis**

- **GLY 10**: Glycine-rich loop, flexible backbone for NAD+ phosphate binding
- **ALA 12**: Hydrophobic environment near nicotinamide ring
- **SER 36**: Hydroxyl group for H-bonding with NAD+ ribose/phosphate
- **SER 187**: Hydroxyl group for NAD+ coordination
- **THR 257**: Hydroxyl group for NAD+ binding
- **VAL 259**: Hydrophobic support near NAD+ binding pocket

**Key Feature**: Rossmann fold motif (GxGxxG/A) for nucleotide binding

### 2. UDP-Glucose Substrate Binding Site
**Substrate recognition and positioning**

- **ASP 34**: Negative charge for coordination with substrate
- **PHE 233**: Aromatic stacking with glucose ring
- **ARG 235**: Positive charge for UDP phosphate coordination
- **SER 165**: Potential H-bonding with glucose hydroxyl groups

### 3. Catalytic Residues
**Direct involvement in oxidation reaction**

- **ASP 34**: Proton abstraction/general base
- **ARG 235**: Stabilization of transition state
- **SER 139**: Potential involvement in catalytic mechanism

**Note**: The actual catalytic cysteine might be missing from this list (residue 281 is absent from the structure)

### 4. Hydrophobic Pocket
**Structural integrity and substrate recognition**

- **VAL 58**: Creates hydrophobic environment
- **PHE 80**: Aromatic ring for substrate stacking
- **ILE 106**: Hydrophobic side chain
- **LEU 135**: Hydrophobic pocket formation
- **VAL 161**: Structural support
- **ILE 163**: Hydrophobic environment
- **ILE 189**: Substrate binding pocket
- **ALA 211**: Small hydrophobic residue

### 5. Structural/Flexibility Elements
**Maintain active site geometry and dynamics**

- **GLY 102**: Loop flexibility near active site
- **GLY 137**: Flexible hinge region
- **PRO 104**: Rigid proline for structural constraint

## Spatial Organization

### Region 1: N-terminal Domain (Residues 10-80)
- NAD+ binding initiation
- Contains early Rossmann fold elements
- Substrate entry point

### Region 2: Central Core (Residues 102-165)
- Catalytic center
- Substrate and cofactor convergence
- High concentration of polar residues

### Region 3: C-terminal Domain (Residues 187-259)
- Additional NAD+ stabilization
- Substrate specificity determination
- Domain closure elements

## Chemical Properties Summary

| Property | Count | Residues |
|----------|-------|----------|
| **Polar** | 8 | SER(4), THR(1), ASP(1), ARG(1), GLY(1)* |
| **Hydrophobic** | 12 | VAL(3), ILE(3), ALA(2), LEU(1), PHE(2), PRO(1) |
| **Aromatic** | 2 | PHE(2) |
| **Charged (+)** | 1 | ARG(1) |
| **Charged (-)** | 1 | ASP(1) |
| **Small/Flexible** | 3 | GLY(3) |

*GLY counted separately due to unique structural role

## Key Observations

1. **High Serine Content (17.4%)**: Critical for NAD+ binding and potential catalytic roles
2. **Balanced Hydrophobic/Polar**: ~52% hydrophobic, ~35% polar - optimal for substrate binding
3. **Limited Charged Residues**: Only ASP and ARG - highly specific substrate coordination
4. **Glycine Flexibility**: 3 glycines provide conformational flexibility for induced fit
5. **Missing Residue 281**: Likely critical catalytic residue (possibly Cysteine) not resolved in structure

## Comparison with UDH Family

Based on the alignment results with other UDH structures:
- **DH35**: 18/23 residues conserved (78%)
- **DHP6**: 20/23 residues conserved (87%)
- **MATA**: 13/23 residues conserved (57%)

The conservation pattern suggests:
- NAD+ binding residues are highly conserved
- Substrate binding pocket shows more variability
- Structural elements (GLY, PRO) are well-conserved

## References

1. **Primary Structure**: PDB 3RFV
2. **Literature**: Campbell et al., Plant Cell, 1997
3. **Alignment Data**: Generated from find_udh_active_sites.py analysis
4. **Residue Definitions**: From atudh_active_site.txt

## Notes

- Residue 281 is defined in the active site list but not present in the PDB structure (likely disordered or missing density)
- All coordinates are from Chain A
- Active site definition may include both direct catalytic residues and structural support residues within ~4Å
- Some residues may serve multiple functional roles

---

**Analysis Date**: 2025-11-20
**Tool**: BioPython PDBParser
**Reference**: pdb3rfv_chainA.pdb
