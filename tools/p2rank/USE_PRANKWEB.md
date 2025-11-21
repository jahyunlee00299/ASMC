# Using PrankWeb (Online P2RANK) for AtUdh Active Site Prediction

## Why Use PrankWeb?

- No local installation needed
- Same P2RANK algorithm
- Free to use
- Results in minutes

## Step-by-Step Guide

### 1. Go to PrankWeb
https://prankweb.cz/

### 2. Upload Structure
- Click "Analyze Structure"
- Upload file: `C:\Users\Jahyun\PycharmProjects\ASMC\pdb_structures\pdb3rfv_chainA.pdb`
- Or enter PDB ID: `3RFV` (chain A)

### 3. Run Prediction
- Click "Predict Binding Sites"
- Wait ~1-2 minutes for results

### 4. Download Results
- Download CSV file with predicted pocket residues
- Save as: `p2rank_3rfv_results.csv`

### 5. Compare with Our Known Active Sites

Then run:
```bash
cd pdb_structures
python compare_p2rank_results.py p2rank_3rfv_results.csv
```

## What to Expect

P2RANK typically predicts:
- **Top 1-3 pockets** per protein
- **Pocket score** (0-100, higher = more confident)
- **Residue-level scores** for each position

### Our Known Active Sites (22 residues):
```
10, 12, 34, 36, 58, 80, 102, 104, 106, 135, 137, 139,
161, 163, 165, 187, 189, 211, 233, 235, 257, 259
```

### Expected Performance:
- **F1-score: 50-70%** (literature benchmark)
- Much better than our geometric heuristics (7-19%)
- May identify pockets we didn't annotate

## Manual Comparison (Without Script)

1. Look at PrankWeb visualization
2. Check which predicted pocket overlaps with:
   - NAD+ binding site (Rossmann fold)
   - UDP-glucose binding site
   - Catalytic residues

3. Cross-reference predicted residues with our list

## Alternative: Download P2RANK Results

PrankWeb provides:
- **3D visualization** (interactive)
- **Residue list** for each pocket
- **Confidence scores**
- **Downloadable PDB** with pocket annotations

This is the BEST way to validate our active site definitions!
