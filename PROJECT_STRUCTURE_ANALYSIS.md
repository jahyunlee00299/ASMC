# ASMC Project Structure Analysis

## Current vs. Official GitHub Structure

### Official ASMC Structure (GitHub: labgem/ASMC)
```
ASMC/
├── asmc/                 # Main package
├── docs/                 # Documentation
├── tests/                # Tests
├── pyproject.toml        # Project config
├── requirements.txt      # Dependencies
├── env.yml              # Conda environment
├── README.md
├── LICENSE.txt
├── MANIFEST.in
└── VERSION
```

### Current Project Structure
```
ASMC/
├── asmc/                      # ✓ Main package (correct)
├── docs/                      # ✓ Documentation (correct)
├── tests/                     # ✓ Tests (correct)
├── pdb_structures/           # ✗ NEW: Not in official repo
│   ├── AtUdh_3rfv.pdb
│   ├── find_udh_active_sites.py
│   ├── predict_pockets_*.py
│   └── results/
├── pyproject.toml            # ✓ Correct
├── requirements.txt          # ✓ Correct
├── env.yml                   # ✓ Correct
├── README.md                 # ✓ Correct
├── LICENSE.txt               # ✓ Correct
├── MANIFEST.in               # ✓ Correct
├── VERSION                   # ✓ Correct
│
├── # Test/Analysis Scripts (should be in examples/ or tests/)
├── analyze_*.py              # ✗ Should move to examples/
├── extract_*.py              # ✗ Should move to examples/
├── run_*.py                  # ✗ Should move to examples/
├── 실행_예제.py               # ✗ Should move to examples/
│
├── # Data Files (should be in data/ or test_data/)
├── alphafold_models.txt      # ✗ Should move to data/
├── pocket_alphafold.txt      # ✗ Should move to data/
├── references*.txt           # ✗ Should move to data/
├── models.txt                # ✗ Should move to data/
├── sequences.fasta           # ✗ Should move to data/
├── udh_*.fasta               # ✗ Should move to data/
│
└── # Output Directories (generated, should be in .gitignore)
    ├── output_*/
    ├── udh_asmc_*/
    ├── 3rfv_superposition_results/
    └── asmc_alphafold_output/
```

## Issues Identified

### 1. **Root Directory Clutter**
- Too many analysis/test scripts in root
- Should be organized into subdirectories

### 2. **pdb_structures/ Directory**
- Not part of official ASMC
- Contains custom UDH active site detection tools
- Should be reorganized

### 3. **Missing Directories**
- `examples/` - For demonstration scripts
- `data/` - For input data files

### 4. **P2RANK Integration**
According to official ASMC:
- P2RANK should be installed separately
- ASMC calls P2RANK via subprocess
- Results parsed from P2RANK CSV output

## Recommended Structure

```
ASMC/
├── asmc/                      # Core package (don't touch)
├── docs/                      # Documentation
├── tests/                     # Unit tests
│
├── examples/                  # NEW: Example workflows
│   ├── udh_analysis/
│   │   ├── run_udh_asmc.py
│   │   ├── analyze_udh.py
│   │   ├── extract_udh_active_sites.py
│   │   └── README.md
│   ├── alphafold_workflow/
│   │   ├── run_asmc_with_alphafold.py
│   │   ├── create_models_txt_for_alphafold.py
│   │   └── README.md
│   └── active_site_detection/
│       ├── find_udh_active_sites.py  # Structure-based
│       ├── predict_pockets_*.py      # P2RANK alternatives
│       └── README.md
│
├── data/                      # NEW: Input data
│   ├── pdb_structures/
│   │   ├── AtUdh_3rfv.pdb
│   │   ├── pdb3rfv_chainA.pdb
│   │   ├── atudh_active_site.txt
│   │   └── README.md
│   ├── alphafold_models.txt
│   ├── pocket_alphafold.txt
│   ├── references_alphafold.txt
│   └── udh_active_sites.fasta
│
├── tools/                     # NEW: External tool integration
│   ├── p2rank/
│   │   ├── install_p2rank.sh
│   │   ├── INSTALL_P2RANK.md
│   │   └── USE_PRANKWEB.md
│   └── README.md
│
├── output/                    # Git-ignored output directory
│   └── .gitkeep
│
├── pyproject.toml
├── requirements.txt
├── env.yml
├── README.md
├── LICENSE.txt
├── MANIFEST.in
└── VERSION
```

## Key Principles from Official ASMC

### 1. P2RANK Usage
```python
# From asmc/run_asmc.py:80-108
def run_prank(ds, outdir):
    """Run p2rank via subprocess"""
    P2RANK = "prank"
    command = f"{P2RANK} predict {ds} -o {outdir}"
    result = subprocess.run(command.split(), capture_output=True)
    return result

# From asmc/asmc.py:72-123
def extract_pocket(outdir):
    """Extract pocket positions from P2RANK predictions.csv"""
    prediction = [f for f in outdir.iterdir() if f.match("*predictions*")][0]
    pred_arr = np.loadtxt(prediction, skiprows=1, delimiter=",", dtype=str)

    # Parse column 9 which contains: "A_10 A_12 A_34 ..."
    res_str = pred_arr[i][9]
    res_list = [(elem.split("_")[0], elem.split("_")[1])
                for elem in res_str.split()]
```

### 2. Pocket File Format
```
# pocket.txt format (from build_pocket_text)
pdb_file[TAB]chain[TAB]residue_numbers

# Example:
AtUdh_3rfv.pdb	A	10,12,34,36,58,80,102,104,106,135
```

### 3. Workflow Integration
```bash
# Official ASMC workflow
asmc run -r references.txt -s sequences.fasta -p pocket.txt -o output/

# Where:
# - references.txt: paths to PDB files
# - sequences.fasta: sequences to model
# - pocket.txt: active site definitions
# - If no -p: runs P2RANK automatically
```

## Action Items

1. **Organize pdb_structures/**
   - Move to `data/pdb_structures/`
   - Keep analysis tools in `examples/active_site_detection/`

2. **Create examples/ directory**
   - Move all `run_*.py`, `analyze_*.py` to appropriate subdirs
   - Add README for each example

3. **Create data/ directory**
   - Move all input files (`.txt`, `.fasta`)
   - Organize by analysis type

4. **P2RANK Integration**
   - Document proper P2RANK installation
   - Create wrapper for PrankWeb results
   - Show how to use P2RANK output with ASMC

5. **Clean .gitignore**
   - Ignore all `output*/`, `*_results/`, `*_asmc_*/`
   - Keep data files but ignore large outputs

## Next Steps

Would you like me to:
1. **Reorganize the directory structure** as recommended above?
2. **Create a proper P2RANK workflow** for UDH analysis?
3. **Document the integration** between custom tools and ASMC?
