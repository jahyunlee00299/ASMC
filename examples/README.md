# ASMC Examples

This directory contains example workflows and analysis scripts for using ASMC.

## Directory Structure

```
examples/
├── README.md                    # This file
├── udh_analysis/                # UDH (Uronate Dehydrogenase) analysis examples
├── alphafold_workflow/          # AlphaFold integration examples
├── active_site_detection/       # Active site detection tools
├── run_asmc_demo.py            # General ASMC demo
├── run_asmc_test.py            # ASMC testing script
├── setup_and_run.py            # Setup and run helper
└── 실행_예제.py                 # Korean language example
```

## Quick Start Examples

### 1. Basic ASMC Demo
```bash
cd examples
python run_asmc_demo.py
```

### 2. UDH Analysis
See `udh_analysis/README.md` for detailed UDH-specific workflows.

### 3. AlphaFold Integration
See `alphafold_workflow/README.md` for using AlphaFold models with ASMC.

### 4. Active Site Detection
See `active_site_detection/README.md` for structure-based active site identification.

## Data Requirements

All examples expect data files to be in `../data/`:
- PDB structures: `../data/pdb_structures/`
- FASTA sequences: `../data/*.fasta`
- Reference files: `../data/references*.txt`
- Pocket definitions: `../data/*_pocket.txt`

## Output

All examples write output to `../output/` by default.

## Contributing

When adding new examples:
1. Create a descriptive subdirectory
2. Add a README.md explaining the workflow
3. Include sample input/output
4. Document dependencies
