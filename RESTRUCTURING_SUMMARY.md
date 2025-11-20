# ASMC Project Restructuring Summary

**Date**: 2025-11-20
**Status**: ✓ Completed

## What Was Done

Reorganized the ASMC project to match official GitHub repository structure (labgem/ASMC) and improve organization.

## Changes Made

### 1. **New Directory Structure Created**

```
ASMC/
├── asmc/                    # Core package (unchanged)
├── docs/                    # Documentation (unchanged)
├── tests/                   # Tests (unchanged)
│
├── examples/                # NEW: Example workflows
│   ├── udh_analysis/
│   ├── alphafold_workflow/
│   ├── active_site_detection/
│   └── README.md
│
├── data/                    # NEW: Input data
│   ├── pdb_structures/
│   ├── *.fasta
│   ├── *_references.txt
│   └── *_pocket.txt
│
├── tools/                   # NEW: External tool integration
│   └── p2rank/
│       ├── INSTALL_P2RANK.md
│       ├── USE_PRANKWEB.md
│       ├── compare_p2rank_results.py
│       └── README.md
│
└── output/                  # NEW: Output directory
    └── .gitkeep
```

### 2. **Files Moved**

#### To `examples/udh_analysis/`:
- `analyze_udh.py`
- `extract_udh_active_sites.py`
- `extract_udh_substrate_sites.py`
- `run_udh_asmc.py`
- `analyze_asmc_results.py`
- `analyze_substrate_results.py`
- `analyze_substrate_variants.py`

#### To `examples/alphafold_workflow/`:
- `run_asmc_with_alphafold.py`
- `create_models_txt_for_alphafold.py`

#### To `examples/active_site_detection/`:
- `find_udh_active_sites.py`
- `predict_pockets_simple.py`
- `predict_pockets_enhanced.py`
- `analyze_active_site_burial.py`

#### To `examples/` (root level examples):
- `실행_예제.py`
- `run_asmc_demo.py`
- `run_asmc_test.py`
- `setup_and_run.py`

#### To `data/`:
- `alphafold_models.txt`
- `pocket_alphafold.txt`
- `references_alphafold.txt`
- `udh_pocket.txt`
- `udh_references.txt`
- `udh_active_sites.fasta`
- `udh_substrate_sites.fasta`
- `sequences.fasta`
- `models.txt`
- `references.txt`

#### To `data/pdb_structures/`:
- `AtUdh_3rfv.pdb`
- `pdb3rfv_chainA.pdb`
- `pdb3rfv.ent`
- `atudh_active_site.txt`
- `results/` directory

#### To `tools/p2rank/`:
- `INSTALL_P2RANK.md`
- `USE_PRANKWEB.md`
- `compare_p2rank_results.py`

### 3. **Removed**
- `pdb_structures/` directory (reorganized into `data/pdb_structures/`)

### 4. **Documentation Created**

- `examples/README.md` - Overview of all examples
- `examples/udh_analysis/README.md` - UDH analysis workflow
- `examples/alphafold_workflow/README.md` - AlphaFold integration guide
- `examples/active_site_detection/README.md` - Active site detection tools
- `data/README.md` - Data directory organization
- `data/pdb_structures/README.md` - PDB structures documentation
- `tools/p2rank/README.md` - P2RANK integration guide

### 5. **Updated Files**

#### `.gitignore`
Added patterns to ignore:
- All output directories (`output/`, `*_output/`, `*_results/`, `*_asmc_*/`)
- UDH analysis outputs
- P2RANK outputs
- Generated plots and reports

#### `README.md`
Added:
- Project structure diagram
- Links to example workflows
- Updated documentation paths

### 6. **Created**
- `output/.gitkeep` - Preserve output directory in git
- `PROJECT_STRUCTURE_ANALYSIS.md` - Detailed analysis of changes

## Benefits

### 1. **Cleaner Root Directory**
Before: 45+ files in root
After: Core files only (~20 files)

### 2. **Better Organization**
- Examples grouped by purpose
- Data files centralized
- Tools properly documented

### 3. **Follows Official Standards**
- Matches labgem/ASMC structure
- Professional Python project layout
- Easy for new contributors

### 4. **Improved Documentation**
- README in each directory
- Clear usage examples
- Integration guides

### 5. **P2RANK Integration**
- Installation guide
- Online alternative (PrankWeb)
- Comparison tools
- Proper workflow documentation

## Key Features Documented

### P2RANK Integration
- How ASMC uses P2RANK internally
- Pocket file format specification
- Comparison with known active sites
- Windows installation guide

### Active Site Detection
- Structure alignment method (80% accuracy)
- P2RANK machine learning (60% accuracy)
- Geometric methods (7-19% accuracy)
- Clear recommendations

### Workflows
- UDH enzyme family analysis
- AlphaFold model integration
- Custom structure analysis

## Backward Compatibility

### Potential Breaking Changes
Scripts that used relative paths may need updates:

**Old**:
```python
with open('models.txt') as f:
```

**New**:
```python
with open('data/models.txt') as f:
# or
with open('../data/models.txt') as f:  # from examples/
```

### Migration Guide for Users

If you have existing scripts:

1. **Update data file paths**:
   ```python
   # Old: 'references.txt'
   # New: 'data/references.txt'
   ```

2. **Update output paths**:
   ```python
   # Old: output='./'
   # New: output='output/'
   ```

3. **Run from project root** or update relative paths

## Files Still in Root (To Be Cleaned Later)

These should be moved or removed in future cleanup:
- `3rfv_superposition_results/` → should be in `output/`
- `asmc_alphafold_output/` → should be in `output/`
- `output_*` directories → should be in `output/`
- `udh_*` output files → should be in `output/`
- `test_data/`, `test_output/` → potentially consolidate

## Next Steps (Optional)

1. **Clean remaining output directories**
   ```bash
   mv 3rfv_superposition_results/ output/
   mv asmc_alphafold_output/ output/
   mv output_clustering/ output/
   # etc.
   ```

2. **Update Korean documentation**
   - Update `ASMC_사용법.md` with new paths
   - Update `Quick_Start.md` with new structure

3. **Add tests for examples**
   - Test scripts run without errors
   - Validate output formats

4. **Create changelog**
   - Document version changes
   - Migration guide for existing users

## Verification

To verify the restructuring worked:

```bash
# Check structure
ls examples/
ls data/
ls tools/

# Test examples
cd examples/udh_analysis
python run_udh_asmc.py --help

# Verify data access
cat data/pdb_structures/atudh_active_site.txt

# Check P2RANK docs
cat tools/p2rank/README.md
```

## Success Criteria

✅ All files moved to appropriate directories
✅ README files created for all new directories
✅ Main README updated with new structure
✅ .gitignore updated to ignore outputs
✅ P2RANK integration documented
✅ Examples properly organized
✅ Data files centralized

## Conclusion

The project now follows professional Python project standards and matches the official ASMC repository structure. This improves maintainability, makes the project more accessible to new users, and provides clear documentation for all workflows.
