# Active Site Detection Tools

Structure-based tools for identifying and analyzing enzyme active sites.

## Scripts

- **find_udh_active_sites.py** - Alignment-based active site detection
- **predict_pockets_simple.py** - Simple geometric pocket prediction  
- **predict_pockets_enhanced.py** - Multi-feature pocket scoring
- **analyze_active_site_burial.py** - Analyze burial depth of known sites

## Recommended Method

**Structure Alignment** (find_udh_active_sites.py) - Best accuracy (~80% F1-score)

```bash
python find_udh_active_sites.py \
    --reference ../../data/pdb_structures/AtUdh_3rfv.pdb \
    --reference-sites ../../data/pdb_structures/atudh_active_site.txt \
    --target protein.pdb \
    --output ../../output/
```

For P2RANK integration, see `../../tools/p2rank/README.md`
