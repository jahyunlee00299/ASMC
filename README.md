<!-- ![ASMC logo](docs/asmc.png) -->
<!-- center and resize the image with html syntax -->
<p align="center">
  <img src="docs/asmc.png" alt="ASMC logo"
       width="140"
       height="178"
    />
</p>

# Active Site Modeling and Clustering (ASMC)

ASMC combines (i) homology modeling of family members (MODELLER), (ii) ligand-binding pocket search (P2RANK), (iii) structural alignment of modeled active sites (USalign) and (iv) density-based spatial clustering of obtained alignments (DBSCAN) in a single command line.

The clustering step can be carried out on either structural or sequence alignment and users can directly analyse their own set of protein 3D structures (e.g. AI-based models) by skipping the homology modeling step.

<!-- ![ASMC workflow](docs/ASMC_workflow.svg) -->
<!-- center and resize the image with html syntax -->
<p align="center">
  <img src="docs/ASMC_workflow.svg" alt="ASMC workflow" />
</p>

## Installation

### Installation with conda and pip

Download the latest GitHub release to obtain the code ([https://github.com/labgem/ASMC/releases](https://github.com/labgem/ASMC/releases)) and extract the code from the archive.

Then, use the following commands from the ASMC/ directory:
```
conda env create -n asmc -f env.yml
pip install ./
```

Conda will install all the python dependencies and two required third-party softwares:
- MODELLER (you still need to request the license key)
- USalign

The pip command is required to create the `asmc` command and use ASMC.

It's also possible to use only the `pip install ./` command, but this will not install any third party software.

### Third-party software dependencies

- P2RANK - for ligand-binding pocket detection ([https://github.com/rdk/p2rank](https://github.com/rdk/p2rank))
- MODELLER - for homology modeling ([https://salilab.org/modeller/](https://salilab.org/modeller/))
- USalign - for structural alignment ([https://github.com/pylelab/USalign](https://github.com/pylelab/USalign))

#### P2RANK setup

Download the p2rank tar.gz file (e.g: p2rank_2.5.tar.gz) and extract the archive.

Create a symbolic link related to the prank script, e.g:
```
ln -s <full_path_to>/p2rank_2.5/prank /usr/bin/prank
```

Modify the previous prank script to work with a symbolic link. At line 22, replace:
```bash
THIS_SCRIPT_DIR_REL_PATH=`dirname "${BASH_SOURCE[0]}"`
```
by
```bash
THIS_SCRIPT_DIR_REL_PATH=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
```

Now, ASMC can use P2RANK to detect ligand binding pockets.

### Installation with Docker

Follow the instructions in the [Docker section of the wiki](https://github.com/labgem/ASMC/wiki/Docker)

## Quick Usage

Run ASMC in a blind way (unknown active site) using a multi fasta file that should contain at least 100 sequences for clustering to be sufficiently relevant.

```
python -m asmc.run_asmc run --log run_asmc.log -t 6 -r reference_file -s sequences.fasta
```

Or if installed with pip, you can also use:
```
asmc run --log run_asmc.log -t 6 -r reference_file -s sequences.fasta
```

`reference_file` should contains the path to the reference(s) structure(s), e.g:
```
<path>/RefA.pdb
<path>/RefB.pdb
```

## Project Structure

```
ASMC/
├── asmc/                    # Core ASMC package
├── docs/                    # Documentation and diagrams
├── tests/                   # Unit tests
│
├── examples/                # Example workflows and analysis scripts
│   ├── udh_analysis/       # UDH enzyme family analysis
│   ├── alphafold_workflow/ # AlphaFold integration examples
│   ├── active_site_detection/  # Structure-based active site tools
│   └── README.md
│
├── data/                    # Input data files
│   ├── pdb_structures/     # PDB files and active site definitions
│   ├── *.fasta             # Sequence files
│   ├── *_references.txt    # Reference structure lists
│   └── *_pocket.txt        # Active site definitions
│
├── tools/                   # External tool integration
│   └── p2rank/             # P2RANK pocket detection integration
│
└── output/                  # Output directory (git-ignored)
```

## Documentation

- **English**: See the [wiki](https://github.com/labgem/ASMC/wiki) for detailed documentation
- **한국어 (Korean)**:
  - `Quick_Start.md` - Quick start guide in 5 minutes
  - `ASMC_사용법.md` - Comprehensive user manual
  - `examples/실행_예제.py` - Interactive execution script with menu

## Examples

### UDH Analysis
```bash
cd examples/udh_analysis
python run_udh_asmc.py --help
```

### AlphaFold Workflow
```bash
cd examples/alphafold_workflow
python run_asmc_with_alphafold.py --help
```

### Active Site Detection
```bash
cd examples/active_site_detection
python find_udh_active_sites.py --help
```

See `examples/README.md` for detailed usage instructions.

## Testing

Run the test suite:
```
python -m pytest tests/ -v
```

Or use example scripts:
```
cd examples
python run_asmc_demo.py
```
