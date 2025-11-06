# CPLEX Setup for CarveWe

CarveWe requires IBM CPLEX with an academic license for optimal performance.

## Quick Setup (3 Steps)

### 1. Install CarveWe
```bash
conda create -n carvewe python=3.12
conda activate carvewe
conda install -c ryanreyn carvewe
```

### 2. Install CPLEX Binaries and CPLEX Python API
Download CPLEX version 22.1.2 from the IBM Academic Initiative:
https://www.ibm.com/academic/technology/data-science

Then install to your conda environment:
```bash
./cplex_studio2212.linux_x86_64.bin -i console \
  -DUSER_INSTALL_DIR=$CONDA_PREFIX/ilog \
  -DLICENSE_ACCEPTED=TRUE
```

Install the publicly available community version of CPLEX into your new environment to access the Python API:
```bash
pip install cplex
```

### 3. Activate Full CPLEX
```bash
# Reactivate your environment to apply CPLEX configuration
conda deactivate
conda activate carvewe

# Verify installation
carvewe-check-cplex
```

## Verification
Run the diagnostic tool to check your setup:
```bash
carvewe-check-cplex
```

This will show:
- ✓ Whether CPLEX binaries are installed
- ✓ Whether you have full or community edition
- ✓ Whether COBRApy can use CPLEX as a solver

## Troubleshooting

**"COMMUNITY EDITION detected"**
- You need to install the full CPLEX binaries from IBM
- Follow step 2 above

**"cplex module not found"**
- This should be installed automatically with CarveWe
- Manually install: `pip install cplex`

**"CPLEX not available as solver"**
- Reactivate your environment: `conda deactivate && conda activate carvewe`
- Run diagnostic: `carvewe-check-cplex`

## How It Works

CarveWe uses a "hybrid" approach:
1. `pip install cplex` provides the Python API (automatically installed with CarveWe)
2. CPLEX binaries in `$CONDA_PREFIX/ilog` provide the full solver engine
3. Environment activation scripts link them together

This gives you:
- ✓ Easy Python API installation (no setup.py hunting!)
- ✓ Full academic license capabilities
- ✓ Clean conda environment isolation