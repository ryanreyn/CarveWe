# CarveWe
## Overview
**CarveWe** is a comprehensive workflow that integrates other open source tools, chiefly the CarveMe model annotation [software](https://github.com/cdanielmachado/carveme), and offers several distinct applications for metabolic modeling. Principally the complete workflow was designed to apply SOM neural networks to feature data of metabolic sensitivities based on *in silico* nutrient limitation tests in order to group marine heterotrophic organisms into functional cohorts. These cohorts were then contextualized with information from a variety of data streams. This effort is described in full detail in this [publication](https://www.biorxiv.org/content/10.1101/2024.05.29.596556v1.abstract). We offer **CarveWe** here as a CLI based tool that allows users to align their own genomes to the genomes used to develop our metabolic clusters, as well as annotate metabolic sensitivities *de novo* with their own selection of genomes. Currently, **CarveWe** supports both amino acid (.faa) and nucleotide (.fna) fasta files, though CarveMe performs optimally with amino acid fasta files. For more information on the workflow, including how to install and use it feel free to refer to the [wiki](https://github.com/ryanreyn/CarveWe/wiki).

## Installation

### Prerequisites
CarveWe requires IBM CPLEX with an academic license. [See CPLEX setup guide](docs/CPLEX_SETUP.md).

### Quick Install
First create a new conda environment and install the carvewe package into this environment:
```bash
conda create -n carvewe python=3.12
conda activate carvewe
conda install -c ryanreyn carvewe
```

Then install the community edition CPLEX that is publicly available to access the Python API
```bash
pip install cplex
```

Then install your CPLEX binaries while the environment is active (to ensure the binaries correctly identify `$CONDA_PREFIX`):
```bash
./cplex_studio_installer.bin -i console \
  -DUSER_INSTALL_DIR=$CONDA_PREFIX/ilog \
  -DLICENSE_ACCEPTED=TRUE
```

Deactivate and reactivate your environment. Your environment will automatically be configured to enable your fully licensed CPLEX which can be verified with `carvewe-check-cplex`:
```bash
conda deactivate
conda activate carvewe

carvewe-check-cplex
```

## Genome alignment
For users interested in a streamlined method to align their genomes of interest to the genomes associated with each of our 8 metabolic clusters, we provide a CLI integrated tool *genome-aligner*. With *genome-aligner* a user can simply provide the path to an individual fasta file or directory of fasta files and generate a report detailing the closest matching genome and its associated metabolic cluster. Currently, this tool supports nucleotide and amino acid fasta sequences. This alignment is done by performing a best hit BLAST from the query genome(s) to a pre-made blast database of the fasta sequences from the CarveWe paper using the appropriate BLAST program for the input file type.

## Metabolic sensitivity prediction
One novel key element of the CarveWe pipeline that we introduced was the notion of performing media predictions on the annotated models with *in silico* conditions meant to mimic a nutrient replete state (unbounded import flux) as well as nutrient deplete states. The **CarveWe** tools also offers a more extensive workflow that enables the user to annotate new CarveMe model ensembles, predict metabolic fluxes and model growth under nutrient replete conditions, and then predict metabolic sensitivities by calculating model growth under a user-prescribed level of nutrient limitation. This remains under development and will be released in full soon with active development to follow.

## Publication reproduction
Reproducible science is essential, so this repository also houses the relevant scripts and publicly available data to reproduce the complete analytical findings from our analysis marine microbial metabolism at scale using GEMs constructed with CarveMe, corresponding to this [publication](https://www.biorxiv.org/content/10.1101/2024.05.29.596556v1.abstract). The *reproduce_publication* directory in this repository contains the necessary codes, principally an R markdown file that is designed to ingest the data housed in the 
*Publication_Data* subdirectory and generate a markdown document that regenerates the results that we report in our paper, including main and supplemental figures.