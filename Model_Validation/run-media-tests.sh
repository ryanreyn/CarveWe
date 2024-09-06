#!/bin/bash
#This is a small bash script to call the python script and run it when the slurm file is submitted
eval "$(conda shell.bash hook)"
conda activate /project/nmherrer_110/tools/.conda/envs/CarveMe

module load python
python rescue-sole-C-media.py
