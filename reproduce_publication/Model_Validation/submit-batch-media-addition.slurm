#!/bin/bash
#SBATCH -p qcb
#SBATCH -A nmherrer_110
#SBATCH -t 7-00:00:00
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=8
#SBATCH --out=sub_out
#SBATCH --error=sub_err
module purge

eval "$(conda shell.bash hook)"
conda activate /project/nmherrer_110/tools/.conda/envs/CarveMe

python sub_python
