#!/bin/bash
#SBATCH -p qcb
#SBATCH -A nmherrer_110
#SBATCH -t 7-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --output=media_out
#SBATCH --error=media_err

./run-media-tests.sh
