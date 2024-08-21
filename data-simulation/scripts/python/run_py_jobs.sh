#!/bin/bash

#SBATCH --job-name=false_positive_test

#SBATCH --account=fc_williamslab

#SBATCH --partition=savio3

#SBATCH --time=00:30:00

module load python
python false_positive.py
