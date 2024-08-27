#!/bin/bash

#SBATCH --job-name=false_positive_test

#SBATCH --account=fc_williamslab

#SBATCH --partition=savio3

#SBATCH --time=20:00:00

module load python
python non_parallel_fp.py
