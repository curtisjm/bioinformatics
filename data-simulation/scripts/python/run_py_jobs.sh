#!/bin/bash

#SBATCH --job-name=fp_test_parallel

#SBATCH --account=fc_williamslab

#SBATCH --partition=savio3_bigmem

#SBATCH --time=02:00:00

module load python
python false_positive.py
