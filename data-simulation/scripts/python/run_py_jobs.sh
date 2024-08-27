#!/bin/bash

#SBATCH --job-name=fp_test_parallel

#SBATCH --account=fc_williamslab

#SBATCH --partition=savio3

#SBATCH --time=02:00:00

export MODULEPATH=${MODULEPATH}:/clusterfs/vector/home/groups/software/sl-7.x86_64/modfiles

module load python
python false_positive.py
