#!/bin/bash

#SBATCH --job-name=full_sample_gen_test

#SBATCH --account=fc_williamslab

#SBATCH --partition=savio3

#SBATCH --time=07:00:00

module load python
python ray_generate_samples.py
