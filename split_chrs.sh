#!/usr/bin/env bash

#SBATCH --job-name=split_chromosomes

#SBATCH --account=fc_williamslab

#SBATCH --partition=savio2_bigmem

#SBATCH --time=40:00:00

mkdir split_chromosomes

# loop through all bed files in the directory and split them by separate chromosomes
for file in *bed; do
    for i in $(seq 1 5); do
        if [ i -eq 1 ]; then
            mkdir split_chromosomes/chr$i
        fi
        grep "Chr$i" $file > split_chromosomes/chr$i/$file
    done
done
