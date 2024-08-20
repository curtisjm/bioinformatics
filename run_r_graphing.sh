#!/usr/bin/env bash

#SBATCH --job-name=run_r_graphing

#SBATCH --account=fc_williamslab

#SBATCH --partition=savio2_bigmem

#SBATCH --time=40:00:00

# TODO: replace path
R_SCRIPT="<path_to_R_script>"
WINDOW_SIZE=250000
FILE_SUFFIX=_all_XXX
Y_SCALE_MAX=100
# add ability to choose colors for the plots later

# parse command line options
while [ $# -gt 0 ]; do
    case "$1" in
        -w) WINDOW_SIZE="$2"; shift;;
        -s) FILE_SUFFIX="$2"; shift;;
        -y) Y_SCALE_MAX="$2"; shift;;
        *) break;;
    esac
    shift
done

for file in *gz; do
    gunzip $file
done