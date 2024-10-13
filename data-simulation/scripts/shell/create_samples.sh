#!/bin/bash

#SBATCH --job-name=create_sample_data_10k_test

#SBATCH --account=fc_williamslab

#SBATCH --partition=savio3

#SBATCH --time=15:00:00

export MODULEPATH=${MODULEPATH}:/clusterfs/vector/home/groups/software/sl-7.x86_64/modfiles

BED_FILE="./real-bed-files/10k_test.bed"
OUT_DIR="./sample-bed-files"
DEPTH=25
NUM_SAMPLES=1
READ_VARIATION=0.15
AVAILABLE_CORES=32

# Optionally, you can specify command line options instead of changing the above values
while [ $# -gt 0 ]
do
    case "$1" in
        -f) BED_FILE="$2"; shift;;
        -o) OUT_DIR="$2"; shift;;
        -d) DEPTH="$2"; shift;;
        -n) NUM_SAMPLES="$2"; shift;;
        -v) READ_VARIATION="$2"; shift;;
        -c) AVAILABLE_CORES="$2"; shift;;
        *) break;;
    esac
    shift
done

# NOTE: need to cite gnu-parallel in any publication that uses it
    # O. Tange (2018): GNU Parallel 2018, Mar 2018, ISBN 9781387509881,
    # DOI https://doi.org/10.5281/zenodo.1146014
module load gnu-parallel
# @book{tange_ole_2018_1146014,
#       author       = {Tange, Ole},
#       title        = {GNU Parallel 2018},
#       publisher    = {Ole Tange},
#       month        = Mar,
#       year         = 2018,
#       ISBN         = {9781387509881},
#       doi          = {10.5281/zenodo.1146014},
#       url          = {https://doi.org/10.5281/zenodo.1146014}
# }

# Export necessary variables for GNU Parallel
export BED_FILE DEPTH READ_VARIATION AVAILABLE_CORES

simulate_reads() {
    local prop=$1
    local uc_count=0
    local mc_count=0

    # Generate a random number between 0 and 1 to determine the number of reads for each cytosine.
    # We will add a random number between -READ_VARIATION * DEPTH and READ_VARIATION * DEPTH to the depth
    # in order to simulate each cytosing being read a different number of times
    local num_reads=$(awk 'BEGIN { srand(); print int(DEPTH + DEPTH * READ_VARIATION * (2 * rand() - 1)) }')

    for ((j = 0; j < num_reads; j++))
    do
        # Perform a weighted coin flip to determine if the cytosine is read as methylated or not by
        # generating a random number between 0 and 1 and checking if it is less than the true proportion of methylation
        rand=$(awk 'BEGIN { srand(); print 100 * rand() }')
        if (( $(echo "$rand < $prop" | bc -l) ))
        then
            mc_count=$((mc_count + 1))
        else
            uc_count=$((uc_count + 1))
        fi
    done

    local sim_prop=$(awk -v mc=$mc_count -v uc=$uc_count 'BEGIN { print 100 * mc / (mc + uc) }')

    # Return the results
    echo "$uc_count\t$mc_count\t$sim_prop"
}

export -f simulate_reads

# Create a given number of new sample bed files
for i in $(seq 1 $NUM_SAMPLES)
do
    # Extract list of true proportions of methylation and process 20 at a time in parallel
    num_blocks=$(echo "$(wc -l < "$BED_FILE") / $AVAILABLE_CORES" | bc -l)
    results=$(awk '{print $6}' "$BED_FILE" | parallel -j $AVAILABLE_CORES --pipe-part --block $num_blocks -k simulate_reads {})

    simulated_data=""

    while read -r line
    do
        simulated_data="${simulated_data}${line}\n"
    done <<< "$results"
    
    out_file=$(basename "${BED_FILE}" .bed)_sample_${i}.bed
    # Merge the original chromosome, start, and end columns with the simulated unmethylated and methylated cytosines and the simulated proportions
    paste -d "\t" <(awk 'BEGIN {OFS="\t"} {print $1, $2, $3}' "$BED_FILE") <(echo -e "$simulated_data") > "${OUT_DIR}/${out_file}"

    echo "Sample $i completed: ${out_file}"
done


