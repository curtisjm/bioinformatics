import ray
import os
import pandas as pd
import numpy as np

BED_FILE = "./real-bed-files/D23_Col0_all_CpG.bed"
OUT_DIR = "./sample-bed-files"
DEPTH = 25
NUM_SAMPLES = 1
READ_VARIATION = 0.15

# Initialize ray to manage parallel tasks
ray.init()

# Load data from bed file into a pandas dataframe
col_labels = ["chr", "start", "end", "uc", "mc", "prop"]
bed_data = pd.read_csv(BED_FILE, sep='\t', names=col_labels, header=None)

# Using a given proportion of methylation, simulate reads of each cytosine and
# mutating the unmethylateted counts, methylated counts, and proportion of methylation
# in the original data frame
@ray.remote
def simulate_reads(prop: float, row: int) -> tuple[int, int, float]:
    np.random.seed()
    uc_count = 0
    mc_count = 0

    # Randomize the number of reads for each cytosine by adding a random number between
    # -READ_VARIATION * DEPTH and READ_VARIATION * DEPTH to the set depth
    num_reads = int(DEPTH + DEPTH * READ_VARIATION * (2 * np.random.rand() - 1))

    # Perform a weighted coin flip to determine if the cytosine is read as methylated or not by
    # generating a random number between 0 and 1 and checking if it is less than the true proportion of methylation
    random_values = 100 * np.random.rand(num_reads)
    mc_count = np.sum(random_values < prop)
    uc_count = num_reads - mc_count

    sim_prop = 100 * mc_count / (mc_count + uc_count)
    
    return (uc_count, mc_count, sim_prop)

# Generate as many samples as specified in a single job
for i in range(NUM_SAMPLES):
    # Use ray to run the simulate_reads function in parallel across all available cores
    futures = [simulate_reads.remote(prop, i) for i, prop in enumerate(bed_data["prop"])]
    # Request the results of the parallel computations
    results = ray.get(futures)

    # Overwrite the real data stored in the dataframe with the sample data
    for row, (uc_count, mc_count, sim_prop) in enumerate(results):
        bed_data.at[row, "uc"] = uc_count
        bed_data.at[row, "mc"] = mc_count
        bed_data.at[row, "prop"] = sim_prop

    # Output the simulated data to a new bed file
    out_file = os.path.join(OUT_DIR, f"{os.path.basename(BED_FILE).replace('.bed', '')}_sample_{i}_ray.bed")
    bed_data.to_csv(os.path.join(out_file), sep='\t', index=False, header=False)
