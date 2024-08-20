import os
import ipyparallel as ipp
import pandas as pd
import numpy as np
from pandas import DataFrame

BED_FILE = "./curtis_testing.bed"
OUT_DIR = "./"
DEPTH = 25
NUM_SAMPLES = 1
READ_VARIATION = 0.15

# Load data from bed file into a pandas dataframe
col_labels = ["chr", "start", "end", "uc", "mc", "prop"]
bed_data = pd.read_csv(BED_FILE, sep='\t', names=col_labels, header=None)

# Use ipyparallel to parallelize independent computations across multiple cores on the node
# mycluster = ipp.Cluster(n = int(os.getenv('SLURM_CPUS_ON_NODE')))
mycluster = ipp.Cluster(n=4)
c = mycluster.start_and_connect_sync()

c[:].use_cloudpickle()

# Use a load-balanced view (sequentially dispatching computational tasks as earlier computational tasks finish
# to run the simulate_reads function in parallel across all available cores 
lview = c.load_balanced_view()

# Cause execution on main process to wait while tasks sent to workers finish
lview.block = True

# Import numpy for use in parallel jobs
# with lview.sync_imports():
#     import numpy as np

# Using a given proportion of methylation, simulate reads of each cytosine and
# mutating the unmethylateted counts, methylated counts, and proportion of methylation
# in the original data frame
# @lview.parallel(block=True)
def simulate_reads(prop: float) -> None:
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
    
    # bed_data.at[row, "uc"] = uc_count
    # bed_data.at[row, "mc"] = mc_count
    # bed_data.at[row, "prop"] = sim_prop
    return (uc_count, mc_count, sim_prop)

# @lview.parallel(block=True)
def wrapper(row: int):
    return (simulate_reads(original_props[row], row, bed_data))


for i in range(NUM_SAMPLES):
    # Pass the simulate_reads function the proportion of methylation and the corresponding row on each call
    # params = zip(bed_data["prop"].copy(), range(bed_data.shape[0]))
    
    # Run the simulate_reads function on each row in the bed file in with ipyparallel
    sim_cytosines = lview.map(simulate_reads, bed_data["prop"])
    for j in range(len(sim_cytosines)):
        uc, mc, prop = sim_cytosines[i]
        bed_data.at[j, "uc"] = uc
        bed_data.at[j, "mc"] = mc
        bed_data.at[j, "prop"] = prop

    # Output the simulated data to a new bed file
    out_file = os.path.join(OUT_DIR, f"{os.path.basename(BED_FILE).replace('.bed', '')}_sample_{i}_py.bed")
    bed_data.to_csv(out_file, sep='\t', index=False, header=False)
