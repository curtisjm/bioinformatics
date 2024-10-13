import os

import dask.dataframe as dd
import numpy as np
import pandas as pd
from dask.delayed import delayed
from dask.distributed import Client, LocalCluster

BED_FILE = "../../real-data/10k_test.bed"
OUT_DIR = "./"
DEPTH = 25
READ_VARIATION = 0.15


@delayed
def simulate_reads(prop: float) -> tuple[int, int, float]:
    np.random.seed()
    num_reads = int(DEPTH + DEPTH * READ_VARIATION * (2 * np.random.rand() - 1))
    random_values = 100 * np.random.rand(num_reads)
    mc_count = np.sum(random_values < prop, dtype=np.int32)
    uc_count = num_reads - mc_count
    sim_prop = 100 * mc_count / (mc_count + uc_count)
    return (uc_count, mc_count, sim_prop)


if __name__ == "__main__":
    # Set up Dask client
    cluster = LocalCluster()
    client = Client(cluster)

    # Read the bed file
    col_labels = ["chr", "start", "end", "uc", "mc", "prop"]
    bed_data = dd.read_csv(BED_FILE, sep="\t", names=col_labels, header=None)

    # Apply simulate_reads to each row in parallel
    sim_data = bed_data["prop"].apply(simulate_reads, meta=("sim_data", "object"))

    # Compute the results
    results = sim_data.compute()

    # Update the DataFrame with the results
    bed_data["uc"] = [r[0] for r in results]
    bed_data["mc"] = [r[1] for r in results]
    bed_data["prop"] = [r[2] for r in results]

    # Write the results to a file
    # out_file = os.path.join(
    #     OUT_DIR, f"{os.path.basename(BED_FILE).replace('.bed', '_simulated.bed')}"
    # )
    # bed_data.to_csv(out_file, sep="\t", index=False, header=False)
    bed_data.to_csv("asdf.bed", sep="\t", index=False, header=False)

    # Close the client
    client.close()
