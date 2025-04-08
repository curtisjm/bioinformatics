import os
import numpy as np
import pandas as pd
import dask.dataframe as dd
from dask import delayed, compute
from typing import Tuple

BED_FILE = "../../real-data/10k_sorted.bed"
OUT_DIR = "./"
DEPTH = 25
READ_VARIATION = 0.15
BATCH_SIZE = 100


def simulate_reads_batch(props: list) -> list[Tuple[int, int, float]]:
    """Process a batch of proportion values at once"""
    np.random.seed()
    results = []

    for prop in props:
        uc_count = 0
        mc_count = 0

        num_reads = int(DEPTH + DEPTH * READ_VARIATION * (2 * np.random.rand() - 1))

        random_values = 100 * np.random.rand(num_reads)
        mc_count = np.sum(random_values < prop, dtype=np.int32)
        uc_count = num_reads - mc_count

        sim_prop = 100 * mc_count / (mc_count + uc_count)

        results.append((uc_count, mc_count, sim_prop))

    return results


if __name__ == "__main__":
    col_labels = ["chr", "start", "end", "uc", "mc", "prop"]
    bed_data = pd.read_csv(BED_FILE, sep="\t", names=col_labels, header=None)

    props = bed_data["prop"].tolist()
    batches = [props[i : i + BATCH_SIZE] for i in range(0, len(props), BATCH_SIZE)]

    delayed_results = [delayed(simulate_reads_batch)(batch) for batch in batches]

    batch_results = compute(*delayed_results)

    sim_data = [result for batch in batch_results for result in batch]

    for row, (uc, mc, pm) in enumerate(sim_data):
        bed_data.at[row, "uc"] = uc
        bed_data.at[row, "mc"] = mc
        bed_data.at[row, "prop"] = pm

    bed_data.to_csv("test.csv", sep="\t", index=False, header=False)
