import os

import dask
import dask.dataframe as dd
import numpy as np
import pandas as pd
from dask.distributed import LocalCluster

BED_FILE = "../../real-data/10k_sorted.bed"
OUT_DIR = "./"
DEPTH = 25
READ_VARIATION = 0.15


def simulate_reads(prop: float) -> tuple[int, int, float]:
    np.random.seed()
    uc_count = 0
    mc_count = 0

    num_reads = int(DEPTH + DEPTH * READ_VARIATION * (2 * np.random.rand() - 1))

    random_values = 100 * np.random.rand(num_reads)
    mc_count = np.sum(random_values < prop, dtype=np.int32)
    uc_count = num_reads - mc_count

    sim_prop = 100 * mc_count / (mc_count + uc_count)

    return (uc_count, mc_count, sim_prop)


if __name__ == "__main__":
    # client = cluster.get_client()
    # cluster = LocalCluster()
    col_labels = ["chr", "start", "end", "uc", "mc", "prop"]
    bed_data = pd.read_csv(BED_FILE, sep="\t", names=col_labels, header=None)
    # dask_bed_data = dd.from_pandas(bed_data, npartitions=bed_data.shape[1])
    # dask_bed_data = dask_bed_data.persist()
    # sim_data = [client.submit(simulate_reads, prop) for prop in bed_data["prop"]]
    # results = sim_data.result()

    sim_data = [simulate_reads(prop) for prop in bed_data["prop"]]

    for row, (uc, mc, pm) in enumerate(sim_data):
        bed_data.at[row, "uc"] = uc
        bed_data.at[row, "mc"] = mc
        bed_data.at[row, "prop"] = pm

    # print(bed_data.iloc[0])
    # print(simulate_reads(bed_data.iloc[0]))

    # bed_data = dd.read_csv(BED_FILE, sep="\t", names=col_labels, header=None)
    # bed_data["uc"], bed_data["mc"], bed_data["prop"] = dask.compute(
    #     bed_data["prop"].map_partitions(simulate_reads)
    # )

    bed_data.to_csv("test.csv", sep="\t", index=False, header=False)
    # out_file = os.path.join(OUT_DIR, f"{os.path.basename(BED_FILE).replace('.bed', '')})
    # out_file = os.path.join(OUT_DIR, "test_fn.bed")
    # bed_data.to_csv(os.path.join(out_file), sep="\t", index=False, header=False)
