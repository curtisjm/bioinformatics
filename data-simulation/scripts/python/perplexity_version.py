import os

import dask.array as da
import dask.dataframe as dd
import numpy as np
import pandas as pd
from dask.distributed import Client, progress

# Set up Dask client
client = Client()

BED_FILE = (
    "/Users/curtis/Documents/bioinformatics/data-simulation/real-data/10k_test.bed"
)
# BED_FILE = (
#     "/Users/curtis/Documents/bioinformatics/data-simulation/real-data/10k_test.bed"
# )
OUT_DIR_DATA = "./"
OUT_DIR_REGIONS = "./"
DEPTH = 25
NUM_SAMPLES = 1
STD_DEV = 0.15
READ_VARIATION = 0.15
ESTIMATED_NUM_DMRS = 1000
MIN_REGION_SIZE = 20
MAX_REGION_SIZE = 3000
# ESTIMATED_NUM_DMRS = 100
# MIN_REGION_SIZE = 5
# MAX_REGION_SIZE = 100
PERCENT_DIFF_TO_BE_CALLED_AS_DMR = 0.4
CHANCE_OF_INCREASE_IN_METHYLATION = 0.9

col_labels = ["chr", "start", "end", "uc", "mc", "prop"]
bed_data = dd.read_csv(BED_FILE, sep="\t", names=col_labels, header=None)


# Using a given proportion of methylation, simulate reads of each cytosine and
# mutating the unmethylateted counts, methylated counts, and proportion of methylation
# in the original data frame
def simulate_reads(prop: float, rng) -> tuple[int, int, float]:
    uc_count = 0
    mc_count = 0

    # Randomize the number of reads for each cytosine by adding a random number between
    # -READ_VARIATION * DEPTH and READ_VARIATION * DEPTH to the set depth
    num_reads = int(DEPTH + DEPTH * READ_VARIATION * (2 * rng.random() - 1))

    # Perform a weighted coin flip to determine if the cytosine is read as methylated or not by
    # generating a random number between 0 and 1 and checking if it is less than the true proportion of methylation
    # TODO: change this to normal distribution and find standard deviation
    random_values = 100 * rng.random(num_reads)
    mc_count = np.sum(random_values < prop, dtype=np.int32)
    uc_count = num_reads - mc_count

    sim_prop = 100 * mc_count / (mc_count + uc_count)

    return (uc_count, mc_count, sim_prop)


# For regions that are not DMRs, simulate the variation in reads
def simulate_reads_for_region(start: int, end: int, bed_data, rng) -> float:
    new_pm = 0

    for row in range(start, end):
        uc_count, mc_count, sim_prop = simulate_reads(bed_data.at[row, "prop"])
        bed_data.at[row, "uc"] = uc_count
        bed_data.at[row, "mc"] = mc_count
        bed_data.at[row, "prop"] = sim_prop

        new_pm += sim_prop

    return new_pm / (end - start)


def percent_diff(original_pm: float, start: int, end: int, bed_data) -> float:
    new_pm = bed_data.loc[start:end, "prop"].mean()
    return abs(new_pm - original_pm)


def produce_dmr_iter_rand(
    start: int, end: int, original_pm: float, inc_or_dec: str, bed_data, rng
) -> None:
    min_percent_diff = PERCENT_DIFF_TO_BE_CALLED_AS_DMR
    max_percent_diff = 1 - original_pm / 100
    goal_percent_diff = 100 * (
        rng.random() * (max_percent_diff - min_percent_diff) + min_percent_diff
    )
    # Select a random cytosine to change the methylation of and increase or decrease it's methlyation until the goal is reached
    while percent_diff(original_pm, start, end) < goal_percent_diff:
        # Select which cytosine we are modifying
        selected_cytosine = rng.integers(start, end, endpoint=True)
        inc_or_dec_multiplier = 1 if inc_or_dec == "+" else -1

        # Determine how much to change the methylation of the cytosine
        min_delta = 0
        max_delta = 100 - bed_data.at[selected_cytosine, "prop"]
        delta = rng.random() * (max_delta - min_delta) + min_delta
        bed_data.at[selected_cytosine, "prop"] += delta * inc_or_dec_multiplier

        # Use the new proportion of methylated reads to calculate the number of methylated and unmethylated reads
        total_num_reads = (
            bed_data.at[selected_cytosine, "uc"] + bed_data.at[selected_cytosine, "mc"]
        )
        bed_data.at[selected_cytosine, "mc"] = round(
            total_num_reads * bed_data.at[selected_cytosine, "prop"] / 100
        )
        bed_data.at[selected_cytosine, "uc"] = (
            total_num_reads - bed_data.at[selected_cytosine, "mc"]
        )

        # Update the proportion of methylated reads to reflected modified counts
        bed_data.at[selected_cytosine, "prop"] = (
            100 * bed_data.at[selected_cytosine, "mc"] / total_num_reads
        )

    return bed_data.loc[start:end, "prop"].mean()


def process_region(region):
    start, end, original_pm, is_dmr, inc_or_dec = region
    if is_dmr:
        new_pm = produce_dmr_iter_rand(start, end, original_pm, inc_or_dec)
    else:
        new_pm = simulate_reads_for_region(start, end)
    return pd.Series({"new_pm": new_pm, "percent_diff": abs(new_pm - original_pm)})


def define_regions_parallel():
    num_rows = bed_data.shape[0].compute()
    rng = np.random.default_rng()

    regions = []
    current_start = 0

    while current_start < num_rows - 1:
        region_size = rng.integers(MIN_REGION_SIZE, MAX_REGION_SIZE)
        current_end = min(current_start + region_size, num_rows - 1)

        region = bed_data.iloc[current_start:current_end].compute()
        original_pm = region["prop"].mean()

        is_dmr = int(rng.random() < CHANCE_OF_DMR)
        inc_or_dec = "+" if rng.random() < CHANCE_OF_INCREASE_IN_METHYLATION else "-"

        regions.append((current_start, current_end, original_pm, is_dmr, inc_or_dec))
        current_start = current_end + 1

    return dd.from_pandas(
        pd.DataFrame(
            regions,
            columns=["start_row", "end_row", "original_pm", "is_dmr", "inc_or_dec"],
        ),
        npartitions=client.ncores(),
    )


def process_chunk(chunk, bed_data, rng):
    results = []
    for region_num in chunk:
        modified_region, modified_bed_data = modification_handler(
            region_num, regions, bed_data, rng
        )
        results.append((modified_region, modified_bed_data))
    return results


if __name__ == "__main__":
    regions = define_regions_parallel()

    # Process regions in parallel
    results = regions.apply(
        process_region, axis=1, meta={"new_pm": float, "percent_diff": float}
    )

    # Compute results
    final_results = results.compute()

    # Combine results with original regions
    final_regions = pd.concat([regions.compute(), final_results], axis=1)

    # Output the simulated data
    output_data_filename = os.path.join(OUT_DIR_DATA, "false_pos_test.bed")
    bed_data.compute().to_csv(output_data_filename, sep="\t", index=False, header=False)

    output_region_filename = os.path.join(OUT_DIR_REGIONS, "false_pos_test_regions.tsv")
    final_regions.to_csv(output_region_filename, sep="\t", index=False, header=True)

    # Close the Dask client
    client.close()
