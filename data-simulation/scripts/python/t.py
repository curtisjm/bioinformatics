import os
import numpy as np
import pandas as pd
import dask.dataframe as dd
from dask import delayed
from dask.distributed import Client, LocalCluster
from numpy.random import Generator, SeedSequence
from pandas import DataFrame, Series

BED_FILE = "../../real-data/10k_sorted.bed"
OUT_DIR_DATA = "./"
OUT_DIR_REGIONS = "./"
DEPTH = 25
NUM_SAMPLES = 1
STD_DEV = 0.15
READ_VARIATION = 0.15
ESTIMATED_NUM_DMRS = 100
MIN_REGION_SIZE = 5
MAX_REGION_SIZE = 100
PERCENT_DIFF_TO_BE_CALLED_AS_DMR = 0.4
CHANCE_OF_INCREASE_IN_METHYLATION = 0.9


def get_regions(bed_data) -> DataFrame:
    rng = np.random.default_rng()
    num_dmrs = 0
    cols = [
        "start_row",
        "end_row",
        "start_coord",
        "end_coord",
        "original_pm",
        "new_pm",
        "percent_diff",
        "is_dmr",
        "inc_or_dec",
    ]
    regions_df = DataFrame(columns=cols)
    current_start = 0
    num_rows = bed_data.shape[0]

    estimated_num_regions = num_rows / ((MAX_REGION_SIZE - MIN_REGION_SIZE) / 2)
    chance_of_dmr = ESTIMATED_NUM_DMRS / estimated_num_regions

    while current_start < num_rows - 1:
        region_size = rng.integers(MIN_REGION_SIZE, MAX_REGION_SIZE)
        current_end = current_start + region_size

        if current_end >= bed_data.shape[0]:
            current_end = bed_data.shape[0] - 1

        while bed_data.at[current_start, "chr"] != bed_data.at[current_end, "chr"]:
            current_end -= 1

        region = bed_data.iloc[current_start:current_end]
        original_pm = region["prop"].mean()

        is_dmr = int(rng.random() < chance_of_dmr)
        num_dmrs += is_dmr

        if is_dmr:
            inc_or_dec = (
                "+" if rng.random() < CHANCE_OF_INCREASE_IN_METHYLATION else "-"
            )
            if (
                inc_or_dec == "+"
                and (100 - original_pm) <= 100 * PERCENT_DIFF_TO_BE_CALLED_AS_DMR
            ):
                inc_or_dec = "-"
            if (
                inc_or_dec == "-"
                and original_pm <= 100 * PERCENT_DIFF_TO_BE_CALLED_AS_DMR
            ):
                inc_or_dec = "+"
        else:
            inc_or_dec = "N/A"

        new_row = DataFrame(
            [
                [
                    current_start,
                    current_end,
                    bed_data.at[current_start, "start"],
                    bed_data.at[current_end, "end"],
                    original_pm,
                    0.0,
                    0.0,
                    is_dmr,
                    inc_or_dec,
                ]
            ],
            columns=cols,
        )
        regions_df = (
            pd.concat([regions_df, new_row], ignore_index=True)
            if not regions_df.empty
            else new_row
        )
        current_start = current_end + 1
    return regions_df


def percent_diff(original_pm: float, start: int, end: int, region: DataFrame) -> float:
    new_pm = region["prop"].mean()
    return abs(new_pm - original_pm)


def simulate_dmr(
    region_data: DataFrame,
    region_info: Series,
    seed: int,
) -> DataFrame:
    # Create a new RNG with the provided seed
    rng = np.random.default_rng(seed)

    # Create a copy of the region data to avoid modifying the original
    region = region_data.copy()

    start = region_info.at["start_row"]
    end = region_info.at["end_row"]
    original_pm = region_info.at["original_pm"]
    inc_or_dec = region_info.at["inc_or_dec"]

    # If not a DMR, return the region unchanged
    if not region_info.at["is_dmr"]:
        return region

    max_percent_diff = 1 - original_pm / 100
    min_percent_diff = PERCENT_DIFF_TO_BE_CALLED_AS_DMR
    goal_percent_diff = 100 * (
        rng.random() * (max_percent_diff - min_percent_diff) + min_percent_diff
    )

    iteration = 0
    max_iterations = MAX_REGION_SIZE * 10

    while (
        percent_diff(original_pm, start, end, region) < goal_percent_diff
        and iteration < max_iterations
    ):
        # Get index relative to the region DataFrame
        region_size = len(region)
        selected_index = rng.integers(0, region_size, endpoint=True)
        inc_or_dec_multiplier = 1 if inc_or_dec == "+" else -1

        min_delta = 0
        max_delta = (
            100 - region.iloc[selected_index]["prop"]
            if inc_or_dec_multiplier == 1
            else region.iloc[selected_index]["prop"]
        )
        delta = rng.random() * (max_delta - min_delta) + min_delta

        # Update the prop value
        region.iloc[selected_index, region.columns.get_loc("prop")] += (
            delta * inc_or_dec_multiplier
        )

        # Update mc and uc based on the new prop
        total_num_reads = (
            region.iloc[selected_index]["uc"] + region.iloc[selected_index]["mc"]
        )
        new_mc = round(total_num_reads * region.iloc[selected_index]["prop"] / 100)
        new_uc = total_num_reads - new_mc

        region.iloc[selected_index, region.columns.get_loc("mc")] = new_mc
        region.iloc[selected_index, region.columns.get_loc("uc")] = new_uc

        # Recalculate prop to ensure consistency
        region.iloc[selected_index, region.columns.get_loc("prop")] = (
            100 * new_mc / total_num_reads if total_num_reads > 0 else 0
        )

        iteration += 1

    return region


def process_region(region_info, bed_data, seed):
    start_row = region_info["start_row"]
    end_row = region_info["end_row"]
    region_data = bed_data.iloc[start_row : end_row + 1].copy()

    # Process the region
    processed_region = simulate_dmr(region_data, region_info, seed)

    # Calculate new metrics for the region
    new_pm = processed_region["prop"].mean()
    percent_diff_val = abs(new_pm - region_info["original_pm"])

    # Update region info
    region_info = region_info.copy()
    region_info["new_pm"] = new_pm
    region_info["percent_diff"] = percent_diff_val

    return processed_region, region_info


if __name__ == "__main__":
    # Set up Dask client
    cluster = LocalCluster(processes=True)
    client = Client(cluster)
    print(f"Dask dashboard available at: {client.dashboard_link}")

    # Read the bed data
    col_labels = ["chr", "start", "end", "uc", "mc", "prop"]
    bed_data = pd.read_csv(BED_FILE, sep="\t", names=col_labels, header=None)

    # Get regions (this is fast enough to run sequentially)
    regions_info = get_regions(bed_data)

    # Create a seed sequence for random number generation
    seed_seq = SeedSequence()
    # Generate seeds for each region
    seeds = seed_seq.spawn(len(regions_info))

    # Process regions in parallel
    delayed_results = []

    for i, (_, region) in enumerate(regions_info.iterrows()):
        # Create a delayed task for each region
        delayed_result = delayed(process_region)(region, bed_data, seeds[i])
        delayed_results.append(delayed_result)

    # Compute all delayed tasks
    results = client.compute(delayed_results)
    results = client.gather(results)

    # Separate processed regions and updated region info
    processed_regions = [result[0] for result in results]
    updated_regions_info = pd.DataFrame([result[1] for result in results])

    # Combine processed regions
    processed_bed_data = pd.concat(processed_regions)

    # Write output files
    output_data_filename = os.path.join(OUT_DIR_DATA, "new_fp.bed")
    processed_bed_data.to_csv(output_data_filename, sep="\t", index=False, header=False)

    output_region_filename = os.path.join(OUT_DIR_REGIONS, "new_fp_regions.tsv")
    updated_regions_info.to_csv(
        output_region_filename, sep="\t", index=False, header=True
    )

    # Close the client
    client.close()
    cluster.close()
