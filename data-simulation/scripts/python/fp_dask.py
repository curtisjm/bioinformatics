import os

import numpy as np
import pandas as pd
from numpy.random import Generator
from pandas import DataFrame, Series

BED_FILE = "../../real-data/10k_sorted.bed"
OUT_DIR_DATA = "./"
OUT_DIR_REGIONS = "./"
DEPTH = 25
NUM_SAMPLES = 1
STD_DEV = 0.15
READ_VARIATION = 0.15
# ESTIMATED_NUM_DMRS = 1000
# MIN_REGION_SIZE = 20
# MAX_REGION_SIZE = 3000
ESTIMATED_NUM_DMRS = 100
MIN_REGION_SIZE = 5
MAX_REGION_SIZE = 100
PERCENT_DIFF_TO_BE_CALLED_AS_DMR = 0.4
CHANCE_OF_INCREASE_IN_METHYLATION = 0.9

# different runs for read depth, region size, percent difference


def get_regions() -> DataFrame:
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
    region: DataFrame,
    region_info: Series,
    rng: Generator,
) -> DataFrame:
    start = region_info.at["start_row"]
    end = region_info.at["end_row"]
    original_pm = region_info.at["original_pm"]
    inc_or_dec = region_info.at["inc_or_dec"]

    max_percent_diff = 1 - original_pm / 100
    min_percent_diff = PERCENT_DIFF_TO_BE_CALLED_AS_DMR
    goal_percent_diff = 100 * (
        rng.random() * (max_percent_diff - min_percent_diff) + min_percent_diff
    )

    iteration = 0

    while percent_diff(original_pm, start, end, region) < goal_percent_diff:
        selected_cytosine = rng.integers(start, end, endpoint=True)
        inc_or_dec_multiplier = 1 if inc_or_dec == "+" else -1

        min_delta = 0
        max_delta = (
            100 - region.at[selected_cytosine, "prop"]
            if inc_or_dec_multiplier == 1
            else region.at[selected_cytosine, "prop"]
        )
        delta = rng.random() * (max_delta - min_delta) + min_delta

        region.at[selected_cytosine, "prop"] += delta * inc_or_dec_multiplier

        # TODO: clean up variables here
        total_num_reads = (
            region.at[selected_cytosine, "uc"] + region.at[selected_cytosine, "mc"]
        )
        region.at[selected_cytosine, "mc"] = round(
            total_num_reads * region.at[selected_cytosine, "prop"] / 100
        )
        region.at[selected_cytosine, "uc"] = (
            total_num_reads - region.at[selected_cytosine, "mc"]
        )

        region.at[selected_cytosine, "prop"] = (
            100 * region.at[selected_cytosine, "mc"] / total_num_reads
        )

        # TODO: do this in a better way
        iteration += 1
        if iteration > MAX_REGION_SIZE * 10:
            break

    return region


def finalize_regions_info() -> None:
    for _, region in regions_info.iterrows():
        new_pm = bed_data.loc[
            region["start_row"] : region["end_row"] + 1, "prop"
        ].mean()
        regions_info.at[region.name, "new_pm"] = new_pm
        regions_info.at[region.name, "percent_diff"] = abs(
            new_pm - region["original_pm"]
        )


if __name__ == "__main__":
    # from dask.distributed import Client, LocalCluster

    # cluster = LocalCluster()
    # client = Client(cluster)

    # TODO: Make child random number generators

    col_labels = ["chr", "start", "end", "uc", "mc", "prop"]
    bed_data = pd.read_csv(BED_FILE, sep="\t", names=col_labels, header=None)

    regions_info = get_regions()

    global_rng = np.random.default_rng()

    bed_data = pd.concat(
        list(
            map(
                lambda region: simulate_dmr(
                    bed_data.iloc[region["start_row"] : region["end_row"] + 1],
                    region,
                    global_rng,
                ),
                [enum_region[1] for enum_region in regions_info.iterrows()],
            )
        )
    )

    finalize_regions_info()

    output_data_filename = os.path.join(OUT_DIR_DATA, "new_fp.bed")
    bed_data.to_csv(output_data_filename, sep="\t", index=False, header=False)

    output_region_filename = os.path.join(OUT_DIR_REGIONS, "new_fp_regions.tsv")
    regions_info.to_csv(output_region_filename, sep="\t", index=False, header=True)

    # x = client.submit(simulate_reads, 0.5, global_rng)
    # print(x.result())

    # futures = client.map(
    #     lambda region: get_simulated_region(
    #         bed_data.iloc[region["start_row"] : region["end_row"] + 1],
    #         region,
    #         global_rng,
    #     ),
    #     [enum_region[1] for enum_region in regions_info.iterrows()],
    # )

    # futures = []
    # for enum_region in regions_info.iterrows():
    #     futures.append(
    #         client.submit(
    #             get_simulated_region,
    #             bed_data.iloc[
    #                 enum_region[1]["start_row"] : enum_region[1]["end_row"] + 1
    #             ],
    #             enum_region[1],
    #             global_rng,
    #         )
    #     )

    # futures = client.map(
    #     lambda region: region,
    #     [enum_region[1] for enum_region in regions_info.iterrows()],
    # )

    # results = client.gather(futures)
    # print(results)

    # df = pd.concat(
    #     list(
    #         map(
    #             lambda region: get_simulated_region(
    #                 bed_data.iloc[region["start_row"] : region["end_row"] + 1],
    #                 region,
    #                 global_rng,
    #             ),
    #             [enum_region[1] for enum_region in regions_info.iterrows()],
    #         )
    #     )
    # )

    # finalize_regions_info()
    # print(df.iloc[5:15])

    # out_file = os.path.join(OUT_DIR, f"{os.path.basename(BED_FILE).replace('.bed', '')}_sample_{i}_ray.bed")

    # output_data_filename = os.path.join(OUT_DIR_DATA, "fp.bed")
    # bed_data.to_csv(output_data_filename, sep="\t", index=False, header=False)
    #
    # output_region_filename = os.path.join(OUT_DIR_REGIONS, "fp_regions.tsv")
    # regions_info.to_csv(output_region_filename, sep="\t", index=False, header=True)

    # client.close()
