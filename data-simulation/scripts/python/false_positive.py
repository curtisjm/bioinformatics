import os

os.environ["RAY_DEDUP_LOGS"] = "0"

from ray.util import ActorPool
from typing import Generator
from pandas import DataFrame
import ray
import pandas as pd
import numpy as np

# BED_FILE = "../../real-data/D23_Col0_all_CpG.bed"
BED_FILE = (
    "/Users/curtis/Documents/bioinformatics/data-simulation/real-data/10k_test.bed"
)
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


@ray.remote
class SimulationActor:
    def __init__(self, global_state: object):
        self.global_state = global_state

    def driver(self, first_region: int, last_region: int) -> None:
        for region_num in range(first_region, last_region + 1):
            # print("In driver loop")
            self.modification_handler(region_num)

    def modification_handler(self, region_num: int) -> None:
        regions = ray.get(self.global_state.get_regions.remote())

        start, end, original_pm, inc_or_dec = (
            regions.at[region_num, "start_row"],
            regions.at[region_num, "end_row"],
            regions.at[region_num, "original_pm"],
            regions.at[region_num, "inc_or_dec"],
        )
        if regions.at[region_num, "is_dmr"]:
            print(
                f"Producing DMR for region {start} to {end} with original pm {original_pm}"
            )
            new_pm = self.produce_dmr_iter_rand(start, end, original_pm, inc_or_dec)
        else:
            print(
                f"Simulating reads for region {start} to {end} with original pm {original_pm}"
            )
            new_pm = self.simulate_reads_for_region(start, end)
        print(f"\t New pm is {new_pm}")
        self.global_state.update_regions_entry.remote(region_num, "new_pm", new_pm)
        self.global_state.update_regions_entry.remote(
            region_num, "percent_diff", abs(new_pm - original_pm)
        )

    def produce_dmr_iter_rand(
        self, start: int, end: int, original_pm: float, inc_or_dec: str
    ) -> None:

        min_percent_diff = PERCENT_DIFF_TO_BE_CALLED_AS_DMR
        max_percent_diff = 1 - original_pm / 100
        goal_percent_diff = 100 * (
            ray.get(self.global_state.rand.remote())
            * (max_percent_diff - min_percent_diff)
            + min_percent_diff
        )
        # Select a random cytosine to change the methylation of and increase or decrease it's methlyation until the goal is reached
        percent_diff = self.percent_diff(original_pm, start, end)
        while percent_diff < goal_percent_diff:
            # print(f"In percent diff loop: PD: {percent_diff}, GPD: {goal_percent_diff}")

            # Select which cytosine we are modifying
            selected_cytosine = ray.get(
                self.global_state.rand_int.remote(start, end, True)
            )
            inc_or_dec_multiplier = 1 if inc_or_dec == "+" else -1

            # Determine how much to change the methylation of the cytosine
            min_delta = 0
            max_delta = 100 - ray.get(
                self.global_state.get_bed_data_entry.remote(selected_cytosine, "prop")
            )
            delta = (
                ray.get(self.global_state.rand.remote()) * (max_delta - min_delta)
                + min_delta
            )
            new_prop = (
                ray.get(
                    self.global_state.get_bed_data_entry.remote(
                        selected_cytosine, "prop"
                    )
                )
                + delta * inc_or_dec_multiplier
            )

            # Use the new proportion of methylated reads to calculate the number of methylated and unmethylated reads
            total_num_reads = ray.get(
                self.global_state.get_bed_data_entry.remote(selected_cytosine, "uc")
            ) + ray.get(
                self.global_state.get_bed_data_entry.remote(selected_cytosine, "mc")
            )
            new_mc = round(total_num_reads * new_prop / 100)
            new_uc = total_num_reads - new_mc
            self.global_state.update_bed_data_entry.remote(
                selected_cytosine, "mc", new_mc
            )
            self.global_state.update_bed_data_entry.remote(
                selected_cytosine, "uc", new_uc
            )

            # Update the proportion of methylated reads to reflected modified counts
            new_prop = 100 * new_mc / total_num_reads
            self.global_state.update_bed_data_entry.remote(
                selected_cytosine, "prop", new_prop
            )

            percent_diff = self.percent_diff(original_pm, start, end)

        return ray.get(self.global_state.get_average_methylation.remote(start, end))

    # For regions that are not DMRs, simulate the variation in reads
    def simulate_reads_for_region(self, start: int, end: int) -> float:
        bed_data = ray.get(self.global_state.get_bed_data.remote())
        new_pm = 0

        for row in range(start, end):
            # print("In simulate reads loop")
            uc_count, mc_count, sim_prop = self.simulate_reads(bed_data.at[row, "prop"])
            self.global_state.update_bed_data_entry.remote(row, "uc", uc_count)
            self.global_state.update_bed_data_entry.remote(row, "mc", mc_count)
            self.global_state.update_bed_data_entry.remote(row, "prop", sim_prop)

            new_pm += sim_prop

        return new_pm / (end - start)

    # Using a given proportion of methylation, simulate reads of each cytosine and
    # mutating the unmethylateted counts, methylated counts, and proportion of methylation
    # in the original data frame
    def simulate_reads(self, prop: float) -> tuple[int, int, float]:
        uc_count = 0
        mc_count = 0

        # Randomize the number of reads for each cytosine by adding a random number between
        # -READ_VARIATION * DEPTH and READ_VARIATION * DEPTH to the set depth
        num_reads = int(
            DEPTH
            + DEPTH
            * READ_VARIATION
            * (2 * ray.get(self.global_state.rand.remote()) - 1)
        )

        # Perform a weighted coin flip to determine if the cytosine is read as methylated or not by
        # generating a random number between 0 and 1 and checking if it is less than the true proportion of methylation
        # TODO: change this to normal distribution and find standard deviation
        random_values = 100 * ray.get(self.global_state.rand_array.remote(num_reads))
        mc_count = np.sum(random_values < prop)
        uc_count = num_reads - mc_count

        sim_prop = 100 * mc_count / (mc_count + uc_count)

        return (uc_count, mc_count, sim_prop)

    def percent_diff(self, original_pm: float, start: int, end: int) -> float:
        new_pm = ray.get(self.global_state.get_average_methylation.remote(start, end))
        # print(f"Original pm: {original_pm}, New pm: {new_pm}")
        return abs(new_pm - original_pm)


@ray.remote
class GlobalStateActor:
    def __init__(self) -> None:
        self.num_dmrs = 0
        # Load data from bed file into a pandas dataframe
        self.col_labels = ["chr", "start", "end", "uc", "mc", "prop"]
        self.bed_data = pd.read_csv(
            BED_FILE, sep="\t", names=self.col_labels, header=None
        )
        self.rng = np.random.default_rng()
        self.regions = self.define_regions()
        self.num_regions = self.regions.shape[0]

    # Divide the bed files into different regions
    def define_regions(self) -> DataFrame:
        # pm stands for percent methylation
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
        num_rows = self.bed_data.shape[0]

        # Calculate the probability of each region being a DMR
        estimated_num_regions = num_rows / ((MAX_REGION_SIZE - MIN_REGION_SIZE) / 2)
        chance_of_dmr = ESTIMATED_NUM_DMRS / estimated_num_regions

        while current_start < num_rows - 1:
            region_size = self.rng.integers(MIN_REGION_SIZE, MAX_REGION_SIZE)
            current_end = current_start + region_size

            # Make sure the region doesn't go beyond length of bed file
            if current_end >= self.bed_data.shape[0]:
                current_end = self.bed_data.shape[0] - 1

            # Make sure the region does not span multiple chromosomes
            while (
                self.bed_data.at[current_start, "chr"]
                != self.bed_data.at[current_end, "chr"]
            ):
                current_end -= 1

            # Calculate average percent methylation for the region
            region = self.bed_data.iloc[current_start:current_end]
            original_pm = region["prop"].mean()

            # Determine if the region is going to be made a DMR
            is_dmr = int(self.rng.random() < chance_of_dmr)
            self.num_dmrs += is_dmr

            # Occasionally decrease methylation instead of increasing it
            if is_dmr:
                inc_or_dec = (
                    "+"
                    if self.rng.random() < CHANCE_OF_INCREASE_IN_METHYLATION
                    else "-"
                )
                # Handle the edge case where the region can't have enough percent difference to hit threshold
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
                        self.bed_data.at[current_start, "start"],
                        self.bed_data.at[current_end, "end"],
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

    def get_num_dmrs(self) -> int:
        return self.num_dmrs

    def get_num_regions(self) -> int:
        return self.num_regions

    def get_bed_data(self) -> DataFrame:
        return self.bed_data

    def get_regions(self) -> DataFrame:
        return self.regions

    def get_rng(self) -> Generator:
        return self.rng

    def get_bed_data_entry(self, row: int, col: str) -> object:
        return self.bed_data.at[row, col]

    def get_regions_entry(self, row: int, col: str) -> object:
        return self.regions.at[row, col]

    def get_range_of_bed_data(self, start: int, end: int, col: str) -> DataFrame:
        return self.bed_data.loc[start:end, col]

    def get_range_of_regions(self, start: int, end: int, col: str) -> DataFrame:
        return self.regions.loc[start:end, col]

    def get_average_methylation(self, start: int, end: int) -> float:
        return self.bed_data.loc[start:end, "prop"].mean()

    def update_bed_data_entry(self, row: int, col: str, value) -> None:
        self.bed_data.at[row, col] = value

    def update_regions_entry(self, row: int, col: str, value) -> None:
        self.regions.at[row, col] = value

    def update_range_of_bed_data(
        self, start: int, end: int, col: str, new_column_segment: DataFrame
    ) -> None:
        self.bed_data.loc[start:end, col] = new_column_segment

    def update_range_of_regions(
        self, start: int, end: int, col: str, new_column_segment: DataFrame
    ) -> None:
        self.regions.loc[start:end, col] = new_column_segment

    def append_row_to_regions(self, new_row: DataFrame) -> None:
        return

    def rand(self) -> float:
        return self.rng.random()

    def rand_array(self, size: int) -> np.ndarray:
        return self.rng.random(size)

    def rand_int(self, start: int, end: int, end_inclusive: bool) -> int:
        return self.rng.integers(start, end, endpoint=end_inclusive)

    def bed_data_to_to_csv(self) -> None:
        # Output the simulated data to a new bed file
        # out_file = os.path.join(OUT_DIR, f"{os.path.basename(BED_FILE).replace('.bed', '')}_sample_{i}_ray.bed")
        output_data_filename = os.path.join(OUT_DIR_DATA, "false_pos_test.bed")
        self.bed_data.to_csv(output_data_filename, sep="\t", index=False, header=False)

    def regions_to_csv(self) -> None:
        output_region_filename = os.path.join(
            OUT_DIR_REGIONS, "false_pos_test_regions.tsv"
        )
        self.regions.to_csv(output_region_filename, sep="\t", index=False, header=True)


if __name__ == "__main__":
    ray.init()

    gs = GlobalStateActor.remote()
    num_regions = ray.get(gs.get_num_regions.remote())
    num_cores = int(ray.available_resources()["CPU"])
    num_regions_per_core = num_regions // num_cores
    region_ranges = [
        (i * num_regions_per_core, (i + 1) * num_regions_per_core - 1)
        for i in range(num_cores)
    ]
    if num_regions % num_cores:
        region_ranges[-1] = (region_ranges[-1][0], num_regions - 1)

    # pool = ActorPool([SimulationActor.remote(gs)])
    # region_ranges = [(0, num_regions - 1)]
    pool = ActorPool([SimulationActor.remote(gs) for _ in range(2 * num_cores)])
    list(
        pool.map(
            lambda actor, region_range: actor.driver.remote(*region_range),
            region_ranges,
        )
    )

    ray.get(gs.bed_data_to_to_csv.remote())
    ray.get(gs.regions_to_csv.remote())

    ray.shutdown()
