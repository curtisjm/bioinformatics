import ray
import os
import pandas as pd
import numpy as np

'''
Notes:
- TODO: Write description of what the script does.
- Does this script need batching in ray? https://docs.ray.io/en/latest/ray-core/patterns/too-fine-grained-tasks.html
- Is the large size of the object in the global state actor a problem? https://docs.ray.io/en/latest/ray-core/patterns/global-variables.html
- Should the rest of the functions be moved into their own actor?
    - Or into the global state actor?
- Make each task larger to reduce overhead of starting a new one.
'''

BED_FILE = "~/Documents/bioinformatics/data-simulation/real-data/curtis_testing.bed"
OUT_DIR_DATA = "~/Documents/bioinformatics/data-simulation/simulated-data"
OUT_DIR_REGIONS = "~/Documents/bioinformatics/data-simulation/simulated-data"
DEPTH = 25
NUM_SAMPLES = 1
#TODO: change this to standard deviation
STD_DEV = 0.15
READ_VARIATION = 0.15
# NUM_DMRS = 1000
# MIN_REGION_SIZE = 20
# MAX_REGION_SIZE = 3000
ESTIMATED_NUM_DMRS = 10
MIN_REGION_SIZE = 2
MAX_REGION_SIZE = 15
PERCENT_DIFF_TO_BE_CALLED_AS_DMR = 0.4
CHANCE_OF_INCREASE_IN_METHYLATION = 0.9

# Using a given proportion of methylation, simulate reads of each cytosine and
# mutating the unmethylateted counts, methylated counts, and proportion of methylation
# in the original data frame
def simulate_reads(prop: float) -> tuple[int, int, float]:
    np.random.seed()
    uc_count = 0
    mc_count = 0

    # Randomize the number of reads for each cytosine by adding a random number between
    # -READ_VARIATION * DEPTH and READ_VARIATION * DEPTH to the set depth
    num_reads = int(DEPTH + DEPTH * READ_VARIATION * (2 * np.random.rand() - 1))

    # Perform a weighted coin flip to determine if the cytosine is read as methylated or not by
    # generating a random number between 0 and 1 and checking if it is less than the true proportion of methylation
    #TODO: change this to normal distribution and find standard deviation
    random_values = 100 * np.random.rand(num_reads)
    mc_count = np.sum(random_values < prop)
    uc_count = num_reads - mc_count

    sim_prop = 100 * mc_count / (mc_count + uc_count)

    return (uc_count, mc_count, sim_prop)
    
# For regions that are not DMRs, simulate the variation in reads
def simulate_reads_for_region(global_state: object, start: int, end: int) -> float:
    new_pm = 0

    for row in range(start, end):
        uc_count, mc_count, sim_prop = simulate_reads(global_state.get_bed_data_entry.remote(row, "prop"))
        global_state.update_bed_date_entry.remote(row, "uc", uc_count)
        global_state.update_bed_date_entry.remote(row, "mc", mc_count)
        global_state.update_bed_date_entry.remote(row, "prop", sim_prop)

        new_pm += sim_prop

    return new_pm / (end - start)

#TODO: make this have randomness in where methlyation is changed
def produce_dmr_increase_all(global_state: object, start: int, end: int, original_pm: float) -> float:
    #TODO: should this be normal?
    #TODO: allow for decrease in methylation
    percent_diff = np.random.uniform(PERCENT_DIFF_TO_BE_CALLED_AS_DMR, 1)
    # Determine if the percent difference should be positive or negative
    if np.random.rand() > CHANCE_OF_INCREASE_IN_METHYLATION:
        percent_diff = -percent_diff
    new_pm = original_pm * (1 + percent_diff)
    # bed_data.loc[start:end, "prop"] = bed_data.loc[start:end, "prop"] + bed_data.loc[start:end, "prop"].multiply(new_pm)
    new_range_of_pm = global_state.get_range_of_bed_data.remote(start, end, "prop").mul(1 + percent_diff)
    global_state.update_range_of_bed_data.remote(start, end, "prop", new_range_of_pm)
    #TODO: modify cytosine counts
    return new_pm

def produce_dmr_iter_rand(start: int, end: int, original_pm: float, new_pm: float) -> None:
    return

@ray.remote
def modification_handler(global_state: object, region_num: int) -> None:
    start, end, original_pm = (global_state.get_regions_entry.remote(region_num, "start"),
                               global_state.get_regions_entry.remote(region_num, "end"),
                               global_state.get_regions_entry.remote(region_num, "original_pm"))
    if global_state.get_regions_entry.remote(region_num, "is_dmr"):
        print(f"Producing DMR for region {start} to {end} with original pm {original_pm}")
        new_pm = produce_dmr_increase_all(global_state, start, end, original_pm)
    else:
        print(f"Simulating reads for region {start} to {end} with original pm {original_pm}")
        new_pm = simulate_reads_for_region(global_state, start, end)
    print(f"\t New pm is {new_pm}")
    global_state.update_regions_entry.remote(region_num, "new_pm", new_pm)

@ray.remote
class GlobalStateActor():
    def __init__(self) -> None:
        self.num_dmrs = 0
        # Load data from bed file into a pandas dataframe
        self.col_labels = ["chr", "start", "end", "uc", "mc", "prop"]
        self.bed_data = pd.read_csv(BED_FILE, sep='\t', names=self.col_labels, header=None)

        self.regions = self.define_regions()
        self.num_regions = self.regions.shape[0]

    # Divide the bed files into different regions 
    def define_regions(self) -> pd.DataFrame:
        # pm stands for percent methylation
        cols = ["start", "end", "original_pm", "new_pm", "is_dmr"]
        regions_df = pd.DataFrame(columns=cols)
        current_start = 0
        num_rows = self.bed_data.shape[0]

        # Calculate the probability of each region being a DMR
        estimated_num_regions = num_rows / ((MAX_REGION_SIZE - MIN_REGION_SIZE) / 2)
        chance_of_dmr = ESTIMATED_NUM_DMRS / estimated_num_regions

        while current_start < num_rows - 1:
            region_size = np.random.randint(MIN_REGION_SIZE, MAX_REGION_SIZE)
            current_end = current_start + region_size

            # Make sure the region doesn't go beyond length of bed file
            if current_end >= self.bed_data.shape[0]:
                current_end = self.bed_data.shape[0] - 1

            # Make sure the region does not span multiple chromosomes
            while self.bed_data.at[current_start, "chr"] != self.bed_data.at[current_end, "chr"]:
                current_end -= 1

            # Calculate average percent methylation for the region
            region = self.bed_data.iloc[current_start:current_end]
            pm = region["prop"].mean()

            # Determine if the region is going to be made a DMR
            is_dmr = int(np.random.rand() < chance_of_dmr)
            self.num_dmrs += is_dmr

            new_col = pd.DataFrame([[current_start, current_end, pm, 0.0, is_dmr]], columns=cols)
            regions_df = pd.concat([regions_df, new_col], ignore_index=True) if not regions_df.empty else new_col
            current_start = current_end
        return regions_df

    def get_num_dmrs(self) -> int:
        return self.num_dmrs

    def get_num_regions(self) -> int:
        return self.num_regions

    def get_bed_data(self) -> pd.DataFrame:
        return self.bed_data
    
    def get_regions(self) -> pd.DataFrame:
        return self.regions
    
    def get_bed_data_entry(self, row: int, col: str) -> object:
        return self.bed_data.at[row, col]
    
    def get_regions_entry(self, row: int, col: str) -> object:
        return self.regions.at[row, col]
    
    def get_range_of_bed_data(self, start: int, end: int, col: str) -> pd.DataFrame:
        return self.bed_data.loc[start:end, col]

    def get_range_of_regions(self, start: int, end: int, col: str) -> pd.DataFrame:
        return self.regions.loc[start:end, col]
    
    def update_bed_date_entry(self, row: int, col: str, value) -> None:
        self.bed_data.at[row, col] = value

    def update_regions_entry(self, row: int, col: str, value) -> None:
        self.regions.at[row, col] = value

    #TODO: check functionality of these methods
    def update_range_of_bed_data(self, start: int, end: int, col: str, new_column_segment: pd.DataFrame) -> None:
        self.bed_data.loc[start:end, col] = new_column_segment
    
    def update_range_of_regions(self, start: int, end: int, col: str, new_column_segment: pd.DataFrame) -> None:
        self.regions.loc[start:end, col] = new_column_segment
    
    def bed_data_to_to_csv(self) -> None:
        # Output the simulated data to a new bed file
        # out_file = os.path.join(OUT_DIR, f"{os.path.basename(BED_FILE).replace('.bed', '')}_sample_{i}_ray.bed")
        output_data_filename = os.path.join(OUT_DIR_DATA, "false_pos_test.bed")
        self.bed_data.to_csv(output_data_filename, sep='\t', index=False, header=False)
    
    def regions_to_csv(self) -> None:
        output_region_filename = os.path.join(OUT_DIR_REGIONS, "false_pos_test_regions.tsv")
        self.regions.to_csv(output_region_filename, sep='\t', index=False, header=True)
       
       
       
       
# Initialize ray to manage parallel tasks
ray.init()

global_state_actor = GlobalStateActor.remote()

futures = [modification_handler.remote(global_state_actor, region_num)
           for region_num in range(ray.get(global_state_actor.get_num_regions.remote()))]
ray.get(futures)

#TODO: set new_pm here?
#TODO: use uniform or normal distribution of percent diff?
#TODO: should is_dmr be set after this step to allow for 0 - 100 percent diff?

#TODO: which version to use?
# global_state_actor.bed_data_to_to_csv.remote()
# global_state_actor.regions_to_csv.remote()
# ray.get(global_state_actor.bed_data_to_to_csv.remote())
# ray.get(global_state_actor.regions_to_csv.remote())

# Output the simulated data to a new bed file
# out_file = os.path.join(OUT_DIR, f"{os.path.basename(BED_FILE).replace('.bed', '')}_sample_{i}_ray.bed")
output_data_filename = os.path.join(OUT_DIR_DATA, "false_pos_test.bed")
ray.get(global_state_actor.get_bed_data.remote()).to_csv(output_data_filename, sep='\t', index=False, header=False)

output_region_filename = os.path.join(OUT_DIR_REGIONS, "false_pos_test_regions.tsv")
ray.get(global_state_actor.get_regions.remote()).to_csv(output_region_filename, sep='\t', index=False, header=True)
