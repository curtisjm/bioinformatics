import ray
import os
import pandas as pd
import numpy as np

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
NUM_DMRS = 10
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
def simulate_reads_for_region(start: int, end: int) -> float:
    new_pm = 0

    for row in range(start, end):
        uc_count, mc_count, sim_prop = simulate_reads(bed_data.at[row, "prop"])
        bed_data.at[row, "uc"] = uc_count
        bed_data.at[row, "mc"] = mc_count
        bed_data.at[row, "prop"] = sim_prop

        new_pm += sim_prop

    return new_pm / (end - start)

# Divide the bed files into different regions 
def define_regions() -> pd.DataFrame:
    # pm stands for percent methylation
    cols = ["start", "end", "original_pm", "new_pm", "is_dmr"]
    regions_df = pd.DataFrame(columns=cols)
    current_start = 0
    num_rows = bed_data.shape[0]

    # Calculate the probability of each region being a DMR
    estimated_num_regions = num_rows / ((MAX_REGION_SIZE - MIN_REGION_SIZE) / 2)
    chance_of_dmr = NUM_DMRS / estimated_num_regions

    while current_start < num_rows - 1:
        region_size = np.random.randint(MIN_REGION_SIZE, MAX_REGION_SIZE)
        current_end = current_start + region_size

        # Make sure the region doesn't go beyond length of bed file
        if current_end >= bed_data.shape[0]:
            current_end = bed_data.shape[0] - 1

        # Make sure the region does not span multiple chromosomes
        while bed_data.at[current_start, "chr"] != bed_data.at[current_end, "chr"]:
            current_end -= 1
        
        # Calculate average percent methylation for the region
        region = bed_data.iloc[current_start:current_end]
        pm = region["prop"].mean()

        # Determine if the region is going to be made a DMR
        is_dmr = int(np.random.rand() < chance_of_dmr)

        new_col = pd.DataFrame([[current_start, current_end, pm, 0, is_dmr]], columns=cols)
        regions_df = pd.concat([regions_df, new_col], ignore_index=True) if not regions_df.empty else new_col
        current_start = current_end
    return regions_df

#TODO: make this have randomness in where methlyation is changed
def produce_dmr_increase_all(start: int, end: int, original_pm: float) -> float:
    #TODO: should this be normal?
    percent_diff = np.random.uniform(PERCENT_DIFF_TO_BE_CALLED_AS_DMR, 1)
    if np.random.rand() < CHANCE_OF_INCREASE_IN_METHYLATION:
        new_pm = original_pm + percent_diff
    else:
        new_pm = original_pm - percent_diff

    bed_data.loc[start:end, "prop"] = bed_data.loc[start:end, "prop"] + bed_data.loc[start:end, "prop"].multiply(new_pm)
    #TODO: modify cytosine counts
    return new_pm

def produce_dmr_iter_rand(start: int, end: int, original_pm: float, new_pm: float) -> None:
    return

@ray.remote
def modification_handler(row: int) -> None:
    if regions.at[row, "is_dmr"]:
        regions.at[row, "new_pm"] = produce_dmr_increase_all(regions.at[row, "start"], regions.at[row, "end"], regions.at[row, "original_pm"])
    else:
       regions.at[row, "new_pm"] = simulate_reads_for_region(regions.at[row, "start"], regions.at[row, "end"])
        

# Initialize ray to manage parallel tasks
ray.init()

# Load data from bed file into a pandas dataframe
col_labels = ["chr", "start", "end", "uc", "mc", "prop"]
bed_data = pd.read_csv(BED_FILE, sep='\t', names=col_labels, header=None)

regions = define_regions()
futures = [modification_handler.remote(i) for i in range(regions.shape[0])]
ray.get(futures)

#TODO: set new_pm here?
#TODO: use uniform or normal distribution of percent diff?
#TODO: should is_dmr be set after this step to allow for 0 - 100 percent diff?

# Output the simulated data to a new bed file
# out_file = os.path.join(OUT_DIR, f"{os.path.basename(BED_FILE).replace('.bed', '')}_sample_{i}_ray.bed")
output_data_filename = os.path.join(OUT_DIR_DATA, "false_pos_test.bed")
bed_data.to_csv(output_data_filename, sep='\t', index=False, header=False)

output_region_filename = os.path.join(OUT_DIR_REGIONS, "false_pos_test_regions.tsv")
regions.to_csv(output_region_filename, sep='\t', index=False, header=True)
