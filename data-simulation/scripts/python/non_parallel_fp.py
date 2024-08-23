import os
import pandas as pd
import numpy as np

# BED_FILE = "../../real-data/100_test.bed"
BED_FILE = "/Users/curtis/Documents/bioinformatics/data-simulation/real-data/100_test.bed"
OUT_DIR_DATA = "./"
OUT_DIR_REGIONS = "./"
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
    uc_count = 0
    mc_count = 0

    # Randomize the number of reads for each cytosine by adding a random number between
    # -READ_VARIATION * DEPTH and READ_VARIATION * DEPTH to the set depth
    num_reads = int(DEPTH + DEPTH * READ_VARIATION * (2 * rng.random() - 1))

    # Perform a weighted coin flip to determine if the cytosine is read as methylated or not by
    # generating a random number between 0 and 1 and checking if it is less than the true proportion of methylation
    #TODO: change this to normal distribution and find standard deviation
    random_values = 100 * rng.random(num_reads)
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


def percent_diff(original_pm: float, start: int, end: int) -> float:
    new_pm = bed_data.loc[start:end, "prop"].mean()
    return abs(new_pm - original_pm)

def produce_dmr_iter_rand(start: int, end: int, original_pm: float) -> None:
    goal_percent_diff = 100 * rng.random() * PERCENT_DIFF_TO_BE_CALLED_AS_DMR
    print(bed_data.loc[start:end])
    # Select a random cytosine to change the methylation of and increase or decrease it's methlyation until the goal is reached
    while percent_diff(original_pm, start, end) < goal_percent_diff:
        # print(f"Percent diff is {percent_diff(original_pm, start, end)}")
        # Select which cytosine we are modifying
        selected_cytosine = rng.integers(start, end)
        
        # Determine how much to change the methylation of the cytosine
        min_delta = 0
        max_delta = 100 - bed_data.at[selected_cytosine, "prop"]
        delta = rng.random() * (max_delta - min_delta) + min_delta
        bed_data.at[selected_cytosine, "prop"] += delta

        # Use the new proportion of methylated reads to calculate the number of methylated and unmethylated reads
        total_num_reads = bed_data.at[selected_cytosine, "uc"] + bed_data.at[selected_cytosine, "mc"]
        bed_data.at[selected_cytosine, "mc"] = round(total_num_reads * bed_data.at[selected_cytosine, "prop"] / 100)
        bed_data.at[selected_cytosine, "uc"] = total_num_reads - bed_data.at[selected_cytosine, "mc"] 

        # Update the proportion of methylated reads to reflected modified counts
        bed_data.at[selected_cytosine, "prop"] = 100 * bed_data.at[selected_cytosine, "mc"] / total_num_reads
        
    return bed_data.loc[start:end, "prop"].mean()

def modification_handler(region_num: int) -> None:
    start, end, original_pm = (regions.at[region_num, "start_row"], regions.at[region_num, "end_row"], regions.at[region_num, "original_pm"])
    if regions.at[region_num, "is_dmr"]:
        print(f"Producing DMR for region {start} to {end} with original pm {original_pm}")
        new_pm = produce_dmr_iter_rand(start, end, original_pm)
    else:
        print(f"Simulating reads for region {start} to {end} with original pm {original_pm}")
        new_pm = simulate_reads_for_region(start, end)
    print(f"\t New pm is {new_pm}")
    regions.at[region_num, "new_pm"] = new_pm

# Divide the bed files into different regions 
def define_regions() -> pd.DataFrame:
    num_dmrs = 0
    # pm stands for percent methylation
    cols = ["start_row", "end_row", "start_coord", "end_coord", "original_pm", "new_pm", "percent_diff", "is_dmr", "inc_or_dec"]
    regions_df = pd.DataFrame(columns=cols)
    current_start = 0
    num_rows = bed_data.shape[0]

    # Calculate the probability of each region being a DMR
    estimated_num_regions = num_rows / ((MAX_REGION_SIZE - MIN_REGION_SIZE) / 2)
    chance_of_dmr = ESTIMATED_NUM_DMRS / estimated_num_regions

    while current_start < num_rows - 1:
        region_size = rng.integers(MIN_REGION_SIZE, MAX_REGION_SIZE)
        current_end = current_start + region_size

        # Make sure the region doesn't go beyond length of bed file
        if current_end >= bed_data.shape[0]:
            current_end = bed_data.shape[0] - 1

        # Make sure the region does not span multiple chromosomes
        while bed_data.at[current_start, "chr"] != bed_data.at[current_end, "chr"]:
            current_end -= 1

        # Calculate average percent methylation for the region
        region = bed_data.iloc[current_start:current_end]
        original_pm = region["prop"].mean()

        # Determine if the region is going to be made a DMR
        is_dmr = int(rng.random() < chance_of_dmr)
        num_dmrs += is_dmr

        if is_dmr:
            inc_or_dec = "+" if rng.random() < CHANCE_OF_INCREASE_IN_METHYLATION else "-"
        else:
            inc_or_dec = "N/A"

        new_row = pd.DataFrame([[current_start, current_end, bed_data.at[current_start, "start"],
                                 bed_data.at[current_end, "end"], original_pm, 0.0, 0.0, is_dmr, inc_or_dec]], columns=cols)
        regions_df = pd.concat([regions_df, new_row], ignore_index=True) if not regions_df.empty else new_row
        current_start = current_end
    return regions_df
        
col_labels = ["chr", "start", "end", "uc", "mc", "prop"]
bed_data = pd.read_csv(BED_FILE, sep='\t', names=col_labels, header=None)

rng = np.random.default_rng()

regions = define_regions()

for region_num in range(regions.shape[0]):
    modification_handler(region_num)

# Output the simulated data to a new bed file
# out_file = os.path.join(OUT_DIR, f"{os.path.basename(BED_FILE).replace('.bed', '')}_sample_{i}_ray.bed")
output_data_filename = os.path.join(OUT_DIR_DATA, "false_pos_test.bed")
bed_data.to_csv(output_data_filename, sep='\t', index=False, header=False)

output_region_filename = os.path.join(OUT_DIR_REGIONS, "false_pos_test_regions.tsv")
regions.to_csv(output_region_filename, sep='\t', index=False, header=True)
