import pandas as pd
import os

in_dir = "./data-simulation/real-data/"
out_dir = "./data-simulation/scripts/sort-test/"

col_labels = ["chr", "start", "end", "uc", "mc", "prop"]


def sort_directory(in_dir, out_dir):
    files = os.listdir(in_dir)

    for file in files:
        if file.endswith(".bed"):
            bed_data = pd.read_csv(
                in_dir + file, sep="\t", names=col_labels, header=None
            )
            bed_data = bed_data.sort_values(by=["chr", "start"])
            bed_data.to_csv(out_dir + file, sep="\t", header=False, index=False)


def sort_single_file(file):
    bed_data = pd.read_csv(file, sep="\t", names=col_labels, header=None)
    bed_data = bed_data.sort_values(by=["chr", "start"])
    bed_data.to_csv(
        f"{os.path.basename(file).replace('.bed', '')}_sorted.bed",
        sep="\t",
        header=False,
        index=False,
    )


sort_single_file("../../real-data/D23_Col0_all_CpG.bed")

