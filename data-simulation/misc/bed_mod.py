import pandas as pd

# col_labels = ["chr", "start", "end", "uc", "mc", "prop"]
# bed_data = pd.read_csv("test.tsv", sep="\t", names=col_labels, header=None)
# print(bed_data)

col_labels = ["chr", "start", "a", "b", "c", "context", "d"]
bed_data = pd.read_csv(
    "all_cytosines_with_cx.txt", sep="\t", names=col_labels, header=None
)
print(bed_data)

bed_data.drop(columns=["a", "b", "c", "d"], inplace=True)

bed_data.insert(2, "col", bed_data["start"] + 1)

bed_data.to_csv("./x.bed", sep="\t", index=False, header=False)
