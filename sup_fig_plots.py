import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import stylia as st
from stylia import NamedColors, ContinuousColorMap, ONE_COLUMN_WIDTH, TWO_COLUMNS_WIDTH
colors = NamedColors()
red = colors.get("red")
blue = colors.get("blue")
purple = colors.get("purple")
yellow = colors.get("yellow")
green = colors.get("green")
orange = colors.get("orange")
gray = colors.get("gray")

import csv
from tqdm import tqdm

file_name = 'subset_snvs_protein_coding_1kGPhg38.tsv'
col_names = ['EAS_AF','EUR_AF','AMR_AF','SAS_AF', "AFR_AF"]
with open(file_name, "r") as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader)
    idxs = [header.index(c) for c in col_names]
    R = []
    for r in tqdm(reader):
        r = [r[idx] for idx in idxs]
        R += [r]

data = pd.DataFrame(R, columns=col_names)
#data = pd.read_csv('subset_snvs_protein_coding_1kGPhg38.tsv', sep='\t', low_memory=False)
print(data.shape)
data = data.apply(pd.to_numeric, errors='coerce')
afr_abundant = data[data["AFR_AF"]>=0.2]
print(afr_abundant.shape)
alleleFreqs = ['EAS_AF','EUR_AF','AMR_AF','SAS_AF']
data['maxFreqs'] = data[alleleFreqs].max(axis=1)
data['AFR_overrepresentation'] = data['AFR_AF'].divide(data['maxFreqs'])
# Replace nan values by 0 (Ex. AF_AFR = 0 and maxFreqs = 0)
data['AFR_overrepresentation_mod'] = data['AFR_overrepresentation'].replace({np.nan:0, np.inf:50})
afr_specific = data[data['AFR_overrepresentation_mod']>=8]
print(afr_specific.shape)

fig, axs = st.create_figure(1,3, width=TWO_COLUMNS_WIDTH, height=TWO_COLUMNS_WIDTH*0.3)
ax = axs.next()
st.label(ax, xlabel="AFR AF", ylabel="Counts", title="")
ax.hist(data['AFR_AF'], bins=10, color=orange, alpha=0.7)
ax.set_xlim(-0.05,1.05)
ax = axs.next()
st.label(ax, xlabel="AFR AF", title="")
ax.set_xlim(-0.05,1.05)
ax.hist(afr_abundant['AFR_AF'], bins=10, color=yellow, alpha=0.7)
ax = axs.next()
st.label(ax, xlabel="Overrepresentation (AFR AF/ Max AF)", title="")
ax.hist(afr_specific['AFR_AF'], bins=10, color=green, alpha=0.7)
ax.set_xlim(-0.05,1.05)
plt.tight_layout()
st.save_figure("supf2.png")
st.save_figure("supf2.pdf")