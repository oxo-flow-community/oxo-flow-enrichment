# libraries
import os
import sys
import pandas as pd
import pycistarget
from pycistarget.input_output import read_hdf5

# Port of upstream workflow/scripts/process_results_pycisTarget.py.
# snakemake.input['motif_hdf5'], snakemake.output['motif_csv'|'hits_csv'|
# 'cistrome_csv'], snakemake.wildcards['region_set'] and
# snakemake.config["pycistarget_parameters"]["annotations_to_use"][0] become
# positional CLI arguments <motif_hdf5> <motif_csv> <hits_csv> <cistrome_csv>
# <region_set_name> <term_col>.

# input
motif_hdf5_path = sys.argv[1]

# output
motif_csv_path = sys.argv[2]
hits_csv_path = sys.argv[3]
cistrome_csv_path = sys.argv[4]
# parameters
region_set_name = sys.argv[5]
term_col = sys.argv[6]

# quit early if file is empty
if os.path.getsize(motif_hdf5_path) == 0:
    open(motif_csv_path, 'w').close()
    open(hits_csv_path, 'w').close()
    open(cistrome_csv_path, 'w').close()
    quit()

# load pycisTarget results from hdf5
results = read_hdf5(motif_hdf5_path)

# extract results
results_df = results[region_set_name].motif_enrichment

# reformat
results_df.index.name = "motif"
results_df.reset_index(inplace=True)
results_df["description"] = results_df["motif"] + "(" + results_df[term_col] + ")"

# save motif enrichments as CSV for downstream processing and plotting
results_df.to_csv(motif_csv_path)

# Save motif hits
motifs_df = pd.DataFrame(results[region_set_name].motif_hits)
motifs_df["database"] = motifs_df["database"].apply(
    lambda x: ";".join(x) if type(x) == list else ""
)
motifs_df["region_set"] = motifs_df["region_set"].apply(
    lambda x: ";".join(x) if type(x) == list else ""
)
motifs_df.to_csv(hits_csv_path)

# Save cistrome hits
cistrome_df = pd.DataFrame(results[region_set_name].cistromes)
cistrome_df["database"] = cistrome_df["database"].apply(
    lambda x: ";".join(x) if type(x) == list else ""
)
cistrome_df["region_set"] = cistrome_df["region_set"].apply(
    lambda x: ";".join(x) if type(x) == list else ""
)
cistrome_df.to_csv(cistrome_csv_path)
