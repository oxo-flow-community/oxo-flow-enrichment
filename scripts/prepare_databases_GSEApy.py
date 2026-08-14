#!/bin/env python
# Port of upstream workflow/scripts/prepare_databases_GSEApy.py.
# snakemake.input[0] / snakemake.output['db_file'] / snakemake.wildcards.database
# become positional CLI arguments <input> <output> <db>.
import json
import gseapy as gp
import os
import shutil
import sys

# input
db_path = sys.argv[1]

# output
results_path = sys.argv[2]

# parameters
db = sys.argv[3]


# if GMT, just copy
if db_path.lower().endswith('.gmt'):
    shutil.copy(db_path, results_path)
elif db_path.lower().endswith('.json'):
    # JSON load and save as GMT
    with open(db_path, 'r') as f:
        data = json.load(f)

    with open(results_path, 'w') as f:
        for key, values in data.items():
            f.write(f"{key}\t\t" + "\t".join(values) + "\n")
else:
    print("Error: Please provide a GMT (*.gmt) or JSON (*.json) database file.")
