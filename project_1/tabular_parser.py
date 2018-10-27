
import argparse
import pandas as pd
import numpy as np
import os

parser = argparse.ArgumentParser(description="Get number of variants from tabular file for each sample.")
parser.add_argument("-t",
					"--tab-dir",
					type=str,
					required=True,
					help="Absolute path to directory containing tabular files files.")
parser.add_argument("-o", 
                    "--out-path",
                    type=str,
                    default=".",
                    help="Path to output directory. Default: current working directory.")
args = parser.parse_args()

tab_dir = args.tab_dir
outpath = os.path.join(args.out_path, "summary.txt")
outfile = open(outpath, "w")

for f in os.listdir(tab_dir):
	if f.endswith(".tabular"):
		with open(os.path.join(tab_dir, f), "r") as infile:
			df = pd.read_table(infile, sep = "\t")
			for col in df:
				if "quality" in col:
					outfile.write(col + "\t" + str(sum(np.where(df[col] != "-", 1, 0))) + "\n")
outfile.close()
