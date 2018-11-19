"""
script to map prokka unique identifiers to MAG names (prefix name passed to prokka). outputs a tsv like so:

PROKKA_ID_1 BIN_NAME_1
PROKKA_ID_2 BIN_NAME_2
PROKKA_ID_3 BIN_NAME_3
PROKKA_ID_4 BIN_NAME_4
		...

USAGE:
python3 ~/scripts/project_2/match_prokka_mag.py -p /home/kchan_mb18/p2/prokka_output -o .
"""


import argparse
import os

parser = argparse.ArgumentParser(description = "Script to map prokka unique identifiers to MAG names (prefix name passed to prokka).")
parser.add_argument("-p", "--prokka-dir", required = True, help = "Path to the pplacer archaeal classifications.")
parser.add_argument("-o", "--output-dir", default = ".", help = "Output directory. Default: current working directory.")
args = parser.parse_args()

prokka_dir = args.prokka_dir
outdir = args.output_dir
prokka_bin_map = {}

for d in os.listdir(prokka_dir):
	if os.path.isdir(d):
		for f in os.listdir(d):
			if f.endswith("tsv"):
				with open(os.path.join(prokka_dir, d, f), "r") as infile:
					prokka_bin_map[infile.readlines()[1].strip().split("\t")[0].split("_")[0]] = os.path.splitext(f)[0]

outfile = open(os.path.join(outdir, "prokka_to_mag_id.txt"), "w")

for k, bin in prokka_bin_map.items():
	outfile.write("{}\n".format("\t".join([k, bin])))
outfile.close()