"""
script to match medium to high quality MAGs with their domain/kingdoms, according to
pplacer placement via gtdbdk.

USAGE:
python3 /home/kchan_mb18/scripts/project_2/match_mag_taxa.py -a /home/kchan_mb18/SaanichInlet_200m/MetaBAT2_SaanichInlet_200m/gtdbtk_output/gtdbtk.ar122.classification_pplacer.tsv -b /home/kchan_mb18/SaanichInlet_200m/MetaBAT2_SaanichInlet_200m/gtdbtk_output/gtdbtk.bac120.classification_pplacer.tsv -m /home/kchan_mb18/SaanichInlet_200m/MetaBAT2_SaanichInlet_200m/MedQPlus_MAGs -o /home/kchan_mb18/p2
"""

import argparse
import os

parser = argparse.ArgumentParser(description = "Script to match medium to high quality MAGs with their domain/kingdoms, according to pplacer placement via gtdbdk.")
parser.add_argument("-a", "--archaea-placements", required = True, help = "Path to the pplacer archaeal classifications.")
parser.add_argument("-b", "--bacteria-placements", required = True, help = "Path to the pplacer bacterial classifications.")
parser.add_argument("-m", "--mags", required = True, help = "Path to the MAGs.")
parser.add_argument("-o", "--output-dir", default = ".", help = "Output directory. Default: current working directory.")
args = parser.parse_args()

archaea_placements = args.archaea_placements
bacteria_placements = args.bacteria_placements
mags = args.mags
output_dir = args.output_dir
id_taxa_map = {}

print("parsing archaeal classifications")
with open(archaea_placements, "rU") as f:
	for i, line in enumerate(f):
		bin_name, tax_string = line.strip().split("\t")
		top_level_taxa = tax_string.split(";")[0].replace("d__", "")
		id_taxa_map[bin_name] = top_level_taxa

print("parsing bacterial classifications")
with open(bacteria_placements, "rU") as f:
	for i, line in enumerate(f):
		bin_name, tax_string = line.strip().split("\t")
		top_level_taxa = tax_string.split(";")[0].replace("d__", "")
		id_taxa_map[bin_name] = top_level_taxa

print("writing tab delimited file")
outpath = os.path.join(output_dir, "id_taxa_map.txt")
outfile = open(outpath, "w")

for file in os.listdir(mags):
	if file.endswith(".fa"):
		bin_name = file.replace(".fa", "")
		outfile.write(bin_name + "\t" + id_taxa_map[bin_name] + "\n")
outfile.close()

print("done. goodbye")
