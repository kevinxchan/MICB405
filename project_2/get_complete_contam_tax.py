"""
USAGE:

python3 scripts/project_2/get_complete_contam_tax.py -c SaanichInlet_200m/MetaBAT2_SaanichInlet_200m/MetaBAT2_SaanichInlet_200m_min1500_checkM_stdout.tsv -a SaanichInlet_200m/MetaBAT2_SaanichInlet_200m/gtdbtk_output/gtdbtk.ar122.classification_pplacer.tsv -b SaanichInlet_200m/MetaBAT2_SaanichInlet_200m/gtdbtk_output/gtdbtk.bac120.classification_pplacer.tsv -t p -o /home/kchan_mb18/
"""

import argparse
import os
import sys

parser = argparse.ArgumentParser(description = "Script to match MAG completion, contamination, and taxonomies from checkM and gtdbdk outputs.")
parser.add_argument("-c", "--checkm-tsv", required = True, help = "Path to the checkM tsv. This contains the completion and contamination estimates.")
parser.add_argument("-a", "--archaea-placements", required = True, help = "Path to the pplacer archaeal classifications.")
parser.add_argument("-b", "--bacteria-placements", required = True, help = "Path to the pplacer bacterial classifications.")
parser.add_argument("-t", "--taxonomy-level", required = True, help = "The desired taxonomy level. Must be one of [d, p, c, o, f, g, s].")
parser.add_argument("-o", "--output-dir", default = ".", help = "Output directory. Default: current working directory.")
args = parser.parse_args()

archaea_placements = args.archaea_placements
bacteria_placements = args.bacteria_placements
checkm_tsv = args.checkm_tsv
taxonomy_level = args.taxonomy_level
output_dir = args.output_dir

if taxonomy_level not in ("d", "p", "c", "o", "f", "g", "s"):
	print("ERROR: taxonomy level %s not recognized." % taxonomy_level)
	sys.exit(1)

bin_to_taxa_map = {}

with open(archaea_placements, "r") as f:
	for i, line in enumerate(f):
		data = line.split("\t")
		tax_string = data[1].split(";")
		for tax_level in tax_string:
			if tax_level.startswith(taxonomy_level + "__"):
				bin_to_taxa_map[data[0]] = tax_level.split("__")[1]

with open(bacteria_placements, "r") as f:
	for i, line in enumerate(f):
		data = line.split("\t")
		tax_string = data[1].split(";")
		for tax_level in tax_string:
			if tax_level.startswith(taxonomy_level + "__"):
				bin_to_taxa_map[data[0]] = tax_level.split("__")[1]

outfile = open(os.path.join(output_dir, "bin_complete_contam_tax.txt"), "w")
outfile.write("bin_name\tcompleteness\tcontamination\ttaxonomy (%s)" % (taxonomy_level) + "\n")
with open(checkm_tsv, "r") as f:
	for i, line in enumerate(f):
		data = line.split("\t")
		bin_id, completeness, contamination = data[0], data[11], data[12]
		if bin_id in bin_to_taxa_map:
			outfile.write("\t".join([bin_id, completeness, contamination, bin_to_taxa_map[bin_id]]) + "\n")

outfile.close()
