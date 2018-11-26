# To whoever is marking this...

The following is a description of the contents of each file in this directory:

- `master.sh`: Shell script used for functional annotation with Prokka, aligning metatranscriptomic reads to binned contigs using BWA and Samtools, and RPKM calculations with the previously mentioned files.

- `match_mag_taxa.py`: A script to match medium to high quality MAGs with their domain/kingdoms, according to
pplacer taxonomic classification via gtdbdk. Used to loop over MAGs when running Prokka to match them with their respective kingdom/domains. 

- `match_prokka_mag.py`: A script to map prokka unique identifiers to MAG names (prefix name passed to prokka). Outputs a tab-delimited txt file. Used as an input to Pathview (in R).

- `get_complete_contam_tax.py`: A script that parses CheckM output for completeness and contamination estimates for each medium to high quality MAG, and combines that with the gtdbdk taxonomic classification for that MAG. Used to input into R to plot completeness vs contamination, with colors corresponding to taxonomic rank.

- `pathview_mag_match.py`: A script which maps KO number from a pathview diagram to Prokka IDs, and then Prokka IDs to bin names and their assigned taxonomies. Useful for profiling distributed metabolic pathways.

- `pie_chart_phyla.R`: R script which generates a pie chart, consisting of proportions of each MAG classified by taxonomy at the phylum level.
