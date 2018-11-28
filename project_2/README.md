# To whoever is marking this...

The following is a description of the contents of each file in this directory:

- `master.sh`: Shell script used for functional annotation with Prokka, aligning metatranscriptomic reads to binned contigs using BWA and Samtools, and RPKM calculations with the previously mentioned files.

- `match_mag_taxa.py`: A script to match medium to high quality MAGs with their domain/kingdoms, according to
pplacer taxonomic classification via gtdbdk. Used to loop over MAGs when running Prokka to match them with their respective kingdom/domains. 

- `match_prokka_mag.py`: A script to map prokka unique identifiers to MAG names (prefix name passed to prokka). Outputs a tab-delimited txt file. Used as an input to Pathview (in R).

- `get_complete_contam_tax.py`: A script that parses CheckM output for completeness and contamination estimates for each medium to high quality MAG, and combines that with the gtdbdk taxonomic classification for that MAG. Used to input into R to plot completeness vs contamination, with colors corresponding to taxonomic rank.

- `pathview_mag_match.py`: A script which maps KO number from a pathview diagram to Prokka IDs, and then Prokka IDs to bin names and their assigned taxonomies. Useful for profiling distributed metabolic pathways.

- `pie_chart_phyla.R`: R script which generates a pie chart, consisting of proportions of each MAG classified by taxonomy at the phylum level.

- `CheckM contamination completeness plot.R`: which contains an R script for a plot that visualizes the output of CheckM quality assessments. Only the medium- to high-quality MAGs were included in this analysis. In addition, an abundance layer representing the mean RPKM values for a given MAG is integrated into the plot. 

- `Saanich Inlet Chemical and environmental gradietns.R`: R script which generates two plots, one for the chemical gradients sampled at different depths and the other for the environmental parameters of Saanich Inlet at different depths.

- `Density scatterplot transformation of RPKMs.R`: R script providing rationale for the logarithmic transformation of RPKM values prior to pathview analysis. A density scatterplot of RPKM values regardless of KO identifier before and after it has undergone a logarithmic transformation.

- `Pathview grouping KOs by taxonomy.R` R script for the data augmentation in preparation for the pathview package visualization of metabolic pathways. 

- `RPKM sulfur bubble plot.R`: R script for the visualization of transcriptional activity of the microbial community at the class taxonomic level against all the genes in the sulfur metabolism KEGG pathway. The transcriptional activity of a given microbe is provided as a function of point size, which represents the Log summed RPKM values.  
