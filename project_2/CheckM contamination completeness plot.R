
# checkM contamination versus completeness

MetaBAT2_checkM <- MetaBAT2_SaanichInlet_200m_min1500_checkM_stdout

MetaBAT2_checkM %>%
  filter(Lineage != "root") %>%
ggplot(aes(x = Completeness, y = Contamination, colour = Lineage)) +
  geom_point(size = 3, alpha = 0.75) +
  theme_bw() +
  labs(x = "Completeness (%)", y = "Contamination (%)", title = "Visualization of checkM results",
       subtitle = "Saanich Inlet: Depth = 200m; Min = 1500") +
  annotate("rect", xmin = 90, xmax = 100, ymin = 0, ymax = 5, colour = "black",
             fill = NA) +
  scale_y_continuous(limits = c(0, 30))

# gtdbdk

gtdbtk_ar122 <- gtdbtk_ar122_classification_pplacer %>%
  select(Bins, L_1, L_2, L_3, L_4, L_5, L_6) %>%
  gather(key = "Lineage", value = "Taxonomies", -Bins)

ggplot(gtdbtk_ar122, aes(x = Bins, y = Taxonomies)) +
  geom_point() +
  coord_flip()

# phylum level

phylum_checkM <- read_tsv(file = "bin_complete_contam_tax_phylum.txt", col_names = T)
phylum_checkM_1 <- phylum_checkM %>% dplyr::rename(Phylum = `taxonomy (p)`)

ggplot(phylum_checkM_1, aes(x = completeness, y = contamination, colour = Phylum)) +
  geom_point(size = 3, alpha = 0.85) +
  annotate("rect", xmin = 90, xmax = 100, ymin = 0, ymax = 5, colour = "black",
           fill = NA)  +
  theme_bw() +
  theme(text = element_text(size = 14)) +
  labs(x = "Completeness (%)", y = "Contamination (%)") 

class_checkM <- read_tsv(file = "bin_complete_contam_tax_class.txt", col_names = T)

ggplot(class_checkM, aes(x = completeness, y = contamination, colour = `taxonomy (c)`)) +
  geom_point(size = 3, alpha = 0.85) +
  annotate("rect", xmin = 90, xmax = 100, ymin = 0, ymax = 5, colour = "black",
           fill = NA)

# Below will be the code for the final contamination-completeness plot
## It will incorporate the mean RPKM value for a given MAG in the form of 
#### A size element for each point. It will also include the bin identifier
##### to allow for easier referencing of the values for a particular bin
#
#
#
# cruise_RPKMs 

SI042_rpkm_data_2 <- read_csv(file = "SI042_200m.rpkm.csv", col_names = F) 
colnames(SI042_rpkm_data_2) <- c("MAG_orf", "SI042_rpkm")

SI048_rpkm_data_2 <- read_csv(file = "SI048_200m.rpkm.csv", col_names = F) 
colnames(SI048_rpkm_data_2) <- c("MAG_orf", "SI048_rpkm")

SI072_rpkm_data_2 <- read_csv(file = "SI072_200m.rpkm.csv", col_names = F) 
colnames(SI072_rpkm_data_2) <- c("MAG_orf", "SI072_rpkm")

SI073_rpkm_data_2 <- read_csv(file = "SI073_200m.rpkm.csv", col_names = F) 
colnames(SI073_rpkm_data_2) <- c("MAG_orf", "SI073_rpkm")

SI074_rpkm_data_2 <- read_csv(file = "SI074_200m.rpkm.csv", col_names = F) 
colnames(SI074_rpkm_data_2) <- c("MAG_orf", "SI074_rpkm")

SI075_rpkm_data_2 <- read_csv(file = "SI075_200m.rpkm.csv", col_names = F) 
colnames(SI075_rpkm_data_2) <- c("MAG_orf", "SI075_rpkm")

## Joining RPKMs into single one

Cruise_RPKM_joined <- left_join(SI042_rpkm_data_2, SI048_rpkm_data_2, by = "MAG_orf") %>%
  left_join(SI072_rpkm_data_2, by = "MAG_orf") %>%
  left_join(SI073_rpkm_data_2, by = "MAG_orf") %>%
  left_join(SI074_rpkm_data_2, by = "MAG_orf") %>%
  left_join(SI075_rpkm_data_2, by = "MAG_orf")

# join id_taxa_map to prokka_mag (based on )

id_taxa_map <- read_tsv(file = "id_taxa_map.txt", col_names = T)
prokka_to_MAG <- read_tsv(file = "prokka_to_mag_id.txt", col_names = T)

joined_taxa_prokka <- left_join(prokka_to_MAG, id_taxa_map, by = "MAG") 

df_joined_taxa_prokka <- as.data.frame(joined_taxa_prokka)

selected_joined_taxa_prokka <- df_joined_taxa_prokka %>%
  dplyr::select(Unique_ID, MAG, Domain)
colnames(selected_joined_taxa_prokka) <- c("MAG_ID", "Bin", "Domain")

# Map the joined cruise RPKMs to the MAG identifiers from prokka (and domain)

RPKM_prokka <- Cruise_RPKM_joined %>%
  separate(MAG_orf, into = c("MAG_ID", "orf_ID")) %>%
  left_join(selected_joined_taxa_prokka, by = "MAG_ID") %>%
  unite(MAG_orf, MAG_ID, orf_ID)

# Map RPKM prokka to KO based on MAG_orf

MAGs_ORFs_KO_2 <- read_tsv(file = "SaanichInlet_200m_all_MAGs_ORFs_ko.cleaned.txt", col_names = F)
colnames(MAGs_ORFs_KO_2) <- c("MAG_orf", "ko")

RPKM_prokka_KO <- left_join(MAGs_ORFs_KO_2, RPKM_prokka, by = "MAG_orf")

# taxonomy RPKM data

taxonomic_file_load <- read_tsv(file = "both_domains_classification_pplacer.txt", 
                                col_names = T)

taxonomies_all <- taxonomic_file_load %>%
  separate(Bin_ID, into = c("Name", "Bin"), sep = "m.") %>% 
  dplyr::rename(`Lowest classification` = Lowest_level_of_classification) %>%
  dplyr::select(-Name) %>%
  gather(key = "Taxonomic_classification", value = "RPKMs", -Bin, -Domain) %>%
  filter(!is.na(RPKMs)) 

taxonomies_all$Bin <- as.integer(taxonomies_all$Bin)

RPKM_all_taxonomies <- left_join(RPKM_prokka_KO, taxonomies_all,
                                     by = "Bin")

phylum_checkM <- read_tsv(file = "bin_complete_contam_tax_phylum.txt", col_names = T)

checkM_bins <- phylum_checkM %>%
  separate(bin_name, into = c("Name", "Bin"), sep = "m.") %>%
  dplyr::select(-Name, -`taxonomy (p)`)

checkM_bins$Bin <- as.integer(checkM_bins$Bin)

RPKM_checkM_all_taxonomies <- left_join(RPKM_all_taxonomies, checkM_bins, by = "Bin")

#Adding a size element to completeness contamination plot

bin_test <- RPKM_checkM_all_taxonomies %>%
  dplyr::select(-Domain.y) %>%
  filter(Taxonomic_classification == "Phylum") %>%
  gather(key = "Cruise", value = "RPKM", -MAG_orf, -ko, -Bin, -Domain.x) %>%
  group_by(Bin) 

bin_test$RPKM <- as.numeric(bin_test$RPKM)

bin_test_1 <- bin_test %>%
  filter(RPKM > 0) %>%
  filter(!is.na(RPKM)) %>%
  group_by(Bin) %>%
  summarize(count = n(),
            mean = mean(RPKM))

bin_test_plot <- left_join(bin_test_1, checkM_bins, by = "Bin")

taxonomies_class <- taxonomies_all %>%
  filter(Taxonomic_classification == "Phylum")

contam_complete_bin <- left_join(taxonomies_class, bin_test_plot, by = "Bin")

## Final contamination-completion plot with MAG mean RPKM as well as number
contam_complete_bin %>%
  separate(RPKMs, into = c("taxon", "Phylum"), sep= "__") %>%
ggplot(aes(x = completeness, y = contamination)) +
  geom_point(aes(size = mean, colour = Phylum), alpha = 0.5) +
  scale_size(range = c(1, 20), name = "Mean RPKM for\na given MAG") +
  geom_text(aes(label = Bin), size = 2) +
  annotate("rect", xmin = 90, xmax = 100, ymin = 0, ymax = 5, 
                fill = "transparent", colour = "black") +
  theme(text = element_text(size = 14), 
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = seq(0, 10, 2)) +
  labs(x = "Completeness (%)", y = "Contamination (%)") 
















