# This doc will be looking to group KO numbers via some given taxonomic rank
## While the lower level taxonomies are ideal, they are less complete;
###lowest level of classification for each bin would be interesting, but runs into same 
####problem as having 50 MAGs (50 unqiue lowest levels or so), and trying to input 
##### that into pathview()

# Same loading procedure for RPKMs/MAGs/KO

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

## Joining RPKMs into single df

Cruise_RPKM_joined <- left_join(SI042_rpkm_data_2, SI048_rpkm_data_2, by = "MAG_orf") %>%
  left_join(SI072_rpkm_data_2, by = "MAG_orf") %>%
  left_join(SI073_rpkm_data_2, by = "MAG_orf") %>%
  left_join(SI074_rpkm_data_2, by = "MAG_orf") %>%
  left_join(SI075_rpkm_data_2, by = "MAG_orf")

# join id_taxa_map to prokka_mag (based on MAG)

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

# File with the classifications of both domains (gtdbdck or whatever it's called)
#for medium to high completeness-contamination bins
# This was constructed by combining the gtbdtk files for archaea and bacteria  

## This one will be for order level

taxonomic_file_load <- read_tsv(file = "both_domains_classification_pplacer.txt", col_names = T) 
# ### This file is just the arc and bac pplacer classifications combined

taxonomic_classifications_order <- taxonomic_file_load %>%
  separate(Bin_ID, into = c("Name", "Bin"), sep = "m.") %>% 
  filter(!is.na(Order)) %>%
  unite(tax_order, Order, Domain, sep = " ") %>%
  dplyr::select(Bin, tax_order) %>%
  separate(tax_order, into = c("o", "Order"), sep = "__") %>%
  dplyr::select(-o)

taxonomic_classifications_order$Bin <- as.integer(taxonomic_classifications_order$Bin)
           
RPKM_prokka_KO_taxonomy <- left_join(RPKM_prokka_KO, taxonomic_classifications_order,
                                   by = "Bin")

# Data formatting for pathview (as per Connor's pathview guide) ORDER

t_RPKM_2 <- RPKM_prokka_KO_taxonomy %>%
  separate(MAG_orf, into=c("prokka_id", "orf_id")) %>%
  dplyr::select(-Bin, -Domain, -prokka_id) %>%
  gather(key = "Cruise", value = "RPKM", -Order, -orf_id, -ko) %>%
  group_by(Order, ko) %>%
  summarize(total = log10(sum(RPKM))) %>%
  filter(total != -Inf) %>%
  spread(key = Order, value = total)

t_RPKM_2_pv_mat <- dplyr::select(t_RPKM_2, -ko)
rownames(t_RPKM_2_pv_mat) <- t_RPKM_2$ko 

## 31 variables... :( as per Luo Weijin (?) - pathview developer - more than 25 samples 
##may overload the pathview graphing capabilities. The pathview figs, KEGG native (doesn't apply to Graphviz)
##is a raster image. So, each rectangle (gene) node span's a fixed number of lines.
##In other words, they cannot be divided into more than 25 pieces - 25 representing an
##upper limit

## Pathview time
library(pathview)

## pathview for orders even though it won't work 

t_RPKM_2_n.pv <- pathview(gene.data = t_RPKM_2_pv_mat,
                          species = "ko",
                          pathway.id="00910",
                          out.suffix = "ko.data_order", 
                          gene.idtype = "KEGG",
                          limit = 5,
                          na.col = "lightgray",
                          # Graphviz still works, but ugly
                          kegg.native = F)
t_RPKM_2_order_df <- t_RPKM_2_n.pv$plot.data.gene
t_RPKM_2.1_order_df <- t_RPKM_2_n.pv$plot.data.cpd

### Looking at the sulfur metabolic pathway

t_RPKM_2_sulfur.pv <- pathview(gene.data = t_RPKM_2_pv_mat,
                          species = "ko",
                          pathway.id="00920",
                          out.suffix = "ko.data_order_sulfur", 
                          gene.idtype = "KEGG",
                          limit = c(-1, 3),
                          # Graphviz still works, but ugly
                          kegg.native = F)

t_RPKM_2_order_sulfur_df <- t_RPKM_2_sulfur.pv$plot.data.gene
t_RPKM_2.1_order_sulfur_df <- t_RPKM_2_sulfur.pv$plot.data.cpd

# The final pathview figures will be according to the class level, as it is relatively
# while being the lowest level that does not overload pathview function
#
# Note: R and pathview package are updated to current versions.
#
#
# Grouping by class 

taxonomic_classifications_class <- taxonomic_file_load %>%
  separate(Bin_ID, into = c("Name", "Bin"), sep = "m.") %>% 
  filter(!is.na(Class)) %>%
  unite(tax_order, Class, Domain, sep = " ") %>%
  dplyr::select(Bin, tax_order) %>%
  separate(tax_order, into = c("c", "Class"), sep = "__") %>%
  dplyr::select(-c)

taxonomic_classifications_class$Bin <- as.integer(taxonomic_classifications_class$Bin)

RPKM_prokka_KO_tax_class <- left_join(RPKM_prokka_KO, taxonomic_classifications_class,
                                     by = "Bin")

# Data formatting for pathview (as per Connor's pathview guide) CLASS

t_RPKM_3 <- RPKM_prokka_KO_tax_class %>%
  separate(MAG_orf, into=c("prokka_id", "orf_id")) %>%
  dplyr::select(-Bin, -Domain, -prokka_id) %>%
  gather(key = "Cruise", value = "RPKM", -Class, -orf_id, -ko) %>%
  group_by(Class, ko) %>%
  summarize(total = log10(sum(RPKM))) %>%
  filter(total != -Inf) %>%
  spread(key = Class, value = total)

t_RPKM_3_pv_mat <- dplyr::select(t_RPKM_3, -ko)
rownames(t_RPKM_3_pv_mat) <- t_RPKM_3$ko 

###NOTE: log10() transformation will introduce negative values. So, the final pathview plots
#### will visualize some of the ORFs along the negative gradient

#Also making an untransformed version
t_RPKM_4 <- RPKM_prokka_KO_tax_class %>%
  separate(MAG_orf, into=c("prokka_id", "orf_id")) %>%
  dplyr::select(-Bin, -Domain, -prokka_id) %>%
  gather(key = "Cruise", value = "RPKM", -Class, -orf_id, -ko) %>%
  group_by(Class, ko) %>%
  summarize(total = sum(RPKM)) %>%
  spread(key = Class, value = total)

t_RPKM_4_pv_mat <- dplyr::select(t_RPKM_4, -ko)
rownames(t_RPKM_4_pv_mat) <- t_RPKM_4$ko 

#write.csv(t_RPKM_3_pv_mat, file = "t_RPKM_3_pv_mat.csv")

## pathview for class 

### nitrogen - 00910

t_RPKM_3_n.pv <- pathview(gene.data = t_RPKM_3_pv_mat,
                          species = "ko",
                          pathway.id="00910",
                          out.suffix = "nitrogen", 
                          gene.idtype = "KEGG",
                          limit = list(gene = c(-1,3)),
                          low = list(gene = "#fdbb2d"), 
                          mid = list(gene = "#b21f1f"),
                          high = list(gene = "#1a2a6c"),
                          kegg.native = T,
                          same.layer= T)

t_RPKM_3_n_layer.pv <- pathview(gene.data = t_RPKM_3_pv_mat,
                          species = "ko",
                          pathway.id="00910",
                          out.suffix = "nitrogen_layer", 
                          gene.idtype = "KEGG",
                          limit = list(gene = c(-1,3)),
                          low = list(gene = "#fdbb2d"), 
                          mid = list(gene = "#b21f1f"),
                          high = list(gene = "#1a2a6c"),
                          kegg.native = T,
                          same.layer= F)


t_RPKM_3_class_df <- t_RPKM_3_n.pv$plot.data.gene
t_RPKM_3.1_class_df <- t_RPKM_3_n.pv$plot.data.cpd


### Sulphur 00920
#### log transformed 

t_RPKM_3_s.pv <- pathview(gene.data = t_RPKM_3_pv_mat,
                          species = "ko",
                          pathway.id="00920",
                          out.suffix = "number of different pathways class", 
                          gene.idtype = "KEGG",
                          limit = list(gene = c(-1,3)),
                          low = list(gene = "#fdbb2d"), 
                          mid = list(gene = "#b21f1f"),
                          high = list(gene = "#1a2a6c"),
                          kegg.native = T,
                          same.layer= T)

sulfur_class_df_gene <- t_RPKM_3_s.pv$plot.data.gene

sulfur_df <- as.data.frame(sulfur_class_df_gene) %>%
  dplyr::select(kegg.names:`<NA>`)

sulfur_df_xls <- sulfur_data_frame_coerced %>% 
  dplyr::select(kegg.names:`<NA>`) %>%
  dplyr::select(-type, -x, -y, -width, -height) %>%
  filter(!is.na(all.mapped)) %>%
  melt(id = c("kegg.names", "labels", "all.mapped")) %>%
  filter(value != "NA") %>%
  dplyr::rename(Class = variable) %>%
  group_by(kegg.names) %>%
  arrange(kegg.names)

KOs_for_join_ <- data.frame(unique(sulfur_df_xls$labels))

sulfur_df_xls$value <- as.numeric(sulfur_df_xls$value)

sulfur_class_df_cpd <- t_RPKM_3_s.pv$plot.data.cpd

####Sulfur non-transformed
t_RPKM_4_s.pv <- pathview(gene.data = t_RPKM_4_pv_mat,
                          species = "ko",
                          pathway.id="00920",
                          out.suffix = "number of different pathways class", 
                          gene.idtype = "KEGG",
                          low = list(gene = "#fdbb2d"), 
                          mid = list(gene = "#b21f1f"),
                          high = list(gene = "#1a2a6c"),
                          kegg.native = T,
                          same.layer= T)

nontransformed_sulfur_class_df_gene <- t_RPKM_4_s.pv$plot.data.gene

sulfur_df_normal_xls <- Sulfur_pathview_gene_df_coerced %>% 
  dplyr::select(kegg.names:`<NA>`) %>%
  dplyr::select(-type, -x, -y, -width, -height) %>%
  filter(!is.na(all.mapped)) %>%
  melt(id = c("kegg.names", "labels", "all.mapped")) %>%
  filter(value != "NA") %>%
  dplyr::rename(Class = variable) %>%
  group_by(kegg.names) %>%
  arrange(kegg.names)


t_RPKM_3_s_layer.pv <- pathview(gene.data = t_RPKM_3_pv_mat,
                          species = "ko",
                          pathway.id="00920",
                          out.suffix = "number of different pathways class layer", 
                          gene.idtype = "KEGG",
                          limit = list(gene = c(-1,3)),
                          low = list(gene = "#fdbb2d"), 
                          mid = list(gene = "#b21f1f"),
                          high = list(gene = "#1a2a6c"),
                          kegg.native = T,
                          same.layer= F)


sulfur_class_df_gene_layer <- t_RPKM_3_s_layer.pv$plot.data.gene
sulfur_class_df_cpd_layer <- t_RPKM_3_s_layer.pv$plot.data.cpd

#write.csv(sulfur_df_xls, file = "sulfur_metabolic_pathway_KO_mapping.txt")


## Below are other potential pathways of interest (or not)
##
##
##
### Glycolysis just as a random pathway 00010

t_RPKM_3_gylc.pv <- pathview(gene.data = t_RPKM_3_pv_mat,
                          species = "ko",
                          pathway.id="00010",
                          out.suffix = "number of different pathways class", 
                          gene.idtype = "KEGG",
                          limit = 2,
                          kegg.native = T,
                          same.layer= T)
t_RPKM_3_gylc.pv <- pathview(gene.data = t_RPKM_3_pv_mat,
                             species = "ko",
                             pathway.id="00010",
                             out.suffix = "number of different pathways class", 
                             gene.idtype = "KEGG",
                             limit = 2,
                             kegg.native = T,
                             same.layer= T)

### Methane - 00680

t_RPKM_3_meth.pv <- pathview(gene.data = t_RPKM_3_pv_mat,
                          species = "ko",
                          pathway.id="01200",
                          out.suffix = "number of different pathways class", 
                          gene.idtype = "KEGG",
                          limit = 2,
                          kegg.native = T,
                          same.layer= T)
t_RPKM_3_meth.pv <- pathview(gene.data = t_RPKM_3_pv_mat,
                             species = "ko",
                             pathway.id="01200",
                             out.suffix = "number of different pathways class", 
                             gene.idtype = "KEGG",
                             limit = 2,
                             kegg.native = T,
                             same.layer= T)
### Propanoate - 00640

t_RPKM_3_prop.pv <- pathview(gene.data = t_RPKM_3_pv_mat,
                          species = "ko",
                          pathway.id="01200",
                          out.suffix = "number of different pathways class", 
                          gene.idtype = "KEGG",
                          limit = 2,
                          kegg.native = T,
                          same.layer= T)

### Inositol phosphate - 00562

t_RPKM_3_phosp.pv <- pathview(gene.data = t_RPKM_3_pv_mat,
                             species = "ko",
                             pathway.id="00562",
                             out.suffix = "number of different pathways class", 
                             gene.idtype = "KEGG",
                             limit = 2,
                             kegg.native = T,
                             same.layer= T)

### carbon (doesn't seem to work in KEGG.native) - 01200

t_RPKM_3_c.pv <- pathview(gene.data = t_RPKM_3_pv_mat,
                          species = "ko",
                          pathway.id="01200",
                          out.suffix = "ko.data_class", 
                          gene.idtype = "KEGG",
                          limit = 2,
                          kegg.native = F,
                          same.layer= T)

t_RPKM_3_c.pv <- pathview(gene.data = t_RPKM_3_pv_mat,
                          species = "ko",
                          pathway.id="01200",
                          out.suffix = "ko.data_class", 
                          gene.idtype = "KEGG",
                          limit = 2,
                          kegg.native = F,
                          same.layer= F)


t_RPKM_3_class_df_c <- t_RPKM_3_c.pv$plot.data.gene
