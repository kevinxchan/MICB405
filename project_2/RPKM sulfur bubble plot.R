# Bubble plot sulphur 

sulfur_bubbs <- sulfur_data_frame_coerced %>%
  dplyr::select(-X__1, -type, -x, -y, -width, -height) %>%
  dplyr::select(kegg.names:`<NA>`) %>%
  gather(key = "Class_identity", value = "Log_sum_RPKMs",
         -kegg.names, -all.mapped, -labels)

sulfur_bubbs$Log_sum_RPKMs <- as.numeric(sulfur_bubbs$Log_sum_RPKMs)

sulfur_bubbs_1 <- sulfur_bubbs %>%
  separate(kegg.names, into = c("K", "kegg_ID"), sep = "K") %>%
  dplyr::select(-K)

sulfur_bubbs_1$kegg_ID <- as.numeric(sulfur_bubbs_1$kegg_ID)
  
  
ggplot(sulfur_bubbs_1, aes(y = kegg_ID, x = Class_identity)) +
  geom_point(aes(size = Log_sum_RPKMs)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = Log_sum_RPKMs), size = 1)


#RPKMs for the bubble analysis
sulfur_bubbs_normal <- Sulfur_pathview_gene_df_coerced %>%
  dplyr::select(-X__1, -type, -x, -y, -width, -height) %>%
  dplyr::select(kegg.names:`<NA>`) %>%
  gather(key = "Class_identity", value = "sum_RPKMs",
         -kegg.names, -all.mapped, -labels)

sulfur_bubbs_normal_1 <- sulfur_bubbs_normal


sulfur_bubbs_normal_1$kegg_ID <- as.numeric(sulfur_bubbs_normal_1$kegg_ID)
sulfur_bubbs_normal_1$sum_RPKMs <- as.numeric(sulfur_bubbs_normal_1$sum_RPKMs)

sulfur_bubbs_normal_2 <- sulfur_bubbs_normal_1 %>%
  separate(Class_identity, into = c("Class", "Domain"), by = " ") # The warnings are for getting rid of the parentheses around domain name, so no worries there




ggplot(sulfur_bubbs_normal_2, aes(y = kegg.names, x = Class, colour = Domain)) +
  geom_point(aes(size = log10(sum_RPKMs)), alpha = 0.35, fill = "black") +
  scale_size(range = c(1, 10), name = "Log sum\nof RPKMs") +
  coord_flip() +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Class", y = "Gene (KO identifier)", subtitle = "Sulfur metabolic pathway") +
  geom_text(aes(label = round(log10(sum_RPKMs), digits = 2)), colour = "black", size = 2)

#
#
#
#
# Since the KO identifiers are not human readable, they arent meaningful.
## So, using inputing a list of the KO identifiers into https://www.genome.jp/kegg/tool/map_brite1.html
### Gave us a list of gene names, which we will now join and plot as a function of taxonomic rank (or coord_flip)

library(readxl)
KO_gene_names <- read_excel("Sulfur pathway KO gene identifiers.xlsx", 
                      col_names = FALSE)

colnames(KO_gene_names) <- c("kegg.names", "gene.name", "higher_order_grouping", "remove_this")

KO_gene_names_selected <- KO_gene_names %>%
  dplyr::select(-remove_this)

sulfur_bubbs_2 <- Sulfur_pathview_gene_df_coerced %>%
  dplyr::select(-X__1, -type, -x, -y, -width, -height) %>%
  dplyr::select(kegg.names:`<NA>`) %>%
  dplyr::rename(`Class level unidentified` = `<NA>`)

sulfur_bubbs_joined <- sulfur_bubbs_2 %>%
  left_join(KO_gene_names_selected, 
                                 by = "kegg.names") %>%
  gather(key = "Class_identity", value = "sum_RPKMs",
         -kegg.names, -all.mapped, -labels, -gene.name, -higher_order_grouping) %>%
  separate(Class_identity, into = c("Class", "Domain"), by = " ") %>%
  mutate(Class = gsub("Class", "Class level unidentified", Class),
         Domain = gsub("level", "Unidentified", Domain)) 

sulfur_bubbs_joined$sum_RPKMs <- as.numeric(sulfur_bubbs_joined$sum_RPKMs)
  
# NA composed of Bacteroidota and Nanoarchaeota that were not identified at the Class level"
sulfur_bubbs_joined_final <- sulfur_bubbs_joined %>%
  group_by(Class, gene.name, Domain) %>% 
  summarize(sum_Class_RPKMs = sum(sum_RPKMs),
            log_sum_Class_RPKMs = log10(sum_Class_RPKMs)) %>%
  mutate(log_sum_Class_RPKMs = gsub("-Inf", "NA", log_sum_Class_RPKMs)) 
  
sulfur_bubbs_joined_final$log_sum_Class_RPKMs <- as.numeric(sulfur_bubbs_joined_final$log_sum_Class_RPKMs)
  
ggplot(sulfur_bubbs_joined_final, aes(x = gene.name, y = Class, colour = Domain)) +
  geom_point(aes(size = log_sum_Class_RPKMs), alpha = 0.35, fill = "black") +
  scale_size(range = c(1, 10), name = "Log sum\nof RPKMs") +
  scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank()) +
  labs(y = "Class", x = "Gene name", subtitle = "Sulfur metabolic pathway") +
  geom_text(aes(label = round(log_sum_Class_RPKMs, digits = 2)), colour = "black", size = 2)

# Aya and Tina's manual pathway grouping

Adding_higher_order_grouping_ <- read_excel("Adding higher order grouping .xlsx")

sulfur_bubbs_hi_grouping_processing <- Sulfur_pathview_gene_df_coerced %>%
  dplyr::select(-X__1, -type, -x, -y, -width, -height) %>%
  dplyr::select(kegg.names:`<NA>`) %>%
  dplyr::rename(`Class level unidentified` = `<NA>`) %>%
  gather(key = "Class", value = "sum_RPKMs", 
         -kegg.names, -labels, -all.mapped) 
  #filter(sum_RPKMs != "NA") %>%
  #filter(sum_RPKMs != "0") %>%
  #mutate(sum_RPKMs = log10(as.numeric(sum_RPKMs))) %>%
  #spread(key = "Class", value = "sum_RPKMs")

sulfur_bubbs_hi_grouping_processing$sum_RPKMs <- as.numeric(sulfur_bubbs_hi_grouping_processing$sum_RPKMs)

hi_grouping_processed <- Adding_higher_order_grouping_ %>%
  dplyr::rename(kegg.names = `KO Number`,
                gene.name = Gene) 

sulfur_bubbs_hi_grouping <- left_join(sulfur_bubbs_hi_grouping_processing, hi_grouping_processed, 
                                      by = c("kegg.names", "Class")) %>%
  separate(Class, into = c("Class", "Domain"), by = " ") 

sulfur_bubbs_hi_grouping_1 <- sulfur_bubbs_hi_grouping

######
#
#
###
#
#
######

KO_gene_names <- read_excel("Sulfur pathway KO gene identifiers.xlsx", 
                            col_names = FALSE)

colnames(KO_gene_names) <- c("kegg.names", "gene.name", "higher_order_grouping", "remove_this")

KO_gene_names_selected <- KO_gene_names %>%
  dplyr::select(-remove_this)

hi_grouping_processed <- Adding_higher_order_grouping_ %>%
  dplyr::rename(all.mapped = `KO Number`,
                gene.name_1 = Gene) %>%
  dplyr::select(-Class, -`RPKM value`)

KO_gene_hi_order <- right_join(hi_grouping_processed, KO_gene_names_selected,
                               by = "kegg.names") %>%
  dplyr::select(-gene.name_1, -higher_order_grouping)
#
sulfur_bubbs_3 <- Sulfur_pathview_gene_df_coerced %>%
  dplyr::select(-X__1, -type, -x, -y, -width, -height) %>%
  dplyr::select(kegg.names:`<NA>`) %>%
  dplyr::rename(`Class level unidentified` = `<NA>`) 

sulfur_bubbs_hi_joined <- sulfur_bubbs_3 %>%
  left_join(hi_grouping_processed, 
            by = "all.mapped") %>%
  gather(key = "Class_identity", value = "sum_RPKMs",
         -kegg.names, -all.mapped, -labels) %>%
  separate(Class_identity, into = c("Class", "Domain"), by = " ") %>%
  mutate(Class = gsub("Class", "Class level unidentified", Class),
         Domain = gsub("level", "Unidentified", Domain)) 

sulfur_bubbs_hi_joined$sum_RPKMs <- as.numeric(sulfur_bubbs_hi_joined$sum_RPKMs)












