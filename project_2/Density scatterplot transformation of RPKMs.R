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


# Mapping the phylum 

phylum_checkM <- read_tsv(file = "bin_complete_contam_tax_phylum.txt", col_names = T)
phylum_checkM_1 <- phylum_checkM %>% dplyr::rename(Phylum = `taxonomy (p)`) %>%
  separate(bin_name, into = c("Name", "Bin"), sep = "m.")
phylum_checkM_1$Bin <- as.integer(phylum_checkM_1$Bin)

# Map RPKM prokka to KO based on MAG_orf

MAGs_ORFs_KO_2 <- read_tsv(file = "SaanichInlet_200m_all_MAGs_ORFs_ko.cleaned.txt", col_names = F)
colnames(MAGs_ORFs_KO_2) <- c("MAG_orf", "ko")

RPKM_prokka_KO <- left_join(MAGs_ORFs_KO_2, RPKM_prokka, by = "MAG_orf")
RPKM_prokka_KO_inner_join <- inner_join(MAGs_ORFs_KO_2, RPKM_prokka, by = "MAG_orf")
reverse_join_RPKM_prokka_KO <- inner_join(RPKM_prokka, MAGs_ORFs_KO_2, by = "MAG_orf")

# Mapping RPKM prokka onto phylum

RPKM_prokka_KO_phylum <- left_join(RPKM_prokka_KO, phylum_checkM_1,
                                   by = "Bin")
#
#
#
#
#
#
##### Rationale for log-transform of RPKM values

comparing_RPKM_log_RPKM <- RPKM_prokka_KO_phylum %>%
  dplyr::select(-Domain, -Bin, -Name, -Phylum, -completeness, -contamination) %>%
  separate(MAG_orf, into = c("mag", "orf_ID")) %>%
  gather(key = "Cruise", value = "RPKM", -mag, -orf_ID, -ko) %>%
  group_by(mag, ko) %>%
  summarize(total = sum(RPKM),
            log_total = log10(total)) %>%
  filter(log_total != -Inf) %>%
  melt(id = c("mag", "ko")) %>%
  separate(ko, into = c("K", "id"), sep = "K") %>%
  dplyr::select(-K)

comparing_RPKM_log_RPKM$id <- as.numeric(comparing_RPKM_log_RPKM$id)

ggplot(comparing_RPKM_log_RPKM, aes(x = id, y = value)) +
  geom_point(alpha = 0.15, colour = "seagreen", size = 2, pch = 19) +
  geom_point(alpha = 0.15, colour = "black", size = 2, pch = 1) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        panel.background = element_rect(fill = "white")) +
  facet_wrap(~ variable, scales = "free_y")

myColor <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
myColor_scale_fill <- scale_fill_gradientn(colours = myColor)

# There's quite a bit of overplotting in the previous figures, and alpha adjustment doesn't do a whole lot.
#so incorporating a density element to the
# scatterplot may provide a better sense of the data 

library(cowplot)

#Density 

normal_RPKM_figure <- comparing_RPKM_log_RPKM %>%
  filter(variable == "total") %>%
  ggplot(aes(x = id, y = value, fill = ..density..)) +
  stat_binhex(bins = 64, alpha = 0.95) +
  theme_bw() +
  scale_fill_viridis() +
  theme(axis.text.x = element_blank(),
        axis.text = element_text(size = 14),
        text = element_text(size = 14), 
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = seq(0, 35000, 5000)) +
  labs(x = "Open reading frame identifier (KO)", y = "RPKM",
       subtitle = "Density scatterplot of RPKM values")

## Normality

### non-transformed

normal_RPKM_values <- comparing_RPKM_log_RPKM %>%
  filter(variable == "total")

RPKM_slope <- diff(quantile(normal_RPKM_values$value, c(0.25, 0.75), na.rm = T)) / 
  diff(qnorm(c(0.25, 0.75)))

RPKM_intercept <- quantile(normal_RPKM_values$value, 0.25, na.rm = T) - 
  RPKM_slope * qnorm(0.25)

normal_RPKM_plot_1 <- ggplot(normal_RPKM_values, aes(sample = value)) +
  stat_qq(pch = 1, colour = "black",  fill = "black", alpha = 0.35, size = 4) +
  stat_qq(pch = 19, colour = "black",  fill = "black", alpha = 0.15, size = 4) +
  geom_abline(aes(slope = RPKM_slope, intercept = RPKM_intercept),
              size = 1, col = "red", alpha = 0.5, lty = 1) +
  labs(y = "RPKM", x = "Theoretical distribution")  +
  scale_x_continuous(breaks = seq(-4, 4, 1)) +
  scale_y_continuous(breaks = seq(0, 30000, 5000)) +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        text = element_text(size = 14), 
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank()) 

normal_RPKM_plot_2 <- axis_canvas(normal_RPKM_plot_1, axis = "y") +
  geom_vridgeline(data = normal_RPKM_values, aes(y = value, x = 0, width = ..density..),
                  stat = "ydensity", alpha = 0.7, size = 0.75, trim = F) 

normal_RPKM_dens_qq_plot <- insert_yaxis_grob(normal_RPKM_plot_1, normal_RPKM_plot_2, grid::unit(0.2, "null"),
                                              position = "right")

normal_RPKM_qq_figure <- ggdraw(normal_RPKM_dens_qq_plot)

### histogram

ggplot(normal_RPKM_values, aes())

## log-transformed

log_RPKM_figure <- comparing_RPKM_log_RPKM %>%
  filter(variable == "log_total") %>%
  ggplot(aes(x = id, y = value, fill = ..density..)) +
  stat_binhex(bins = 64, alpha = 0.95) +
  scale_fill_viridis() +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text = element_text(size = 14),
        text = element_text(size = 14), 
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = seq(-2, 4, 1)) +
  labs(x = "Open reading frame identifier (KO)", y = "Log RPKM",
       subtitle = "Density scatterplot of Log RPKM values") +
  annotation_logticks(sides = "l")

### DEnsity plot log transformed

log_RPKM_values <- comparing_RPKM_log_RPKM %>%
  filter(variable == "log_total")

log_RPKM_slope <- diff(quantile(log_RPKM_values$value, c(0.25, 0.75), na.rm = T)) / 
  diff(qnorm(c(0.25, 0.75)))

log_RPKM_intercept <- quantile(log_RPKM_values$value, 0.25, na.rm = T) - 
  log_RPKM_slope * qnorm(0.25)

log_RPKM_plot_1 <- ggplot(log_RPKM_values, aes(sample = value)) +
  stat_qq(pch = 1, colour = "black",  fill = "black", alpha = 0.35, size = 4) +
  stat_qq(pch = 19, colour = "black",  fill = "black", alpha = 0.15, size = 4) +
  geom_abline(aes(slope = log_RPKM_slope, intercept = log_RPKM_intercept),
              size = 1, col = "red", alpha = 0.5, lty = 1) +
  labs(y = "Log RPKM", x = "Theoretical distribution")  +
  scale_x_continuous(breaks = seq(-4, 4, 1)) +
  scale_y_continuous(breaks = seq(-1, 4, 1)) +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        text = element_text(size = 14), 
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank()) +
  annotation_logticks(sides = "l")

log_RPKM_plot_2 <- axis_canvas(log_RPKM_plot_1, axis = "y") +
  geom_vridgeline(data = log_RPKM_values, aes(y = value, x = 0, width = ..density..),
                  stat = "ydensity", alpha = 0.7, size = 0.75, trim = F, fill = "darkgrey") 

log_RPKM_dens_qq_plot <- insert_yaxis_grob(log_RPKM_plot_1, log_RPKM_plot_2, grid::unit(0.2, "null"),
                                              position = "right")

log_RPKM_qq_figure <- ggdraw(log_RPKM_dens_qq_plot)


plot_grid(normal_RPKM_figure, log_RPKM_figure, 
          normal_RPKM_qq_figure, log_RPKM_qq_figure, labels=c("A", "B", "C", "D"), align="hv", axis="tb", ncol = 2)
