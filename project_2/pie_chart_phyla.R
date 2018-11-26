library(ggplot2)
library(tidyverse)
library(ggpubr)
library(cowplot)

raw_dat_bac <- read_tsv("gtdbtk.bac120.classification_pplacer.tsv")
raw_dat_ar <- read_tsv("gtdbtk.ar122.classification_pplacer.tsv")

phyla_dat_bac <- mutate(raw_dat_bac, phylum=sapply((str_split(raw_dat_bac[[2]], ";")), "[[", 2))
phyla_dat_bac <- mutate(phyla_dat_bac, phylum=gsub("p__", "", phylum))
phyla_dat_ar <- mutate(raw_dat_ar, phylum=sapply((str_split(raw_dat_ar[[2]], ";")), "[[", 2))
phyla_dat_ar <- mutate(phyla_dat_ar, phylum=gsub("p__", "", phylum))


df <- as.data.frame(table(select(phyla_dat_bac, phylum)))
df <- rename(df, Phylum = "Var1")
df.2 <- as.data.frame(table(select(phyla_dat_ar, phylum)))
df.2 <- rename(df.2, Phylum = "Var1")

labs <- vector(mode="character", length=length(df$Phylum))
labs.2 <- vector(mode="character", length=length(df.2$Phylum))

bac <- ggpie(df, "Freq", label = labs, fill = "Phylum")
arc <- ggpie(df.2, "Freq", label = labs.2, lab.font = c(2, "bold", "black"), fill = "Phylum")

ggsave("bac_phyla_pie_chart.png", bac, width=10, height=6)
ggsave("arc_phyla_pie_chart.png", arc, width=10, height=6)