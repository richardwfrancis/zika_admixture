library(tidyverse)
setwd("/Users/rfrancis_adm/Documents/TKI/GeneticsHealth/Timo/zika/scientific_data_paper/admixture")

# set colours
require("RColorBrewer")

# 20130606_g1k.ped comes from 1KG and is in /data/workspace/richard/zika/admixture/1KG
# zika.ped.txt was made by me based on Copy of Clinical_data_exome_subset.xlsx but using the IDs from the vcf
# (cat 20130606_g1k.ped | cut -f 1,2,3,4,5,6,7; cat zika.ped.txt) > 20130606_g1k.zika.ped
# perl -pi -e 's!\r!!' 20130606_g1k.zika.ped
# PED file for 1KG and zika
# filter the unwanted samples
PED <- read.table("20130606_g1k.zika.ped", header = TRUE, skip = 0, sep = "\t", stringsAsFactors = FALSE) %>%
  filter(! Individual.ID %in% c("JN-132", "JN-150", "JN-233"))

# PCA data for 1KG and zika
eigenvec <- read.table("plink.eigenvec.zika", header = FALSE, skip=0, sep = " ", stringsAsFactors = FALSE) %>%
  left_join(PED, by = c("V2"="Individual.ID")) %>%
  column_to_rownames("V2") %>%
  select(-"V1")
colnames(eigenvec) <- c(paste("Principal Component ", c(1:20), sep = ""),
                         "Family.ID", "Paternal.ID", "Maternal.ID", "Gender",
                         "Phenotype", "Population")

eigenvec <- eigenvec %>%
  rownames_to_column("rowname") %>%
  mutate(region=case_when(
    Population %in% c("ACB","ASW","ESN","GWD","LWK","MSL","YRI") ~ "African",
    Population %in% c("CLM","MXL","PEL","PUR") ~ "Hispanic",
    Population %in% c("CDX","CHB","CHS","JPT","KHV") ~ "East-Asian",
    Population %in% c("CEU","FIN","GBR","IBS","TSI") ~ "Caucasian",
    Population %in% c("BEB","GIH","ITU","PJL","STU") ~ "South Asian",
    Population %in% c("ZIKA") ~ "Zika",
    TRUE ~ "Other"
    )
    ) %>%
  filter(! Population == "Other") %>%
  column_to_rownames("rowname")

# generate PCA bi-plots
project.pca <- eigenvec[,1:20]
summary(project.pca)

ggplot(data = project.pca, aes(x=`Principal Component 1`, y = `Principal Component 2`, colour=eigenvec$region)) +
  geom_point() +
  theme_bw() +
  theme(
    legend.position = c(0.5, 0.15),
    legend.direction="horizontal",
    legend.text=element_text(size=15,family="Arial"),
    axis.title.x=element_text(size=20,family="Arial"),
    axis.title.y=element_text(size=20,family="Arial"),
    axis.text.x=element_text(size=16,family="Arial"),
    axis.text.y=element_text(size=16,family="Arial")
    ) +
  scale_colour_discrete(name="") +
  labs(x="Principal Component 1", y = "Principal Component 2")
