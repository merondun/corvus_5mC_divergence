.libPaths('~/mambaforge/envs/r/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_202504')
library(tidyverse)
library(UpSetR)
options(scipen=999)
all = read_tsv('Full-HZ-CG.DSS.BP-10x.Classified-Annotated-SNPless_20250410.txt')

ap = all %>% select(site,DMP.HZ,DMP.CG)
ap = ap %>% mutate(across(starts_with("DMP"), ~ ifelse(. == "MISSING", "Unknown", .)))
names(ap) = c('site','HybZon','ComGar')

#only include sites which at least 1 taxon occurrence 
taxon_data <- ap %>%
  filter(HybZon == "Taxon" | ComGar == "Taxon")

# binary presence/absence matrix 
binary_matrix <- taxon_data %>%
  mutate(HybZon = ifelse(HybZon == "Taxon", 1, 0),
         ComGar = ifelse(ComGar == "Taxon", 1, 0)) %>%
  select(HybZon, ComGar)

# sum up the binary flags to get counts for each combination
aggregated_counts <- binary_matrix %>%
  group_by(HybZon, ComGar) %>%
  summarise(count = n(), .groups = 'drop')

# 'exression' data format for upsetR 
upset_data <- c(
  HybZon = aggregated_counts$count[aggregated_counts$HybZon == 1 & aggregated_counts$ComGar == 0],
  ComGar = aggregated_counts$count[aggregated_counts$HybZon == 0 & aggregated_counts$ComGar == 1],
  "HybZon&ComGar" = aggregated_counts$count[aggregated_counts$HybZon == 1 & aggregated_counts$ComGar == 1]
)

up = upset(fromExpression(upset_data), order.by = "freq")
up

pdf('~/symlinks/hzone/Methylation_Analyses/figures/20250416_Upset_Taxon_Plot.pdf',height=4,width=6)
up
dev.off()

