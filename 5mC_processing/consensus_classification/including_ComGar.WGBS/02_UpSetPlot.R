.libPaths('~/mambaforge/envs/r/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_Sept')
library(tidyverse)
library(UpSetR)
options(scipen=999)
all = read_tsv('Full.DSS.BP-10x.Classified-Annotated_20250107.txt')

ap = all %>% select(site,DMP.wgbs,DMP.hz,DMP.cg)
ap = ap %>% mutate(across(starts_with("DMP"), ~ ifelse(. == "MISSING", "Unknown", .)))
names(ap) = c('site','ComGar.WGBS','HybZon.RRBS','ComGar.RRBS')

#only include sites which at least 1 taxon occurrence 
taxon_data <- ap %>%
  filter(ComGar.WGBS == "Taxon" | HybZon.RRBS == "Taxon" | ComGar.RRBS == "Taxon")

# binary presence/absence matrix 
binary_matrix <- taxon_data %>%
  mutate(ComGar_WGBS = ifelse(ComGar.WGBS == "Taxon", 1, 0),
         HybZon_RRBS = ifelse(HybZon.RRBS == "Taxon", 1, 0),
         ComGar_RRBS = ifelse(ComGar.RRBS == "Taxon", 1, 0)) %>%
  select(ComGar_WGBS, HybZon_RRBS, ComGar_RRBS)

# sum up the binary flags to get counts for each combination
aggregated_counts <- binary_matrix %>%
  group_by(ComGar_WGBS, HybZon_RRBS, ComGar_RRBS) %>%
  summarise(count = n(), .groups = 'drop')

# 'exression' data format for upsetR 
upset_data <- c(
  ComGar_WGBS = aggregated_counts$count[aggregated_counts$ComGar_WGBS == 1 & aggregated_counts$HybZon_RRBS == 0 & aggregated_counts$ComGar_RRBS == 0],
  HybZon_RRBS = aggregated_counts$count[aggregated_counts$ComGar_WGBS == 0 & aggregated_counts$HybZon_RRBS == 1 & aggregated_counts$ComGar_RRBS == 0],
  ComGar_RRBS = aggregated_counts$count[aggregated_counts$ComGar_WGBS == 0 & aggregated_counts$HybZon_RRBS == 0 & aggregated_counts$ComGar_RRBS == 1],
  "ComGar_WGBS&HybZon_RRBS" = aggregated_counts$count[aggregated_counts$ComGar_WGBS == 1 & aggregated_counts$HybZon_RRBS == 1 & aggregated_counts$ComGar_RRBS == 0],
  "ComGar_WGBS&ComGar_RRBS" = aggregated_counts$count[aggregated_counts$ComGar_WGBS == 1 & aggregated_counts$HybZon_RRBS == 0 & aggregated_counts$ComGar_RRBS == 1],
  "HybZon_RRBS&ComGar_RRBS" = aggregated_counts$count[aggregated_counts$ComGar_WGBS == 0 & aggregated_counts$HybZon_RRBS == 1 & aggregated_counts$ComGar_RRBS == 1],
  "All" = aggregated_counts$count[aggregated_counts$ComGar_WGBS == 1 & aggregated_counts$HybZon_RRBS == 1 & aggregated_counts$ComGar_RRBS == 1]
)

up = upset(fromExpression(upset_data), order.by = "freq")
up

pdf('20250107_Upset_Taxon_Plot.pdf',height=4,width=7)
up
dev.off()

