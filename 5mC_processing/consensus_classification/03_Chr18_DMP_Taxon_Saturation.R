setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_202504')
.libPaths('~/mambaforge/envs/r/lib/R/library')
### Are DMPs enriched on chr18?
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

hz <- read.table('HZ_DSS.10x.txt',header=TRUE)
hz_cands <- hz %>% filter(fdrs_Chr18_Hybrid_Index.HZ < 0.01) %>% select(chr,start) %>% mutate(Experiment = 'HybZon')

cg <- read.table('CG_DSS.10x.txt',header=TRUE)
cg_cands <- cg %>% filter(fdrs_SpeciesS.CG < 0.01) %>% select(chr,start) %>% mutate(Experiment = 'ComGar')

taxcands <- rbind(hz_cands, cg_cands) %>% filter(!grepl('scaf',chr)) %>% as_tibble
write_zero <- taxcands %>% mutate(end = start, start = start -1 ) %>% select(chr,start,end,Experiment)
write.table(write_zero,file='CG-HZ_DMPsTaxon.bed',quote=F,sep='\t',row.names=F,col.names = F)

chr_ord <- taxcands %>% select(chr) %>% unique %>%  
  mutate(chrs = gsub('[chr.A]', '', chr),   # Remove 'chr', '.', and 'A'
         chrs = as.numeric(gsub('Z', '99', chrs))) %>% arrange(chrs)
taxcands$chr <- factor(taxcands$chr,levels=chr_ord$chr)

# Unadjusted 
count_plot <- taxcands %>% 
  ggplot(aes(y=chr,fill=Experiment))+
  geom_bar()+
  xlab('Count of Taxon DMPs')+ylab('')+
  scale_fill_manual(values=c('#33638d','#b8de29'))+
  theme_bw(base_size=8)
count_plot

taxcounts <- taxcands %>% 
  count(chr,Experiment)


## correct for counts by the number of DMPs assayed by chromosome
cg_counts <- cg %>% count(chr) %>% filter(!grepl('MT|scaf',chr)) %>% arrange(desc(n)) %>% mutate(Experiment = 'ComGar')
hz_counts <- hz %>% count(chr) %>% filter(!grepl('MT|scaf',chr)) %>% arrange(desc(n)) %>% mutate(Experiment = 'HybZon')
xp_counts <- rbind(cg_counts,hz_counts)
all_chrs <- unique(xp_counts$chr)

# init df
binomial_results <- data.frame()

# loop xp & chr
for (exp in unique(taxcounts$Experiment)) {
  exp_data <- taxcounts %>% 
    filter(Experiment == exp) %>% 
    complete(chr = all_chrs, Experiment, fill = list(n = 0)) %>% left_join(.,xp_counts %>% dplyr::rename(chr_total = n)) %>% 
    mutate(gw_total = sum(chr_total),
           expected_prop = chr_total / gw_total)
  exp_total_dmps <- sum(exp_data$n)
  
  for (chr in unique(exp_data$chr)) {
    chr_dmps <- exp_data %>% filter(chr == !!chr) %>% pull(n)
    chr_expected_prop <- exp_data %>% filter(chr == !!chr) %>% pull(expected_prop)
    
    # bin test
    bin_test <- binom.test(chr_dmps, exp_total_dmps, chr_expected_prop)
    
    result_row <- data.frame(
      Experiment = exp,
      Chromosome = chr,
      chr_DMPs = chr_dmps,
      Total_DMPs = exp_total_dmps,
      Expected_proportion = chr_expected_prop,
      p_value = bin_test$p.value,
      lower_CI = bin_test$conf.int[1],
      upper_CI = bin_test$conf.int[2],
      observed_proportion = bin_test$estimate
    )
    
    binomial_results <- bind_rows(binomial_results, result_row)
  }
}

binomial_results <- binomial_results %>%
  mutate(significant = (lower_CI > Expected_proportion) )
binomial_results$Chromosome <- factor(binomial_results$Chromosome,levels=chr_ord$chr)
write.table(binomial_results,file='~/symlinks/hzone/Methylation_Analyses/figures/20250418_Binomial_Results.txt',quote=F,sep='\t',row.names=F)

binom_plot <- binomial_results %>%
  ggplot(aes(y = Chromosome, color = Experiment)) +
  
  # Observed proportion and 95% CI
  geom_point(aes(x = observed_proportion), size = 2) +
  geom_errorbarh(aes(xmin = lower_CI, xmax = upper_CI), height = 0.2) +
  
  # Expected proportion
  geom_point(aes(x = Expected_proportion), shape = 4, size = 2, stroke = 0.5, color = "black") +
  
  # asterisk for significant deviations
  geom_text(data = . %>% filter(significant),
            aes(x = max(upper_CI, Expected_proportion) + 0.005, label = "*"),  # nudge slightly right
            color = "black", size = 3) +
  
  facet_wrap(~ Experiment, scales = "free_y") +
  scale_color_manual(values=c('#33638d','#b8de29'))+
  labs(
    x = "Proportion of DMPs",
    y = "Chromosome"
  ) +
  theme_bw(base_size = 8)
binom_plot

ggsave(ggarrange(count_plot,binom_plot,common.legend = TRUE),file='~/symlinks/hzone/Methylation_Analyses/figures/20250418_DMP_Counts.pdf',height=5,width=6)


