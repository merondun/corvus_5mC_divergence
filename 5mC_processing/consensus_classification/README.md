# Overlap among experiments

Require an FDR of 1%, a primary effect divergence of 25%, non-significance for other effects. Conserved sites exhibit less than a 10% difference in methylation range. 

After the first script, run this:

And now intersect with our REGION, DXY, FST, PI, data: 

```bash
rgndir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Regions
gc=$rgndir/GCF_000738735.5_ASM73873v5_genomic.GC-Content.500bp.bed
rgn=$rgndir/GCF_000738735.5_ASM73873v5_genomic.GenomicFeaturesCDS.bed
popdat=$rgndir/GCF_000738735.5_ASM73873v5_genomic.PopGen2023JAN.bed

#create sites file 
sed '1d' Full.DSS.BP-10x.Classified-20250107.txt | awk '{OFS="\t"}{print $2, $3, $3, $1}' > Sites.bed

#create master overlapped genetic file 
bedtools intersect -wao -a Sites.bed -b $rgn | bedtools intersect -wao -a - -b $gc | bedtools intersect -header -wao -a - -b $popdat | awk '{OFS="\t"}{print $1, $16, $17, $4, $2, $8, $13, $18, $19, $20, $21, $22}' > All-Overlapped.bed
```

and re-merge with the full dataset

```r
library(tidyverse)
options(scipen=999)
init <- read.table('Full.DSS.BP-10x.Classified-20250107.txt',header=T)
anno <- read.table('All-Overlapped.bed',header=F)
names(anno) <- c('chr','GenStart','GenEnd','site','start', 'Region','GC','HapD','FuLiD','TajimaD','Dxy','FST')
all <- full_join(init,anno)
all$Region <- gsub('\\.','Missing',all$Region)
nrow(all)
write.table(all,'Full.DSS.BP-10x.Classified-Annotated_20250107.txt',quote=F,sep='\t',row.names=F)
```





And merge with our REGION, DXY, FST, PI, data: 

```bash
rgndir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/crow_scratch/2021-07_HybridZone/Nov/Regions
gc=$rgndir/GCF_000738735.5_ASM73873v5_genomic.GC-Content.500bp.bed
rgn=$rgndir/GCF_000738735.5_ASM73873v5_genomic.GenomicFeaturesCDS.bed
popdat=$rgndir/GCF_000738735.5_ASM73873v5_genomic.PopGen2023JAN.bed

#create sites file 
sed '1d' Full.DSS.BP-10x.Classified-20250107.txt | awk '{OFS="\t"}{print $2, $3, $3, $1}' > Sites_Year.bed

#create master overlapped genetic file 
bedtools intersect -wao -a Sites_Year.bed -b $rgn | bedtools intersect -wao -a - -b $gc | bedtools intersect -header -wao -a - -b $popdat | awk '{OFS="\t"}{print $1, $16, $17, $4, $2, $8, $13, $18, $19, $20, $21, $22}' > All-Overlapped-Year.bed
```

and re-merge with the full dataset

```r
library(tidyverse)
options(scipen=999)
init <- read.table('Full.DSS.BP-10x.Classified-23AUG30_YEAR.txt',header=T)
anno <- read.table('All-Overlapped-Year.bed',header=F)
names(anno) <- c('chr','GenStart','GenEnd','site','start', 'Region','GC','HapD','FuLiD','TajimaD','Dxy','FST')
all <- full_join(init,anno)
all$Region <- gsub('\\.','Missing',all$Region)
nrow(all)
all = all %>% dplyr::rename(Class = DMP_Strict)
write.table(all,'Full.DSS.BP-10x.Classified-Annotated-YEAR_23AUG30.txt',quote=F,sep='\t',row.names=F)
```


