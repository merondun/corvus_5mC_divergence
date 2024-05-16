## Consensus Classification

Require an FDR of 10%, a primary effect divergence of 25%, non-significance for other effects. Conserved sites exhibit less than a 10% difference in methylation range. 

```R
.libPaths('~/miniconda3/envs/randomforest/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/crow_scratch/2021-07_HybridZone/Nov/Classification_Sept')
library(dplyr)
library(tidyverse)
library(data.table)

#thresholds
pval <- 0.1
div <- 0.25
dd <- 0.1

#### CG
datsCG <- read.table('CG.DSS.BP-10x.txt',header=TRUE)
datsCG <- datsCG %>% dplyr::rename(fdrs_Species.cg = fdrs_SpeciesS.cg, fdrs_Sex.cg = fdrs_SexM.cg)

#with interacting effects, exclude variables significant for other variables
CG <- datsCG %>% mutate(DMP = ifelse(fdrs_Species.cg < pval & abs(divergence_Species.cg) > div &
                                       fdrs_Sex.cg > pval & fdrs_StageCHK.cg > pval & fdrs_StageYRL.cg > pval, 'Taxon',  #species
                                     ifelse(fdrs_Sex.cg < pval & abs(divergence_Sex.cg) > div &
                                              fdrs_Species.cg > pval & fdrs_StageCHK.cg > pval & fdrs_StageYRL.cg > pval, 'Sex', #sex
                                            ifelse((fdrs_StageCHK.cg < pval | fdrs_StageYRL.cg < pval) & (abs(divergence_StageCHK.cg) > div | abs(divergence_StageYRL.cg) > div) &
                                                     (fdrs_Sex.cg > pval & fdrs_Species.cg > pval),'Age',  #stage
                                                   ifelse(Divergence.cg < dd,  'Conserved','Unknown')))))

CG %>% dplyr::count(DMP)

#### WGBS
datsWGBS <- read.table('WGBS.DSS.BP-5x.txt',header=TRUE)

#with interacting effects, exclude variables significant for other variables
WGBS <- datsWGBS %>% mutate(DMP = ifelse(fdrs_Species.wgbs < pval & abs(divergence_Species.wgbs) > div &
                                           fdrs_Sex.wgbs > pval & fdrs_TissueM.wgbs > pval & fdrs_TissueL.wgbs > pval, 'Taxon',  #species
                                         ifelse(fdrs_Sex.wgbs < pval & abs(divergence_Sex.wgbs) > div &
                                                  fdrs_Species.wgbs > pval & fdrs_TissueM.wgbs > pval & fdrs_TissueL.wgbs > pval, 'Sex', #sex
                                                ifelse((fdrs_TissueM.wgbs < pval | fdrs_TissueL.wgbs < pval) & (abs(divergence_TissueM.wgbs) > div | abs(divergence_TissueL.wgbs) > div) &
                                                         (fdrs_Sex.wgbs > pval & fdrs_Species.wgbs > pval),'Tissue',  #tissue
                                                       ifelse(Divergence.wgbs < dd,  'Conserved','Unknown')))))


WGBS %>% dplyr::count(DMP)

#merged
dmr1 <- merge(WGBS,CG,by=c('site','chr','start'),
              suffixes=c('.wgbs','.cg'),all=TRUE)
#merge with HZ data..
hz <- read.table('HZ.DSS.BP-10x.txt',header=T)
hz = hz %>% mutate(DMP.hz = ifelse(fdrs_Chr18_Hybrid_Index.hz < pval & fdrs_Distance.hz > pval & fdrs_Temperature.hz > pval, 'Taxon',
                                   ifelse(fdrs_Chr18_Hybrid_Index.hz > pval & fdrs_Distance.hz < pval & fdrs_Temperature.hz > pval, 'Distance',
                                          ifelse(fdrs_Chr18_Hybrid_Index.hz > pval & fdrs_Distance.hz > pval & fdrs_Temperature.hz < pval, 'Temperature',
                                                 ifelse(Divergence.hz < dd, 'Conserved','Unknown')))))
hz %>% dplyr::count(DMP.hz)
dmr2 <- merge(dmr1,hz,by=c('site','chr','start'),all=TRUE)

#count NAs by experiment
nac = dmr2 %>% select(site,DMP.wgbs,DMP.cg,DMP.hz)
nac$NAs = rowSums(is.na(nac))
nac %>% dplyr::count(NAs)

#replace NA with missing
dmr3 = dmr2 %>% replace_na(list(DMP.wgbs = 'MISSING',DMP.cg ='MISSING',DMP.hz = 'MISSING'))
#now, classify based on consensus on the common garden: 
dmr4 = dmr3 %>% mutate(DMP_Strict = ifelse(DMP.wgbs == 'Tissue' & 
                                             ((DMP.cg == 'Conserved' & DMP.hz == 'MISSING') | 
                                                (DMP.hz == 'Conserved' & DMP.cg == 'MISSING') |
                                                (DMP.cg == 'Conserved' & DMP.hz == 'Conserved')), 'Tissue',
                                           ifelse(DMP.cg == 'Age' & 
                                                    ((DMP.wgbs == 'Conserved' & DMP.hz == 'MISSING') | 
                                                       (DMP.hz == 'Conserved' & DMP.wgbs == 'MISSING' ) |
                                                       (DMP.wgbs == 'Conserved' & DMP.hz == 'Conserved')), 'Age',
                                                  ifelse(DMP.hz == 'Conserved' & 
                                                           ((DMP.cg == 'Sex' & DMP.wgbs == 'Sex') | 
                                                              ( DMP.cg == 'MISSING' & DMP.wgbs == 'Sex') |
                                                              (DMP.cg == 'Sex' & DMP.wgbs == 'MISSING' )), 'Sex',
                                                         ifelse(DMP.hz == 'Taxon' & 
                                                                  ((DMP.cg == 'Taxon' & DMP.wgbs == 'Taxon') | 
                                                                     ( DMP.cg == 'MISSING' & DMP.wgbs == 'Taxon') |
                                                                     (DMP.cg == 'Taxon' & DMP.wgbs == 'MISSING' )), 'Taxon',
                                                                ifelse(DMP.hz == 'Temperature' & 
                                                                         ((DMP.cg == 'Conserved' & DMP.wgbs == 'MISSING' ) | 
                                                                            (DMP.wgbs == 'Conserved' & DMP.cg == 'MISSING') |
                                                                            (DMP.cg == 'Conserved' & DMP.wgbs == 'Conserved')), 'Temperature',
                                                                       ifelse(DMP.hz == 'Distance' & 
                                                                                ((DMP.cg == 'Conserved' & DMP.wgbs == 'MISSING' ) | 
                                                                                   (DMP.wgbs == 'Conserved' & DMP.cg == 'MISSING') |
                                                                                   (DMP.cg == 'Conserved' & DMP.wgbs == 'Conserved')), 'Distance',
                                                                              ifelse(DMP.cg == 'Conserved' & DMP.wgbs == 'Conserved', 'Conserved',
                                                                                     'Unclassified'))))))))
dmr4 = dmr4 %>% mutate(DMP_Single = ifelse(DMP.wgbs == 'Tissue' & 
                                             (DMP.cg == 'Conserved' | DMP.cg == 'Unknown' | DMP.cg == 'MISSING') &
                                             (DMP.hz == 'Conserved' | DMP.hz == 'Unknown' | DMP.hz == 'MISSING'), 'Tissue',
                                           ifelse(DMP.cg == 'Age' & 
                                                    (DMP.wgbs == 'Conserved' | DMP.wgbs == 'Unknown' | DMP.wgbs == 'MISSING') &
                                                    (DMP.hz == 'Conserved' | DMP.hz == 'Unknown' | DMP.hz == 'MISSING'), 'Age',
                                                  ifelse((DMP.hz == 'Conserved' | DMP.hz == 'Unknown' | DMP.hz == 'MISSING') & 
                                                           ((DMP.cg == 'Sex' & DMP.wgbs == 'Sex') | 
                                                              (( DMP.cg == 'MISSING' | DMP.cg == 'Unknown') & DMP.wgbs == 'Sex') |
                                                              (DMP.cg == 'Sex' & ( DMP.wgbs == 'MISSING' | DMP.wgbs == 'Unknown'))), 'Sex',
                                                         ifelse((DMP.hz == 'Taxon' | DMP.hz == 'Unknown' | DMP.hz == 'MISSING') & 
                                                                  ((DMP.cg == 'Taxon' & DMP.wgbs == 'Taxon') | 
                                                                     (( DMP.cg == 'MISSING' | DMP.cg == 'Unknown') & DMP.wgbs == 'Taxon') |
                                                                     (DMP.cg == 'Taxon' & ( DMP.wgbs == 'MISSING' | DMP.wgbs == 'Unknown'))), 'Taxon',
                                                                ifelse(DMP.hz == 'Temperature' & 
                                                                         (DMP.cg == 'Conserved' | DMP.cg == 'Unknown' | DMP.cg == 'MISSING') &
                                                                         (DMP.wgbs == 'Conserved' | DMP.wgbs == 'Unknown' | DMP.wgbs == 'MISSING'),'Temperature',
                                                                       ifelse(DMP.hz == 'Distance' & 
                                                                                (DMP.cg == 'Conserved' | DMP.cg == 'Unknown' | DMP.cg == 'MISSING') &
                                                                                (DMP.wgbs == 'Conserved' | DMP.wgbs == 'Unknown' | DMP.wgbs == 'MISSING'), 'Distance',
                                                                              ifelse((DMP.cg == 'Conserved' & (DMP.wgbs == 'Conserved' | DMP.wgbs == 'MISSING')) |
                                                                                       (DMP.wgbs == 'Conserved' & (DMP.cg == 'Conserved' | DMP.cg == 'MISSING')) , 'Conserved',
                                                                                     'Unclassified'))))))))

dmr4 %>% dplyr::count(DMP_Strict)
dmr4 %>% dplyr::count(DMP_Single)
dmr4 %>% group_by(DMP_Strict,DMP.wgbs,DMP.cg,DMP.hz) %>% count() %>% data.frame

dmr5 = dmr4 
dmr5$Mean5mC.cg = rowMeans(dmr5[grepl('^D_.*cg|^S_.*cg',names(dmr5))],na.rm=TRUE)
dmr5$Mean5mC.wgbs = rowMeans(dmr5[grepl('^D_.*wgbs|^S_.*wgbs',names(dmr5))],na.rm=TRUE)
dmr5$Mean5mC.hz = rowMeans(dmr5[grepl('^D_.*hz|^S_.*hz|^H_.*hz',names(dmr5))],na.rm=TRUE)                                                     

write.table(dmr5,'Full.DSS.BP-10x.Classified-23FEB09.txt',quote=F,sep='\t',row.names=FALSE)
```

### Merge Datasets

And now intersect with our REGION, DXY, FST, PI, data: 

```bash
rgndir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/crow_scratch/2021-07_HybridZone/Nov/Regions
gc=$rgndir/GCF_000738735.5_ASM73873v5_genomic.GC-Content.500bp.bed
rgn=$rgndir/GCF_000738735.5_ASM73873v5_genomic.GenomicFeaturesCDS.bed
popdat=$rgndir/GCF_000738735.5_ASM73873v5_genomic.PopGen2023JAN.bed

#create sites file 
sed '1d' Full.DSS.BP-10x.Classified-23FEB09.txt | awk '{OFS="\t"}{print $2, $3, $3, $1}' > Sites.bed

#create master overlapped genetic file 
bedtools intersect -wao -a Sites.bed -b $rgn | bedtools intersect -wao -a - -b $gc | bedtools intersect -header -wao -a - -b $popdat | awk '{OFS="\t"}{print $1, $16, $17, $4, $2, $8, $13, $18, $19, $20, $21, $22}' > All-Overlapped.bed
```

and re-merge with the full dataset

```r
library(tidyverse)
options(scipen=999)
init <- read.table('Full.DSS.BP-10x.Classified-23FEB09.txt',header=T)
anno <- read.table('All-Overlapped.bed',header=F)
names(anno) <- c('chr','GenStart','GenEnd','site','start', 'Region','GC','HapD','FuLiD','TajimaD','Dxy','FST')
all <- full_join(init,anno)
all$Region <- gsub('\\.','Missing',all$Region)
nrow(all)
all = all %>% dplyr::rename(Class = DMP_Strict)
write.table(all,'Full.DSS.BP-10x.Classified-Annotated_23FEB09.txt',quote=F,sep='\t',row.names=F)
```

## Lone Exp. Support: Focal Region

Taking a single-experimental approach: what's happening with hybrids in the focal region? Do hybrids that are more carrion crow-like exhibit patterns different than hooded-crow like, within the focal reigon compared to the autosomal background?

```R
.libPaths('~/miniconda3/envs/r/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/crow_scratch/2021-07_HybridZone/Nov/Classification_Sept')
library(dplyr)
library(maditr)
library(tidyverse)
library(data.table)

#thresholds
pval <- 0.1
div <- 0.25
dd <- 0.1

#### CG
datsCG <- read.table('CG.DSS.BP-10x.txt',header=TRUE)
datsCG <- datsCG %>% dplyr::rename(fdrs_Species.cg = fdrs_SpeciesS.cg, fdrs_Sex.cg = fdrs_SexM.cg)

#with interacting effects, exclude variables significant for other variables
CG <- datsCG %>% mutate(DMP = ifelse(fdrs_Species.cg < pval & abs(divergence_Species.cg) > div &
                                       fdrs_Sex.cg > pval & fdrs_StageCHK.cg > pval & fdrs_StageYRL.cg > pval, 'Taxon',  #species
                                     ifelse(fdrs_Sex.cg < pval & abs(divergence_Sex.cg) > div &
                                              fdrs_Species.cg > pval & fdrs_StageCHK.cg > pval & fdrs_StageYRL.cg > pval, 'Sex', #sex
                                            ifelse((fdrs_StageCHK.cg < pval | fdrs_StageYRL.cg < pval) & (abs(divergence_StageCHK.cg) > div | abs(divergence_StageYRL.cg) > div) &
                                                     (fdrs_Sex.cg > pval & fdrs_Species.cg > pval),'Age',  #stage
                                                   ifelse(Divergence.cg < dd,  'Conserved','Unknown')))))

CG %>% dplyr::count(DMP)

#### WGBS
datsWGBS <- read.table('WGBS.DSS.BP-5x.txt',header=TRUE)

#with interacting effects, exclude variables significant for other variables
WGBS <- datsWGBS %>% mutate(DMP = ifelse(fdrs_Species.wgbs < pval & abs(divergence_Species.wgbs) > div &
                                           fdrs_Sex.wgbs > pval & fdrs_TissueM.wgbs > pval & fdrs_TissueL.wgbs > pval, 'Taxon',  #species
                                         ifelse(fdrs_Sex.wgbs < pval & abs(divergence_Sex.wgbs) > div &
                                                  fdrs_Species.wgbs > pval & fdrs_TissueM.wgbs > pval & fdrs_TissueL.wgbs > pval, 'Sex', #sex
                                                ifelse((fdrs_TissueM.wgbs < pval | fdrs_TissueL.wgbs < pval) & (abs(divergence_TissueM.wgbs) > div | abs(divergence_TissueL.wgbs) > div) &
                                                         (fdrs_Sex.wgbs > pval & fdrs_Species.wgbs > pval),'Tissue',  #tissue
                                                       ifelse(Divergence.wgbs < dd,  'Conserved','Unknown')))))


WGBS %>% dplyr::count(DMP)

#merged
dmr1 <- merge(WGBS,CG,by=c('site','chr','start'),
              suffixes=c('.wgbs','.cg'),all=TRUE)
#merge with HZ data..
hz <- read.table('HZ.DSS.BP-10x.txt',header=T)
hz = hz %>% mutate(DMP.hz = ifelse(fdrs_Chr18_Hybrid_Index.hz < pval & fdrs_Distance.hz > pval & fdrs_Temperature.hz > pval, 'Taxon',
                                   ifelse(fdrs_Chr18_Hybrid_Index.hz > pval & fdrs_Distance.hz < pval & fdrs_Temperature.hz > pval, 'Distance',
                                          ifelse(fdrs_Chr18_Hybrid_Index.hz > pval & fdrs_Distance.hz > pval & fdrs_Temperature.hz < pval, 'Temperature',
                                                 ifelse(Divergence.hz < dd, 'Conserved','Unknown')))))
hz %>% dplyr::count(DMP.hz)
dmr2 <- merge(dmr1,hz,by=c('site','chr','start'),all=TRUE)

#count NAs by experiment
nac = dmr2 %>% select(site,DMP.wgbs,DMP.cg,DMP.hz)
nac$NAs = rowSums(is.na(nac))
nac %>% dplyr::count(NAs)

#replace NA with missing
dmr3 = dmr2 %>% replace_na(list(DMP.wgbs = 'MISSING',DMP.cg ='MISSING',DMP.hz = 'MISSING'))
#Without relying on significance without the HybZon experiment, just assign taxon / tissue and conserved sites 
dmr4 = dmr3 %>% mutate(TaxonLone = ifelse((DMP.cg == 'MISSING' & DMP.wgbs == 'Taxon') | 
                                            (DMP.wgbs == 'MISSING' & DMP.cg == 'Taxon') | 
                                            (DMP.wgbs == 'Taxon' & DMP.cg == 'Taxon'),'Taxon',
                                          ifelse((DMP.cg == 'MISSING' & DMP.wgbs == 'Tissue') | 
                                                   (DMP.wgbs == 'Tissue' & DMP.cg == 'Conserved'),'Tissue',
                                                 ifelse((DMP.cg == 'Age' & DMP.wgbs == 'MISSING') | 
                                                          (DMP.wgbs == 'Conserved' & DMP.cg == 'Age'),'Age',
                                                        ifelse((DMP.cg == 'MISSING' & DMP.wgbs == 'Sex') | 
                                                                 (DMP.wgbs == 'MISSING' & DMP.cg == 'Sex') | 
                                                                 (DMP.wgbs == 'Sex' & DMP.cg == 'Sex'),'Sex',
                                                              ifelse(DMP.wgbs == 'Conserved' & DMP.cg == 'Conserved','Conserved','Unknown'))))))

#Grab important columns
hz = dmr4 %>% dplyr::select(chr,start,TaxonLone,contains(c('.hz'))) %>% dplyr::select(-contains(c('fdrs','DMP','Div','Var','div')))
hz = hz %>% mutate(Focal_Region = ifelse(chr == 'chr18' & start > 8.07e6 & start < 10.07e6,'Focal','Undiverged'))
#only keep sites which are taxon/tissue/conserved 
hzm = hz %>% filter(TaxonLone != 'Unknown') 
hzm =  hzm[rowSums(is.na(hzm[grepl('^S_|^H_|^D_', names(hzm))])) <= 3, ] #remove sites with 3 or more individuals missing data 
#count sites
hzm %>% dplyr::count(TaxonLone,Focal_Region)
hzm2 = hzm %>% 
  pivot_longer(!c(chr,start,TaxonLone,Focal_Region))
hzm2 = hzm2 %>% mutate(name = gsub('_BL_CHK_M.hz','',name)) %>% na.omit
mdat = read.table('All-Metadata.txt',header=TRUE)
hzm3 = left_join(hzm2,mdat %>% select(Short_ID,Chr18_Hybrid_Index,Species) %>% dplyr::rename(name = Short_ID))
hzm3 = hzm3 %>% mutate(HI = ifelse(Species == 'C.corone','C.c. corone',
                                   ifelse(Species == 'C.cornix','C.c. cornix',
                                          ifelse(Chr18_Hybrid_Index < 0.25,'C.c. corone-like',
                                                 ifelse(Chr18_Hybrid_Index > 0.75,'C.c. cornix-like','Hybrid')))))

fullboots = NULL
options(dplyr.summarise.inform = FALSE)
set.seed(123)
B <- 999 # number of bootstrap replicates

for (grp in unique(hzm3$HI)) {
  for (vr in unique(hzm3$TaxonLone)) {
    cat('Working on Experiment & Variable: ',grp,' ',vr,'\n')
    bd = hzm3 %>% filter(HI == grp & TaxonLone == vr)
    n = bd %>% filter(Focal_Region == 'Focal') %>% nrow
    bdat=NULL
    ## Run the bootstrap procedure:
    for (b in 1:B) {
      cat('Iteration: ',b,'\n')
      d1 = bd %>% group_by(Focal_Region) %>% slice_sample(n=n,replace=TRUE) %>% summarize(mean = mean(value,na.rm=TRUE))
      bdat = rbind(bdat,d1 )
    }
    bdat$HI = grp
    bdat$Variable = vr
    names(bdat) = c('Focal_Region','Boots','HI','Variable')
    fullboots = rbind(fullboots,bdat)
  }
}

write.table(fullboots,file='HypoBoots_bg.txt',quote=F,sep='\t',row.names=F)
bdatplot = read.table('HypoBoots_bg.txt',header=TRUE,sep='\t')

#Plot 
bdatplot$HI = factor(bdatplot$HI,levels=c('C.c. corone','C.c. corone-like','Hybrid','C.c. cornix-like','C.c. cornix'))
#calculate CIs
cis = bdatplot %>% 
  group_by(HI,Variable,Focal_Region) %>% 
  summarize(LowerBound = quantile(Boots,0.025), #can change based on your significance threshold 
            UpperBound = quantile(Boots,0.975))
cis = cis %>% group_by(HI,Variable) %>% 
  mutate(Signif = ifelse(LowerBound[Focal_Region == 'Focal'] > UpperBound[Focal_Region == 'Undiverged'],'*',
                         ifelse(UpperBound[Focal_Region == 'Focal'] < LowerBound[Focal_Region == 'Undiverged'],'*','n.s.')))
bp = bdatplot %>% 
  ggplot(aes(x=HI,fill=Focal_Region,y=Boots))+
  geom_boxplot(alpha=0.5)+
  #geom_point(data=pvals,aes(x=HI,group=HI,y=Observed),pch=21,size=5,fill='darkviolet',position=position_dodge(width=1))+
  geom_text(data=cis,aes(x=HI,group=HI,y=Inf,label=Signif,size=Signif),vjust=1.5,position=position_dodge(width=1))+
  scale_size_manual(values=c(5,4))+
  theme_bw()+
  geom_hline(yintercept=0,lty=2)+
  scale_fill_manual(values=c('darkorchid2','grey60'))+
  ylim(c(0,1))+
  facet_grid(Variable ~., scales='free')+
  xlab('')+ylab('Relative Methylation Status')
bp

pdf('Methylation_levels_Bootstraps.pdf',height=6,width=8)
bp
dev.off()


#HI, plot raw 
hzm3$HI = factor(hzm3$HI,levels=c('C.c. corone','C.c. corone-like','Hybrid','C.c. cornix-like','C.c. cornix'))
mcp = hzm3 %>% 
  mutate(ord = fct_reorder(name,Chr18_Hybrid_Index)) %>% 
  filter(TaxonLone == 'Taxon' & Focal_Region == 'Focal') %>% 
  ggplot(aes(x=ord,y=value,fill=HI))+
  geom_violin(alpha=0.5)+
  geom_jitter()+
  facet_grid(TaxonLone ~ HI,scales='free')+
  scale_fill_grey()+xlab('')+ylab('Methylation Level (%) In Focal Region CpGs')+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
mcp
pdf('5mC_DistributionsRaw.pdf',height=5,width=8)
mcp 
dev.off()

```



## Classification with Year

```R
.libPaths('~/mambaforge/envs/randomforest/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/crow_scratch/2021-07_HybridZone/Nov/Classification_Sept')
library(dplyr)
library(tidyverse)
library(data.table)

#thresholds
pval <- 0.1
div <- 0.25
dd <- 0.1

#### CG
datsCG <- read.table('CG.DSS.BP-10x.txt',header=TRUE)
datsCG <- datsCG %>% dplyr::rename(fdrs_Species.cg = fdrs_SpeciesS.cg, fdrs_Sex.cg = fdrs_SexM.cg)

#with interacting effects, exclude variables significant for other variables
CG <- datsCG %>% mutate(DMP = ifelse(fdrs_Species.cg < pval & abs(divergence_Species.cg) > div &
                                       fdrs_Sex.cg > pval & fdrs_StageCHK.cg > pval & fdrs_StageYRL.cg > pval, 'Taxon',  #species
                                     ifelse(fdrs_Sex.cg < pval & abs(divergence_Sex.cg) > div &
                                              fdrs_Species.cg > pval & fdrs_StageCHK.cg > pval & fdrs_StageYRL.cg > pval, 'Sex', #sex
                                            ifelse((fdrs_StageCHK.cg < pval | fdrs_StageYRL.cg < pval) & (abs(divergence_StageCHK.cg) > div | abs(divergence_StageYRL.cg) > div) &
                                                     (fdrs_Sex.cg > pval & fdrs_Species.cg > pval),'Age',  #stage
                                                   ifelse(Divergence.cg < dd,  'Conserved','Unknown')))))

CG %>% dplyr::count(DMP)

#### WGBS
datsWGBS <- read.table('WGBS.DSS.BP-5x.txt',header=TRUE)

#with interacting effects, exclude variables significant for other variables
WGBS <- datsWGBS %>% mutate(DMP = ifelse(fdrs_Species.wgbs < pval & abs(divergence_Species.wgbs) > div &
                                           fdrs_Sex.wgbs > pval & fdrs_TissueM.wgbs > pval & fdrs_TissueL.wgbs > pval, 'Taxon',  #species
                                         ifelse(fdrs_Sex.wgbs < pval & abs(divergence_Sex.wgbs) > div &
                                                  fdrs_Species.wgbs > pval & fdrs_TissueM.wgbs > pval & fdrs_TissueL.wgbs > pval, 'Sex', #sex
                                                ifelse((fdrs_TissueM.wgbs < pval | fdrs_TissueL.wgbs < pval) & (abs(divergence_TissueM.wgbs) > div | abs(divergence_TissueL.wgbs) > div) &
                                                         (fdrs_Sex.wgbs > pval & fdrs_Species.wgbs > pval),'Tissue',  #tissue
                                                       ifelse(Divergence.wgbs < dd,  'Conserved','Unknown')))))


WGBS %>% dplyr::count(DMP)

#merged
dmr1 <- merge(WGBS,CG,by=c('site','chr','start'),
              suffixes=c('.wgbs','.cg'),all=TRUE)
#merge with HZ data..
hz <- read.table('HZ.DSS.BP-10x_YEAR.txt',header=T)
hz = hz %>% mutate(DMP.hz = ifelse(fdrs_Chr18_Hybrid_Index.hz < pval & fdrs_Distance.hz > pval & fdrs_Temperature.hz > pval & fdrs_Year2013.hz > pval & fdrs_Year2014.hz > pval, 'Taxon',
                                   ifelse(fdrs_Chr18_Hybrid_Index.hz > pval & fdrs_Distance.hz < pval & fdrs_Temperature.hz > pval & fdrs_Year2013.hz > pval & fdrs_Year2014.hz > pval, 'Distance',
                                          ifelse(fdrs_Chr18_Hybrid_Index.hz > pval & fdrs_Distance.hz > pval & fdrs_Temperature.hz < pval & fdrs_Year2013.hz > pval & fdrs_Year2014.hz > pval, 'Temperature',
                                                 ifelse(fdrs_Year2013.hz < pval & fdrs_Year2014.hz < pval & fdrs_Chr18_Hybrid_Index.hz > pval & fdrs_Distance.hz > pval & fdrs_Temperature.hz > pval, 'Season',
                                                 ifelse(Divergence.hz < dd, 'Conserved','Unknown'))))))
hz %>% dplyr::count(DMP.hz)
dmr2 <- merge(dmr1,hz,by=c('site','chr','start'),all=TRUE)

#count NAs by experiment
nac = dmr2 %>% select(site,DMP.wgbs,DMP.cg,DMP.hz)
nac$NAs = rowSums(is.na(nac))
nac %>% dplyr::count(NAs)

#replace NA with missing
dmr3 = dmr2 %>% replace_na(list(DMP.wgbs = 'MISSING',DMP.cg ='MISSING',DMP.hz = 'MISSING'))
#now, classify based on consensus on the common garden: 
dmr4 = dmr3 %>% mutate(DMP_Strict = ifelse(DMP.wgbs == 'Tissue' & 
                                             ((DMP.cg == 'Conserved' & DMP.hz == 'MISSING') | 
                                                (DMP.hz == 'Conserved' & DMP.cg == 'MISSING') |
                                                (DMP.cg == 'Conserved' & DMP.hz == 'Conserved')), 'Tissue',
                                           ifelse(DMP.cg == 'Age' & 
                                                    ((DMP.wgbs == 'Conserved' & DMP.hz == 'MISSING') | 
                                                       (DMP.hz == 'Conserved' & DMP.wgbs == 'MISSING' ) |
                                                       (DMP.wgbs == 'Conserved' & DMP.hz == 'Conserved')), 'Age',
                                                  ifelse(DMP.hz == 'Conserved' & 
                                                           ((DMP.cg == 'Sex' & DMP.wgbs == 'Sex') | 
                                                              ( DMP.cg == 'MISSING' & DMP.wgbs == 'Sex') |
                                                              (DMP.cg == 'Sex' & DMP.wgbs == 'MISSING' )), 'Sex',
                                                         ifelse(DMP.hz == 'Taxon' & 
                                                                  ((DMP.cg == 'Taxon' & DMP.wgbs == 'Taxon') | 
                                                                     ( DMP.cg == 'MISSING' & DMP.wgbs == 'Taxon') |
                                                                     (DMP.cg == 'Taxon' & DMP.wgbs == 'MISSING' )), 'Taxon',
                                                                ifelse(DMP.hz == 'Temperature' & 
                                                                         ((DMP.cg == 'Conserved' & DMP.wgbs == 'MISSING' ) | 
                                                                            (DMP.wgbs == 'Conserved' & DMP.cg == 'MISSING') |
                                                                            (DMP.cg == 'Conserved' & DMP.wgbs == 'Conserved')), 'Temperature',
                                                                       ifelse(DMP.hz == 'Season' & 
                                                                                ((DMP.cg == 'Conserved' & DMP.wgbs == 'MISSING' ) | 
                                                                                   (DMP.wgbs == 'Conserved' & DMP.cg == 'MISSING') |
                                                                                   (DMP.cg == 'Conserved' & DMP.wgbs == 'Conserved')), 'Season',
                                                                       ifelse(DMP.hz == 'Distance' & 
                                                                                ((DMP.cg == 'Conserved' & DMP.wgbs == 'MISSING' ) | 
                                                                                   (DMP.wgbs == 'Conserved' & DMP.cg == 'MISSING') |
                                                                                   (DMP.cg == 'Conserved' & DMP.wgbs == 'Conserved')), 'Distance',
                                                                              ifelse(DMP.cg == 'Conserved' & DMP.wgbs == 'Conserved', 'Conserved',
                                                                                     'Unclassified')))))))))
dmr4 = dmr4 %>% mutate(DMP_Single = ifelse(DMP.wgbs == 'Tissue' & 
                                             (DMP.cg == 'Conserved' | DMP.cg == 'Unknown' | DMP.cg == 'MISSING') &
                                             (DMP.hz == 'Conserved' | DMP.hz == 'Unknown' | DMP.hz == 'MISSING'), 'Tissue',
                                           ifelse(DMP.cg == 'Age' & 
                                                    (DMP.wgbs == 'Conserved' | DMP.wgbs == 'Unknown' | DMP.wgbs == 'MISSING') &
                                                    (DMP.hz == 'Conserved' | DMP.hz == 'Unknown' | DMP.hz == 'MISSING'), 'Age',
                                                  ifelse((DMP.hz == 'Conserved' | DMP.hz == 'Unknown' | DMP.hz == 'MISSING') & 
                                                           ((DMP.cg == 'Sex' & DMP.wgbs == 'Sex') | 
                                                              (( DMP.cg == 'MISSING' | DMP.cg == 'Unknown') & DMP.wgbs == 'Sex') |
                                                              (DMP.cg == 'Sex' & ( DMP.wgbs == 'MISSING' | DMP.wgbs == 'Unknown'))), 'Sex',
                                                         ifelse((DMP.hz == 'Taxon' | DMP.hz == 'Unknown' | DMP.hz == 'MISSING') & 
                                                                  ((DMP.cg == 'Taxon' & DMP.wgbs == 'Taxon') | 
                                                                     (( DMP.cg == 'MISSING' | DMP.cg == 'Unknown') & DMP.wgbs == 'Taxon') |
                                                                     (DMP.cg == 'Taxon' & ( DMP.wgbs == 'MISSING' | DMP.wgbs == 'Unknown'))), 'Taxon',
                                                                ifelse(DMP.hz == 'Temperature' & 
                                                                         (DMP.cg == 'Conserved' | DMP.cg == 'Unknown' | DMP.cg == 'MISSING') &
                                                                         (DMP.wgbs == 'Conserved' | DMP.wgbs == 'Unknown' | DMP.wgbs == 'MISSING'),'Temperature',
                                                                       ifelse(DMP.hz == 'Distance' & 
                                                                                (DMP.cg == 'Conserved' | DMP.cg == 'Unknown' | DMP.cg == 'MISSING') &
                                                                                (DMP.wgbs == 'Conserved' | DMP.wgbs == 'Unknown' | DMP.wgbs == 'MISSING'), 'Distance',
                                                                              ifelse((DMP.cg == 'Conserved' & (DMP.wgbs == 'Conserved' | DMP.wgbs == 'MISSING')) |
                                                                                       (DMP.wgbs == 'Conserved' & (DMP.cg == 'Conserved' | DMP.cg == 'MISSING')) , 'Conserved',
                                                                                     'Unclassified'))))))))

dmr4 %>% dplyr::count(DMP_Strict)
dmr4 %>% dplyr::count(DMP_Single)
dmr4 %>% group_by(DMP_Strict,DMP.wgbs,DMP.cg,DMP.hz) %>% count() %>% data.frame

dmr5 = dmr4 
dmr5$Mean5mC.cg = rowMeans(dmr5[grepl('^D_.*cg|^S_.*cg',names(dmr5))],na.rm=TRUE)
dmr5$Mean5mC.wgbs = rowMeans(dmr5[grepl('^D_.*wgbs|^S_.*wgbs',names(dmr5))],na.rm=TRUE)
dmr5$Mean5mC.hz = rowMeans(dmr5[grepl('^D_.*hz|^S_.*hz|^H_.*hz',names(dmr5))],na.rm=TRUE)                                                     

write.table(dmr5,'Full.DSS.BP-10x.Classified-23AUG30_YEAR.txt',quote=F,sep='\t',row.names=FALSE)
```

And merge with our REGION, DXY, FST, PI, data: 

```bash
rgndir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/crow_scratch/2021-07_HybridZone/Nov/Regions
gc=$rgndir/GCF_000738735.5_ASM73873v5_genomic.GC-Content.500bp.bed
rgn=$rgndir/GCF_000738735.5_ASM73873v5_genomic.GenomicFeaturesCDS.bed
popdat=$rgndir/GCF_000738735.5_ASM73873v5_genomic.PopGen2023JAN.bed

#create sites file 
sed '1d' Full.DSS.BP-10x.Classified-23AUG30_YEAR.txt | awk '{OFS="\t"}{print $2, $3, $3, $1}' > Sites_Year.bed

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



## Summarizing Counts

Summarized from table S5:

```R
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/crow_scratch/2021-07_HybridZone/Nov')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(scales)
library(ggpubr)

c = read_tsv('Classification_Numbers.txt')
counts = c %>% select(!c(Class,ComGar.WGBS,ComGar.RRBS,HybZon.RRBS)) %>% 
  pivot_longer(!c(Count)) %>%
  filter(value == 1) %>% group_by(name) %>% summarize(sum = sum(Count))
names(counts) = c('Filter','CpGs')
counts$Filter = factor(counts$Filter,levels=rev(c('All CpGs','Sequenced In At Least 2 Experiments','Passes Necessary Condition','No Conflicted DMP','Passes Sufficient Condition')))

#separate first and the rest
c1 = counts %>% filter(Filter == 'All CpGs') %>% 
  mutate(Filter = gsub('All CpGs','                                               All CpGs',Filter)) %>% 
  ggplot(aes(x = CpGs, y = reorder(Filter, CpGs),label=CpGs)) +
  geom_text(hjust=-0.5)+
  geom_segment(aes(xend = 0, yend = Filter), size = 1.5, color = "grey") +
  geom_point(size = 5, fill= "blue",pch=21) + xlab('')+ylab('')+
  scale_x_continuous(limits=c(0,6.1e6),
                     labels = scales::label_number_si(scale = 1e-6, suffix = "M"),
                     breaks = scales::pretty_breaks(n = 5)) +
  theme_minimal()
c1
c2 = counts %>% filter(Filter != 'All CpGs') %>% 
  mutate(Pct = paste0(CpGs,' (',round(CpGs/63073,3)*100,'%)'),
         Pct = ifelse(Filter == 'Sequenced In At Least 2 Experiments',CpGs,Pct)) %>% 
  ggplot(aes(x = CpGs, y = Filter,label=Pct)) +
  geom_text(hjust=-0.35)+
  geom_segment(aes(xend = 0, yend = Filter), size = 1.5, color = "grey") +
  geom_point(size = 5,pch=21,fill='blue') + xlab('')+ylab('')+
  scale_fill_manual(values=c('red','blue'))+
  scale_x_continuous(limits=c(0,78e4),
                     labels = scales::label_number_si(scale = 1e-3, suffix = "K"),
                     breaks = scales::pretty_breaks(n = 5)) +
  theme_minimal()
c2 

pdf('Multiexperimental_Counts_2023OCT09.pdf',height=6,width=8)
ggarrange(c1,c2,nrow=2,ncol=1,heights = c(0.25,0.75))
dev.off()

```

## UpSet Plot

```bash
.libPaths('~/mambaforge/envs/r/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/crow_scratch/2021-07_HybridZone/Nov/Classification_Sept')
library(tidyverse)
library(UpSetR)
options(scipen=999)
all = read_tsv('Full.DSS.BP-10x.Classified-Annotated_23FEB09.txt')

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

fromExpression(expressionInput)
up = upset(fromExpression(upset_data), order.by = "freq")

pdf('../figures/Upset_Taxon_Plot_2024APR30.pdf',height=4,width=7)
up
dev.off()

#save those 66 sites
targets = ap %>% filter(ComGar.WGBS == 'Taxon' & ComGar.RRBS == 'Taxon') %>% data.frame %>% select(site) %>% separate(site, into=c('chr','start')) %>% mutate(end=start)
write.table(targets,file='ComGar_Overlap_Taxon_2024APR30.txt',quote=F,sep='\t',row.names=F,col.names=F)

#of those 66, is chr18 still enriched?
targets %>% mutate(Focal_Region = ifelse(chr == 'chr18' & start > 8070001 & start < 10070000,'Focal','Undiverged'))

```


