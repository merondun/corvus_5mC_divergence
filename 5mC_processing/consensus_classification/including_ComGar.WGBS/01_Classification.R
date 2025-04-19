.libPaths('~/miniconda3/envs/randomforest/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_Sept/')
library(tidyverse)
library(data.table)
library(viridis)

# Thresholds: FDR corrected pvalue, divergence between groups for factors (e.g. 25%), min% for conserved
pval <- 0.01
div <- 0.25
dd <- 0.1

#### CG
datsCG <- read.table('CG.DSS.BP-10x.txt',header=TRUE)
datsCG <- datsCG %>% dplyr::rename(fdrs_Species.cg = fdrs_SpeciesS.cg, fdrs_Sex.cg = fdrs_SexM.cg)

# With interacting effects, exclude variables significant for other variables
CG <- datsCG %>% mutate(DMP = ifelse(fdrs_Species.cg < pval & abs(divergence_Species.cg) > div &
                                       fdrs_Sex.cg > pval & fdrs_StageCHK.cg > pval & fdrs_StageYRL.cg > pval, 'Taxon',  #species
                                     ifelse(fdrs_Sex.cg < pval & abs(divergence_Sex.cg) > div &
                                              fdrs_Species.cg > pval & fdrs_StageCHK.cg > pval & fdrs_StageYRL.cg > pval, 'Sex', #sex
                                            ifelse((fdrs_StageCHK.cg < pval | fdrs_StageYRL.cg < pval) & (abs(divergence_StageCHK.cg) > div | abs(divergence_StageYRL.cg) > div) &
                                                     (fdrs_Sex.cg > pval & fdrs_Species.cg > pval),'Age',  #stage
                                                   ifelse(Divergence.cg < dd,  'Conserved','Unknown')))))

CG %>% dplyr::count(DMP)

#### WGBS
datsWGBS <- read.table('WGBS.DSS.BP-10x.txt',header=TRUE)

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
hz <- read.table('HZ.DSS.BP-Parentals-YearChr18HI_10x.txt',header=T)
hz = hz %>% mutate(DMP.hz = ifelse(fdrs_Chr18_Hybrid_Index.hz < pval & fdrs_Year2013.hz > pval & fdrs_Year2014.hz > pval, 'Taxon',
                                   ifelse(Divergence.hz < dd, 'Conserved','Unknown')))
hz %>% dplyr::count(DMP.hz)
dmr2 <- merge(dmr1,hz,by=c('site','chr','start'),all=TRUE)

#count NAs by experiment
nac = dmr2 %>% select(site,DMP.wgbs,DMP.cg,DMP.hz)
nac$NAs = rowSums(is.na(nac))
nac %>% dplyr::count(NAs)

#replace NA with missing
dmr3 = dmr2 %>% replace_na(list(DMP.wgbs = 'MISSING',DMP.cg ='MISSING',DMP.hz = 'MISSING'))

# Assign sites based on conflicts between experiments
dmr4 <- dmr3 %>% 
  mutate(
    Class = case_when(
      # Taxon classification
      (DMP.cg == "Taxon" & DMP.hz %in% c("Taxon", "Unknown", "Conserved", "MISSING") & DMP.wgbs %in% c("Taxon", "Unknown", "Conserved", "MISSING")) ~ "Taxon",
      (DMP.hz == "Taxon" & DMP.cg %in% c("Taxon", "Unknown", "Conserved", "MISSING") & DMP.wgbs %in% c("Taxon", "Unknown", "Conserved", "MISSING")) ~ "Taxon",
      (DMP.wgbs == "Taxon" & DMP.hz %in% c("Taxon", "Unknown", "Conserved", "MISSING") & DMP.cg %in% c("Taxon", "Unknown", "Conserved", "MISSING")) ~ "Taxon",
      
      # Tissue classification
      (DMP.wgbs == "Tissue" & DMP.hz %in% c("Unknown", "Conserved", "MISSING") & DMP.cg %in% c("Unknown", "Conserved", "MISSING")) ~ "Tissue",
      
      # Age classification
      (DMP.cg == "Age" & DMP.hz %in% c("Unknown", "Conserved", "MISSING") & DMP.wgbs %in% c("Unknown", "Conserved", "MISSING")) ~ "Age",
      
      # Sex classification
      (DMP.cg == "Sex" & DMP.hz %in% c("Unknown", "Conserved", "MISSING") & DMP.wgbs %in% c("Sex", "Unknown", "Conserved", "MISSING")) ~ "Sex",
      (DMP.wgbs == "Sex" & DMP.hz %in% c("Unknown", "Conserved", "MISSING") & DMP.cg %in% c("Sex", "Unknown", "Conserved", "MISSING")) ~ "Sex",
      
      # Conserved classification
      (DMP.cg == "Conserved" & DMP.hz == "Conserved") ~ "Conserved",
      
      # Conflict
      TRUE ~ "Unknown"
    ),
    # Count experiments supporting each classification
    SupportCount = case_when(
      Class == "Taxon" ~ rowSums(cbind(DMP.cg == "Taxon", DMP.wgbs == "Taxon", DMP.hz == "Taxon")),
      Class == "Tissue" ~ as.integer(DMP.wgbs == "Tissue"),
      Class == "Age" ~ as.integer(DMP.cg == "Age"),
      Class == "Sex" ~ rowSums(cbind(DMP.cg == "Sex", DMP.wgbs == "Sex")),
      Class == "Conserved" ~ 3, # All three experiments must agree
      TRUE ~ 0 # For conflicts or unknown classifications
    )
  )

dmr4 <- dmr4 %>% mutate(Focal_Region = ifelse(chr == 'chr18' & start > 8.07e6 & start < 10.07e6,'Focal','Undiverged'))
dmr4 %>% count(SupportCount,Class)
dmr4 %>% filter(Class == 'Taxon' & SupportCount >= 2 & DMP.hz == 'Taxon')
dmr4 %>% filter(Class == 'Taxon' & SupportCount >= 2 & DMP.hz == 'Taxon') %>% select(matches('site|Focal|fdrs')) %>% select(matches('site|Focal|Species|Hybrid'))


dmr4 %>% filter(chr == 'chr18' & start %in% c(10000910, 9005914, 9005915)) %>% select(site,contains('.hz'))
dmr5 = dmr4 
dmr5$Mean5mC.cg = rowMeans(dmr5[grepl('^D_.*cg|^S_.*cg',names(dmr5))],na.rm=TRUE)
dmr5$Mean5mC.wgbs = rowMeans(dmr5[grepl('^D_.*wgbs|^S_.*wgbs',names(dmr5))],na.rm=TRUE)
dmr5$Mean5mC.hz = rowMeans(dmr5[grepl('^D_.*hz|^S_.*hz|^H_.*hz',names(dmr5))],na.rm=TRUE)                                                     

write.table(dmr5,'Full.DSS.BP-10x.Classified-20250107.txt',quote=F,sep='\t',row.names=FALSE)

