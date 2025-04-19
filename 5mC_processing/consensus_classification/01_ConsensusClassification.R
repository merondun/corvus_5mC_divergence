setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_202504/')
.libPaths('~/mambaforge/envs/randomforest/lib/R/library')
library(tidyverse)
library(data.table)
library(viridis)

# Thresholds: FDR corrected pvalue, divergence between groups for factors (e.g. 25%), min% for conserved
pval <- 0.01
div <- 0.25
dd <- 0.1

#### CG
datsCG <- read.table('CG_DSS.10x.txt',header=TRUE)
datsCG <- datsCG %>% dplyr::rename(fdrs_Species.CG = fdrs_SpeciesS.CG, fdrs_Sex.CG = fdrs_SexM.CG)

# With interacting effects, exclude variables significant for other variables
cg <- datsCG %>% mutate(DMP.CG = ifelse(fdrs_Species.CG < pval & abs(divergence_Species.CG) > div &
					                                          fdrs_Sex.CG > pval & fdrs_StageCHK.CG > pval & fdrs_StageYRL.CG > pval, 'Taxon',  #species
									                                          ifelse(fdrs_Sex.CG < pval & abs(divergence_Sex.CG) > div &
															                                                  fdrs_Species.CG > pval & fdrs_StageCHK.CG > pval & fdrs_StageYRL.CG > pval, 'Sex', #sex
																				                                                 ifelse((fdrs_StageCHK.CG < pval | fdrs_StageYRL.CG < pval) & (abs(divergence_StageCHK.CG) > div | abs(divergence_StageYRL.CG) > div) &
																											                                                        (fdrs_Sex.CG > pval & fdrs_Species.CG > pval),'Age',  #stage
																																	                                                      ifelse(Divergence.CG < dd,  'Conserved','Unknown')))))

cg %>% dplyr::count(DMP.CG)


#merge with HZ data..
hz <- read.table('HZ_DSS.10x.txt',header=T)
hz = hz %>% mutate(DMP.HZ = ifelse(fdrs_Chr18_Hybrid_Index.HZ < pval & fdrs_Year2013.HZ > pval & fdrs_Year2014.HZ > pval, 'Taxon',
				                                      ifelse(Divergence.HZ < dd, 'Conserved','Unknown')))
hz %>% dplyr::count(DMP.HZ)
dmr2 <- merge(cg,hz,by=c('site','chr','start'),all=TRUE)

#count NAs by experiment
nac = dmr2 %>% select(site,DMP.CG,DMP.HZ)
nac$NAs = rowSums(is.na(nac))
nac %>% dplyr::count(NAs)

#replace NA with missing
dmr3 = dmr2 %>% replace_na(list(DMP.CG ='MISSING',DMP.HZ = 'MISSING'))

# Assign sites based on conflicts between experiments
dmr4 <- dmr3 %>% 
	  mutate(
		     Class = case_when(
				             # Taxon classification
				             (DMP.CG == "Taxon" & DMP.HZ %in% c("Taxon", "Unknown", "Conserved", "MISSING")) ~ "Taxon",
					           (DMP.HZ == "Taxon" & DMP.CG %in% c("Taxon", "Unknown", "Conserved", "MISSING")) ~ "Taxon",
					           
					           # Age classification
					           (DMP.CG == "Age" & DMP.HZ %in% c("Unknown", "Conserved", "MISSING"))  ~ "Age",
						         
						         # Sex classification
						         (DMP.CG == "Sex" & DMP.HZ %in% c("Unknown", "Conserved", "MISSING")) ~ "Sex",
						         
						         # Conserved classification
						         (DMP.CG == "Conserved" & DMP.HZ == "Conserved") ~ "Conserved",
							       
							       # Conflict
							       TRUE ~ "Unknown"
							     ),
		     # Count experiments supporting each classification
		     SupportCount = case_when(
					            Class == "Taxon" ~ rowSums(cbind(DMP.CG == "Taxon", DMP.HZ == "Taxon")),
						          Class == "Age" ~ as.numeric(DMP.CG == "Age"),
						          Class == "Sex" ~ as.numeric(DMP.CG == "Sex"),
							        Class == "Conserved" ~ 2, # All three experiments must agree
							        TRUE ~ 0 # For conflicts or unknown classifications
								    )
		       )

dmr4 <- dmr4 %>% mutate(Focal_Region = ifelse(chr == 'chr18' & start > 8.07e6 & start < 10.07e6,'Focal','Undiverged'))
dmr4 %>% count(SupportCount,Class)
dmr4 %>% filter(Class == 'Taxon' & SupportCount >= 2 & DMP.HZ == 'Taxon')
dmr4 %>% filter(Class == 'Taxon' & SupportCount >= 2 & DMP.HZ == 'Taxon') %>% select(matches('site|Focal|fdrs')) %>% select(matches('site|Focal|Species|Hybrid'))


dmr4 %>% filter(chr == 'chr18' & start %in% c(10000910, 9005914, 9005915)) %>% select(site,contains('.HZ'))
dmr5 = dmr4 
dmr5$Mean5mC.CG = rowMeans(dmr5[grepl('^D_.*CG|^S_.*CG',names(dmr5))],na.rm=TRUE)
dmr5$Mean5mC.HZ = rowMeans(dmr5[grepl('^D_.*HZ|^S_.*HZ|^H_.*HZ',names(dmr5))],na.rm=TRUE)                                                     

write.table(dmr5,'Full-HZ-CG.DSS.BP-10x.Classified-20250410.txt',quote=F,sep='\t',row.names=FALSE)

