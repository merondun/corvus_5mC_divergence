# Ordinations

## Methylation Ordinations

5mC: dbRDA 

```R
#### dbRDA
.libPaths('~/mambaforge/envs/randomforest/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_Sept/')
library(dplyr)
library(tidyverse)
library(tidymodels)
library(viridis)
library(ggplot2)
library(scales)
library(vegan)
library(ecodist)
library(ggord)
library(ggtext)
library(ggpubr)
library(GGally)

# Import 5mC data
class <- read.table('Full.DSS.BP-10x.Classified-Annotated_20250107.txt',header=TRUE)
class <- subset(class, Region != 'Missing' & Region != '.')

##### WGBS #####
# signal experiment 
xp <- '.wgbs'

# Select data, round to 1%
pcd <- class %>% select(contains(xp)) %>% select(-contains(c('fdrs','DMP','Mean','divergence','Variance','Divergence'))) %>% na.omit()
pcdd <- t(pcd)
pcr <- round(pcdd,2)

# grab ID
variables <- as.data.frame(names(pcd))
names(variables) <- 'ID'
variables$BAK <- variables$ID
variables$BAK <-  gsub(xp,'',variables$BAK)
variables <- variables %>% separate(BAK, into = c('Taxon','ID1','ID2','Tissue','Age','Sex'))
variables <- variables %>% select(ID,Taxon,Tissue,Sex)
variables <- variables %>% mutate(Taxon = ifelse(Taxon == 'D','C.C.','H.C.'),
                                  Tissue = ifelse(Tissue == 'L','Liver',
                                                  ifelse(Tissue == 'M','Spleen','Blood')),
                                  Sex = ifelse(Sex == 'F','Female','Male'))
#binarize in vegan format
rda_recipe =  # which consists of the formula (outcome ~ predictors)
  recipe(~ Taxon + Tissue + Sex, data = variables) %>%
  # and some pre-processing steps
  step_novel(all_nominal(), -all_outcomes()) %>% # converts all nominal variables to factors and takes care of other issues related to categorical variables
  step_dummy(all_nominal(), -all_outcomes()) %>% # converts our factor column ocean_proximity into numeric binary (0 and 1) variables
  step_zv(all_numeric(), -all_outcomes()) %>% # removes any numeric variables that have zero variance
  step_corr(all_predictors(), threshold = 0.8, method = "spearman") # will remove predictor variables that have large correlations with other predictor variables
proc = rda_recipe %>% # apply the recipe to the training data
  prep(variables) %>% juice()

#identify which distance is superior
rankindex(proc , pcr, indices = c("euc", "man", "gow", "bra", "kul"),stepacross= FALSE, method = "spearman")
# euc       man       gow       bra       kul 
# 0.6293556 0.6497516 0.6452239 0.6147569 0.6133181 

#constrained ordination with step selection 
wgbsnull = dbrda(pcr ~ 1, proc, dist="gow",scaling=TRUE)  # Model with intercept only
wgbsdb = dbrda(pcr ~ ., proc, dist="gow",scaling=TRUE) # Model with all explanatory variables

## With scope present, the default direction is "both"
#ordistep(wgbsnull, scope = formula(wgbsdb), perm.max = 200)

#by terms
tmw <- as.data.frame(anova(wgbsdb, by="terms", permu=10000)) # test for sign. environ. variables
tmw
# Df   SumOfSqs        F Pr(>F)
# Taxon_H.C.     1 0.06361369 1.602512  0.106
# Tissue_Liver   1 0.15342459 3.864965  0.001
# Tissue_Spleen  1 0.13395068 3.374392  0.001
# Sex_Male       1 0.05839425 1.471027  0.100
# Residual       7 0.27787369       NA     NA

# adjusted R^2
R2adjw <- round(RsquareAdj(wgbsdb)$adj.r.squared,3)
pvalw <- anova(wgbsdb) # overall test of the significant of the analysis
plabw <- round(pvalw[1,4],3)
plabw
# [1] 0.001
R2adjw
# [1] 0.365

#plot raw
p1 <- ggord(wgbsdb,variables$Taxon,ellipse=F,veccol="blue",labcol="blue",addcol="black",addsize=0.0000001,addpch=27,coord_fix=FALSE)+xlab('dbRDA1')+ylab('dbRDA2')
p1$layers[[1]] <- NULL
new_dat <- data.frame(p1$data$one,p1$data$two,variables$Taxon,variables$Sex,variables$Tissue)
new_dat$variables.SEX_SP <- paste0(new_dat$variables.Taxon,'_',new_dat$variables.Sex)
new_dat$variables.SEX_SP <- factor(new_dat$variables.SEX_SP, levels=c('C.C._Female','C.C._Male','H.C._Female','H.C._Male'))
names(new_dat) <- gsub('.*\\.','',names(new_dat))

# add points layer manually, redefine some aesthetics
wgbsp <- p1 + geom_point(data = new_dat, show.legend = F,aes(x=one,
                                                             y=two, 
                                                             shape=SEX_SP,
                                                             colour = Tissue,
                                                             fill = Tissue,
                                                             size=10)) +
  scale_color_manual(values=rev(viridis(3,option='turbo')))+
  #geom_text(data=new_dat,aes(x=one,y=two,label=Sex))+
  scale_fill_manual(values=rev(viridis(3,option='turbo')))+
  scale_shape_manual(values=c(1,21,2,24))+
  geom_richtext(aes(x=Inf,y=Inf,label=paste("adj. R<sup>2</sup>: ",round(R2adjw,3))),vjust=1,hjust=1,size=5,label.colour=NA)+
  theme_classic(base_size=12)+
  geom_richtext(aes(x=Inf,y=Inf,label=paste("p = ",plabw)),vjust=2,hjust=1,size=5,label.colour=NA)+
  expand_limits(x = c(min(new_dat[,1])+min(new_dat[,1])*2,
                      max(new_dat[,1])+max(new_dat[,1])*.5),
                y = c(min(new_dat[,2])+min(new_dat[,2])*2,
                      max(new_dat[,2])+max(new_dat[,2])*2))+xlab('dbRDA1')+ylab('dbRDA2')
wgbsp

pdf('20250107_RDA_WGBS-Small.pdf', height=4, width=4.5)
wgbsp
dev.off()  

#legend
shp = unique(new_dat %>% select(c(SEX_SP))) %>% mutate(Shape = c(21,1,24,2), Color = 'black') %>% dplyr::rename(Legend = SEX_SP)
col = unique(new_dat %>% select(c(Tissue))) %>% mutate(Color = rev(viridis(3,option='turbo')),Shape=21) %>% dplyr::rename(Legend = Tissue)
allw = rbind(shp,col) %>% mutate(x=row_number(),y=1)
allw = allw %>% mutate(Legend = gsub('C.C._Female', 'C.C. Female',Legend),
                       Legend = gsub('C.C._Male', 'C.C. Male',Legend),
                       Legend = gsub('H.C._Female', 'H.C. Female',Legend),
                       Legend = gsub('H.C._Male', 'H.C. Male',Legend),
                       Legend = gsub('ADL', 'Adult',Legend))
allw$Legend <- factor(allw$Legend, levels=c('C.C. Female','C.C. Male','H.C. Female','H.C. Male','Blood','Liver','Spleen'))

pdf('20250107_RDA_WGBS-Legend.pdf', height=6, width=8)
plot.new()
legend(x = "top",inset=c(0,-.1),xpd=TRUE,pt.cex=1.5,legend=allw$Legend,pt.bg=allw$Color,pch=allw$Shape,ncol=4,cex=.8)
dev.off()  

##### CG #####
# signal experiment 
xp <- '.cg'

#select data, round to 1%, and change 0% 5mC to 1% 5mC for CCA
pcd <- class %>% select(contains(xp)) %>% select(-contains(c('fdrs','DMP','Mean','divergence','Variance','Divergence'))) %>% na.omit()
pcdd <- t(pcd)
pcr <- round(pcdd,2)

#grab ID
variables <- as.data.frame(names(pcd))
names(variables) <- 'ID'
variables$BAK <- variables$ID
variables$BAK <-  gsub(xp,'',variables$BAK)
variables <- variables %>% separate(BAK, into = c('Taxon','ID1','ID2','Tissue','Age','Sex'))
variables <- variables %>% select(ID,Taxon,Age,Sex)
variables <- variables %>% mutate(Taxon = ifelse(Taxon == 'D','C.C.','H.C.'),
                                  Sex = ifelse(Sex == 'F','Female','Male'))
#binarize in vegan format
rda_recipe =  # which consists of the formula (outcome ~ predictors)
  recipe(~ Taxon + Age + Sex, data = variables) %>%
  # and some pre-processing steps
  step_novel(all_nominal(), -all_outcomes()) %>% # converts all nominal variables to factors and takes care of other issues related to categorical variables
  step_dummy(all_nominal(), -all_outcomes()) %>% # converts our factor column ocean_proximity into numeric binary (0 and 1) variables
  step_zv(all_numeric(), -all_outcomes()) %>% # removes any numeric variables that have zero variance
  step_corr(all_predictors(), threshold = 0.8, method = "spearman") # will remove predictor variables that have large correlations with other predictor variables
proc = rda_recipe %>% # apply the recipe to the training data
  prep(variables) %>% juice()

#identify which distance is superior
rankindex(proc , pcr, indices = c("euc", "man", "gow", "bra", "kul"),stepacross= FALSE, method = "spearman")
# euc       man       gow       bra       kul 
# 0.3956835 0.3697792 0.3743758 0.3488608 0.3584490 

#constrained ordination with step selection 
cgnull <- dbrda(pcr ~ 1, proc, dist="euc",scaling=TRUE)  # Model with intercept only
cgdb = dbrda(pcr ~ ., proc, dist="euc",scaling=TRUE) # Model with all explanatory variables

## With scope present, the default direction is "both"
#ordistep(cgnull, scope = formula(cgdb), perm.max = 200)

#by terms
tmc <- as.data.frame(anova(cgdb, by="terms", permu=10000)) # test for sign. environ. variables
tmc
# Df   Variance         F Pr(>F)
# Taxon_H.C.  1  127.78286 1.9789492  0.001
# Age_CHK     1  165.53272 2.5635743  0.001
# Age_YRL     1   49.35923 0.7644172  0.966
# Sex_Male    1  130.87982 2.0269113  0.001
# Residual   19 1226.85023        NA     NA

# adjusted R^2
R2adjc <- round(RsquareAdj(cgdb)$adj.r.squared,3)
pvalc <- anova(cgdb) # overall test of the significant of the analysis
plabc <- round(pvalc[1,4],3)
plabc
# [1] 0.001
R2adjc
# [1] 0.127

#plot raw
p2 <- ggord(cgdb,variables$Taxon,ellipse=F,veccol="blue",labcol="blue",addcol="black",addsize=0.0000001,addpch=27,coord_fix=FALSE)+xlab('dbRDA1')+ylab('dbRDA2')
p2$layers[[1]] <- NULL
new_dat <- data.frame(p2$data$one,p2$data$two,variables$Taxon,variables$Sex,variables$Age)
new_dat$variables.SEX_SP <- paste0(new_dat$variables.Taxon,'_',new_dat$variables.Sex)
new_dat$variables.SEX_SP <- factor(new_dat$variables.SEX_SP, levels=c('C.C._Female','C.C._Male','H.C._Female','H.C._Male'))
names(new_dat) <- gsub('.*\\.','',names(new_dat))

# add points layer manually, redefine some aesthetics
cgp <- p2 + geom_point(data = new_dat, show.legend = F,aes(x=one,
                                                           y=two, 
                                                           shape=SEX_SP,
                                                           colour = Age,
                                                           fill = Age,
                                                           size=10)) +
  scale_color_manual(values=rev(viridis(3,option='cividis')))+
  scale_fill_manual(values=rev(viridis(3,option='cividis')))+
  #geom_text(data=new_dat,aes(x=one,y=two,label=Sex))+
  scale_shape_manual(values=c(1,21,2,24))+
  geom_richtext(aes(x=Inf,y=Inf,label=paste("adj. R<sup>2</sup>: ",round(R2adjc,3))),vjust=1,hjust=1,size=5,label.colour=NA)+
  theme_classic(base_size=12)+
  geom_richtext(aes(x=Inf,y=Inf,label=paste("p = ",plabc)),vjust=2,hjust=1,size=5,label.colour=NA)+
  expand_limits(x = c(min(new_dat[,1])+min(new_dat[,1])*1,
                      max(new_dat[,1])+max(new_dat[,1])*1.5),
                y = c(min(new_dat[,2])+min(new_dat[,2])*1,
                      max(new_dat[,2])+max(new_dat[,2])*1.5))+xlab('dbRDA1')+ylab('dbRDA2')
cgp

pdf('20250107_RDA_CG-Small.pdf', height=4, width=4.5)
cgp
dev.off()  

#legend
shp = unique(new_dat %>% select(c(SEX_SP))) %>% mutate(order=c(2,1,3,4)) %>% arrange(order) %>% 
  mutate(Shape = c(21,1,24,2), Color = 'black') %>% dplyr::rename(Legend = SEX_SP) %>% select(-order)
col = unique(new_dat %>% select(c(Age))) %>% mutate(Color = rev(viridis(3,option='cividis')),Shape=21) %>% dplyr::rename(Legend = Age)
allc = rbind(shp,col) %>% mutate(x=row_number(),y=1)
allc = allc %>% mutate(Legend = gsub('C.C._Female', 'C.C. Female',Legend),
                       Legend = gsub('C.C._Male', 'C.C. Male',Legend),
                       Legend = gsub('H.C._Female', 'H.C. Female',Legend),
                       Legend = gsub('H.C._Male', 'H.C. Male',Legend),
                       Legend = gsub('ADL', 'Adult',Legend),
                       Legend = gsub('YRL', 'Yearling',Legend),
                       Legend = gsub('CHK', 'Chick',Legend),
                       order = c(2,1,3,4,7,5,6)) %>% 
  arrange(order)
allc$Legend <- factor(allc$Legend, levels=c('C.C. Male','C.C. Female','H.C. Female','H.C. Male','Chick','Adult','Yearling'))

pdf('20250107_RDA_CG-Legend.pdf', height=6, width=8)
plot.new()
legend(x = "top",inset=c(0,-.1),xpd=TRUE,pt.cex=1.5,legend=allc$Legend,pt.bg=allc$Color,pch=allc$Shape,ncol=4,cex=.8)
dev.off()  

##### HZ #####
# signal experiment 
xp <- '.hz'

#select data, round to 1%, and change 0% 5mC to 1% 5mC for CCA
pcd <- class %>% select(contains(xp)) %>% select(-contains(c('fdrs','DMP','Mean','divergence','Variance','Divergence'))) %>% na.omit()
pcdd <- t(pcd)
pcr <- round(pcdd,2)
rownames(pcr) <- gsub('.hz','',rownames(pcr))

#grab ID
variables <- read.table('All-Metadata.txt',header=T)
variables <- subset(variables,Experiment=='HZ') %>% select(-ID)
variables <- variables %>% rename('ID'='Short_ID','Year'='Capture_Year')
variables$Year <- as.factor(variables$Year)
keep <- as.data.frame(row.names(pcr))
names(keep) <- 'ID'
variables <- left_join(variables,keep)
variables <- variables[match(rownames(pcr), variables$ID), ]
rownames(variables) <- variables$ID

#pre-process data 
rda_recipe =  # which consists of the formula (outcome ~ predictors)
  recipe(~ Chr18_Hybrid_Index + Year, data = variables) %>%
  # and some pre-processing steps
  step_novel(all_nominal(), -all_outcomes()) %>% # converts all nominal variables to factors
  step_dummy(all_nominal(), -all_outcomes()) %>% # converts our factor column into numeric binary 
  step_zv(all_numeric(), -all_outcomes()) %>% # removes any numeric variables that have zero variance
  step_scale(all_numeric(), -all_outcomes()) %>% # removes any numeric variables that have zero variance
  step_corr(all_predictors(), threshold = 0.8, method = "spearman") # will remove predictor variables that have large correlations with other predictor variables
proc = rda_recipe %>% # apply the recipe to the training data
  prep(variables) %>% juice()
proc

#identify which distance is superior
rankindex(proc , pcr, indices = c("euc", "man", "gow", "bra", "kul"),stepacross= FALSE, method = "spearman")
# euc       man       gow       bra       kul 
# 0.1196425 0.1492669 0.2045072 0.1842789 0.1853819 

#constrained ordination with step selection 
hznull <- dbrda(pcr ~ 1, proc, dist="euc",scaling=TRUE)  # Model with intercept only
hzdb = dbrda(pcr ~ ., proc, dist="euc",scaling=TRUE) # Model with all explanatory variables

## With scope present, the default direction is "both"
#ordistep(hznull, scope = formula(hzdb), perm.max = 200)

#by terms
tmh <- as.data.frame(anova(hzdb, by="terms", permu=10000)) # test for sign. environ. variables
tmh 
# Df  Variance        F Pr(>F)
# Chr18_Hybrid_Index  1  121.7465 1.189927  0.044
# Year_X2013          1  151.6893 1.482582  0.001
# Year_X2014          1  102.7136 1.003903  0.394
# Residual           18 1841.6562       NA     NA

# adjusted R^2
R2adjh <- round(RsquareAdj(hzdb)$adj.r.squared,3)
pvalh <- anova(hzdb) # overall test of the significant of the analysis
plabh <- round(pvalh[1,4],3)
plabh
# [1] 0.004
R2adjh
# [1] 0.031

#plot raw
p3 <- ggord(hzdb,variables$Taxon,ellipse=F,veccol="blue",labcol="blue",addcol="black",addsize=0.0000001,addpch=27,coord_fix=FALSE)+xlab('dbRDA1')+ylab('dbRDA2')
p3$layers[[1]] <- NULL
new_dat <- data.frame(p3$data$one,p3$data$two,variables$Species,variables$Chr18_Hybrid_Index,variables$Year)
names(new_dat) <- gsub('.*\\.','',names(new_dat))
new_dat

# add points layer manually, redefine some aesthetics
hzp <- p3 + 
  geom_point(data = new_dat, show.legend = F,aes(x=one,
                                                 y=two, 
                                                 shape=Year,
                                                 fill = Chr18_Hybrid_Index,
                                                 size=10)) +
  scale_fill_gradient(low='grey30',high='grey90')+
  scale_shape_manual(values=c(21,24,23))+
  geom_richtext(aes(x=Inf,y=Inf,label=paste("adj. R<sup>2</sup>: ",round(R2adjh,3))),vjust=1,hjust=1,size=5,label.colour=NA)+
  geom_richtext(aes(x=Inf,y=Inf,label=paste("p = ",plabh)),vjust=2,hjust=1,size=5,label.colour=NA)+
  theme_classic(base_size=12)+
  expand_limits(x = c(min(new_dat[,1])+min(new_dat[,1])*1.5,
                      max(new_dat[,1])+max(new_dat[,1])*1.5),
                y = c(min(new_dat[,2])+min(new_dat[,2])*0.5,
                      max(new_dat[,2])+max(new_dat[,2])*1.5))+
  guides(fill=guide_legend(override.aes=list(shape=21)))

hzp

pdf('20250107_RDA_HZ.pdf', height=4, width=4.5)
hzp
dev.off()  

pdf('20250107_RDA_HZ-Legend.pdf', height=6, width=8)
hzp
dev.off() 

```



## Resequencing Ordinations

Resequencing PCA, NMDS: 

```r
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Crow_Population_Resequencing_ASM73873v5/bcftools/phased')
.libPaths('~/miniconda3/envs/snprelate/lib/R/library')
library(SeqArray)
library(SNPRelate)
library(vegan)
library(ggplot2)

# open a vcf
seqVCF2GDS('Background.vcf.gz','Background.gds')
gds <- seqOpen('Background.gds')

#assess missing data
RV <- snpgdsSampMissRate(gds)
summary(RV)

# run PCA on full dataset
RV <- snpgdsPCA(gds,autosome.only =F)

# eigenvalues
head(RV$eigenval)

# make a data.frame
tab <- data.frame(ID = RV$sample.id,
                  PC1 = RV$eigenvect[,1],    # the first eigenvector
                  PC2 = RV$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)
write.table(tab,file='Background_PCs.txt',quote=F,sep='\t',row.names=F)
write.table(RV$eigenval,file='Background_EVs.txt',quote=F,sep='\t',row.names=F)

#LD Prune 
set.seed(1000)
snpset <- snpgdsLDpruning(gds, autosome.only =F, slide.max.bp = 50000, ld.threshold=0.2)
snpset.id <- unlist(snpset, use.names=FALSE)  # get the variant IDs of a LD-pruned set
head(snpset.id)

#create new output
grm_fn <- "Background_prune.gds"
seqSetFilter(gds, variant.id=snpset.id)

# export to a GDS genotype file without annotation data
seqExport(gds, grm_fn, info.var=character(), fmt.var=character(), samp.var=character())

#import LD-pruned file 
gds <- seqOpen('Background_prune.gds')
summary(gds)

# run PCA on LD dataset
RV <- snpgdsPCA(gds,autosome.only =F)

# eigenvalues
head(RV$eigenval)

# make a data.frame
tab <- data.frame(ID = RV$sample.id,
                  PC1 = RV$eigenvect[,1],    # the first eigenvector
                  PC2 = RV$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)
write.table(tab,file='Background_LD-PCs.txt',quote=F,sep='\t',row.names=F)
write.table(RV$eigenval,file='Background_LD-EVs.txt',quote=F,sep='\t',row.names=F)

#read metadat
variables <- read.table('../k28.metdat',header=F)
names(variables) <- c('ID','Locality','Sex','Taxon')

#Plot original PCA 
tab = read.table('Background_PCs.txt',header=TRUE)
RV = read.table('Background_EVs.txt',header=TRUE)
tab <- merge(tab,variables,by='ID')

np <- 
  ggplot() +
  geom_point(data = tab, show.legend = T,aes(x=PC1,
                                             y=PC2, 
                                             shape = Locality,
                                             colour = Taxon,
                                             fill = Taxon,
                                             size=15)) +
  scale_color_manual(values=c('grey10','grey60'))+
  scale_fill_manual(values=c('grey10','grey60'))+
  scale_shape_manual(values=c(16,15,18,17))+
  xlab(paste0('PC1 : ',round(RV$x[1],2),'% Explained'))+
  ylab(paste0('PC2 : ',round(RV$x[2],2),'% Explained'))+
  theme_classic(base_size=14)
np

pdf('Reseq_PCA-Background.pdf',height=3,width=4)
np
dev.off()

gm <- snpgdsGetGeno(gds) 
gmut <- as.data.frame(gm)

#### NMDS ####
library(tidyverse)
set.seed(123)
gmt <- as.data.frame(t(gm))
names(gmt) <- RV$sample.id
gmut <- t(gmt)

nmds = metaMDS(gmut, distance = "bray", trymax = 100)
nmds
data.scores <- as.data.frame(scores(nmds)$sites)
data.scores$ID <- rownames(data.scores)
data.scores <- merge(data.scores,variables,by='ID')
data.scores = data.scores %>% mutate(Taxon = ifelse(grepl('^D_',ID)==TRUE,'C.C.','H.C.'),
                                      Stress = nmds$stress,
                                      Sites = nrow(gmt))
write.table(data.scores,file='Background_NMDS.txt',quote=F,sep='\t',row.names=F)

np <- 
  ggplot() +
  geom_point(data = data.scores, show.legend = T,aes(x=NMDS1,
                                                     y=NMDS2, 
                                                     shape = Locality,
                                                     colour = Taxon,
                                                     fill = Taxon,
                                                     size=15)) +
  scale_color_manual(values=c('grey10','grey60'))+
  scale_fill_manual(values=c('grey10','grey60'))+
  scale_shape_manual(values=c(16,15,18,17))+
  geom_text(aes(x=Inf,y=Inf,label=paste0('Stress: ',round(nmds$stress,3))),vjust=2.5,hjust=1,size=5)+
  theme_classic(base_size=14)
np

pdf('Background_NMDS.pdf', height=5, width=6)
np
dev.off()  

#perform anosim:
ord <- as.data.frame(row.names(gmut))
names(ord) <- 'ID'
vars <- left_join(data.frame(name=ord),data.scores,by="ID")
ano  <- anosim(gmut, vars$Taxon, distance = "bray", permutations = 9999)
anoreq <- data.frame(ano$statistic,ano$signif)
anoreq$variable <- 'Taxon'
names(anoreq) <- c('statistic','pval','variable')
anoreq$experiment <- 'Reseq'
write.table(anoreq,file='Background-ANOSIM.txt',quote=F,sep='\t',row.names=F)

```



