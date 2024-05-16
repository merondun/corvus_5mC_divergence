# Methylation Ordinations

Ordinations for ComGarRRBS(CG), ComGarWGBS(WGBS), and HybZonRRBS(HZ). PCA and NMDS were performed on each. Additional sensitivities incorporating a year effect included below. PDF plots, including small versions included in the main figures, and PDFs used for legends, included. 

5mC: dbRDA + NMDS 

```R
#### PCA, NMDS, db-RDA
.libPaths('~/mambaforge/envs/randomforest/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/crow_scratch/2021-07_HybridZone/Nov/Classification_Sept')
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

#import 5mC data
class <- read.table('Full.DSS.BP-10x.Classified-Annotated_23FEB09.txt',header=TRUE)
class <- subset(class, Region != 'Missing' & Region != '.')

##### WGBS #####
#signal experiment 
xp <- '.wgbs'

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
# 0.5968575 0.6059976 0.6099752 0.5869134 0.5823433 

#constrained ordination with step selection 
wgbsnull = dbrda(pcr ~ 1, proc, dist="gow",scaling=TRUE)  # Model with intercept only
wgbsdb = dbrda(pcr ~ ., proc, dist="gow",scaling=TRUE) # Model with all explanatory variables

## With scope present, the default direction is "both"
#ordistep(wgbsnull, scope = formula(wgbsdb), perm.max = 200)

#by terms
tmw <- as.data.frame(anova(wgbsdb, by="terms", permu=10000)) # test for sign. environ. variables
tmw
# Df   SumOfSqs        F Pr(>F)
# Taxon_H.C.     1 0.05988854 1.515357  0.096
# Tissue_Liver   1 0.13240726 3.350296  0.001
# Tissue_Spleen  1 0.11579391 2.929929  0.001
# Sex_Male       1 0.05365433 1.357613  0.105
# Residual       7 0.27664745       NA     NA

# adjusted R^2
R2adjw <- round(RsquareAdj(wgbsdb)$adj.r.squared,3)
pvalw <- anova(wgbsdb) # overall test of the significant of the analysis
plabw <- round(pvalw[1,4],3)
plabw
# [1] 0.001
R2adjw
# [1] 0.319

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

pdf('RDA_WGBS-Small.pdf', height=4, width=4.5)
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

pdf('RDA_WGBS-Legend.pdf', height=6, width=8)
plot.new()
legend(x = "top",inset=c(0,-.1),xpd=TRUE,pt.cex=1.5,legend=allw$Legend,pt.bg=allw$Color,pch=allw$Shape,ncol=4,cex=.8)
dev.off()  

#also perform an unconstrained nmds
set.seed(123)
nmds = metaMDS(pcr, distance = "gow", trymax = 100)
nmds
data.scores <- as.data.frame(scores(nmds))
data.scores$ID <- rownames(data.scores)
data.scores <- merge(data.scores,variables,by='ID')
data.scores$SEX_SP <- paste0(data.scores$Taxon,'_',data.scores$Sex)
data.scores$SEX_SP <- factor(data.scores$SEX_SP, levels=c('C.C._Female','C.C._Male','H.C._Female','H.C._Male'))

loci <- as.data.frame(nmds$species)

np <- 
  ggplot() +
  geom_point(data = data.scores, show.legend = FALSE,aes(x=NMDS1,
                                                         y=NMDS2, 
                                                         shape = SEX_SP,
                                                         colour = Tissue,
                                                         fill = Tissue,
                                                         size=15)) +
  scale_color_manual(values=rev(viridis(3,option='turbo')))+
  scale_fill_manual(values=rev(viridis(3,option='turbo')))+
  scale_shape_manual(values=c(1,21,2,24))+
  geom_text(aes(x=Inf,y=Inf,label=paste0('Stress: ',round(nmds$stress,3))),vjust=2.5,hjust=1,size=5)+
  theme_classic(base_size=14)
np

#perform anosim:
ord <- as.data.frame(row.names(pcr))
names(ord) <- 'ID'
vars <- left_join(data.frame(name=ord),data.scores,by="ID")
ano  <- anosim(pcr, vars$Taxon, distance = "gow", permutations = 9999)
ano
anowg <- data.frame(ano$statistic,ano$signif)
anowg$variable <- 'Taxon'
names(anowg) <- c('statistic','pval','variable')
anowg$experiment <- 'WGBS'

pdf('NMDS_WGBS.pdf', height=5, width=6)
np
dev.off()  

#and traditional pca
pcd$var <- apply(pcd, 1, var, na.rm=TRUE)
pcv <- subset(pcd,var > 0)
pcz <- pcv %>% select(-c(var)) %>% t(.)
pca <- prcomp(pcz,scale=TRUE)
var_explained <- pca$sdev^2/sum(pca$sdev^2)
vec <- data.frame(pca$x)
vec$ID <- row.names(vec)
vec <- merge(vec,variables,by='ID')
vec$SEX_SP <- paste0(vec$Taxon,'_',vec$Sex)
vec$SEX_SP <- factor(vec$SEX_SP, levels=c('C.C._Female','C.C._Male','H.C._Female','H.C._Male'))

np <- 
  ggplot() +
  geom_point(data = vec, show.legend = F,aes(x=PC1,
                                             y=PC2, 
                                             shape = SEX_SP,
                                             colour = Tissue,
                                             fill = Tissue,
                                             size=15)) +
  scale_color_manual(values=rev(viridis(3,option='turbo')))+
  scale_fill_manual(values=rev(viridis(3,option='turbo')))+
  scale_shape_manual(values=c(1,21,2,24))+
  xlab(paste0('PC1 : ',round(var_explained[1],3)*100,'% Explained'))+
  ylab(paste0('PC2 : ',round(var_explained[2],3)*100,'% Explained'))+
  theme_classic(base_size=14)
np

pdf('WGBS_PCA.pdf',height=4,width=4.5)
np
dev.off()

##### CG #####
#signal experiment 
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

pdf('RDA_CG-Small.pdf', height=4, width=4.5)
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

pdf('RDA_CG-Legend.pdf', height=6, width=8)
plot.new()
legend(x = "top",inset=c(0,-.1),xpd=TRUE,pt.cex=1.5,legend=allc$Legend,pt.bg=allc$Color,pch=allc$Shape,ncol=4,cex=.8)
dev.off()  

#also perform an unconstrained nmds
set.seed(123)
nmds = metaMDS(pcr, distance = "euc", trymax = 100)
nmds
data.scores <- as.data.frame(scores(nmds))
data.scores$ID <- rownames(data.scores)
data.scores <- merge(data.scores,variables,by='ID')
data.scores$SEX_SP <- paste0(data.scores$Taxon,'_',data.scores$Sex)
data.scores$SEX_SP <- factor(data.scores$SEX_SP, levels=c('C.C._Female','C.C._Male','H.C._Female','H.C._Male'))

loci <- as.data.frame(nmds$species)

np <- 
  ggplot() +
  geom_point(data = data.scores, show.legend = FALSE,aes(x=NMDS1,
                                                         y=NMDS2, 
                                                         shape = SEX_SP,
                                                         colour = Age,
                                                         fill = Age,
                                                         size=15)) +
  scale_color_manual(values=rev(viridis(3,option='cividis')))+
  scale_fill_manual(values=rev(viridis(3,option='cividis')))+
  scale_shape_manual(values=c(1,21,2,24))+
  geom_text(aes(x=Inf,y=Inf,label=paste0('Stress: ',round(nmds$stress,3))),vjust=2.5,hjust=1,size=5)+
  theme_classic(base_size=14)
np

pdf('NMDS_CG.pdf', height=5, width=6)
np
dev.off()

#perform anosim:
ord <- as.data.frame(row.names(pcr))
names(ord) <- 'ID'
vars <- left_join(data.frame(name=ord),data.scores,by="ID")
ano  <- anosim(pcr, vars$Taxon, distance = "euc", permutations = 9999)
anocg <- data.frame(ano$statistic,ano$signif)
anocg$variable <- 'Taxon'
names(anocg) <- c('statistic','pval','variable')
anocg$experiment <- 'CG'

#and traditional pca
pcd$var <- apply(pcd, 1, var, na.rm=TRUE)
pcv <- subset(pcd,var > 0)
pcz <- pcv %>% select(-c(var)) %>% t(.)
pca <- prcomp(pcz,scale=TRUE)
var_explained <- pca$sdev^2/sum(pca$sdev^2)
vec <- data.frame(pca$x)
vec$ID <- row.names(vec)
vec$ID <- gsub('_BL_CHK_M.hz','',vec$ID)
vec <- merge(vec,variables,by='ID')
vec$SEX_SP <- paste0(vec$Taxon,'_',vec$Sex)
vec$SEX_SP <- factor(vec$SEX_SP, levels=c('C.C._Female','C.C._Male','H.C._Female','H.C._Male'))

np <- 
  ggplot() +
  geom_point(data = vec, show.legend = F,aes(x=PC1,
                                             y=PC2, 
                                             shape = SEX_SP,
                                             colour = Age,
                                             fill = Age,
                                             size=15)) +
  scale_color_manual(values=rev(viridis(3,option='cividis')))+
  scale_fill_manual(values=rev(viridis(3,option='cividis')))+
  scale_shape_manual(values=c(1,21,2,24))+
  xlab(paste0('PC1 : ',round(var_explained[1],3)*100,'% Explained'))+
  ylab(paste0('PC2 : ',round(var_explained[2],3)*100,'% Explained'))+
  theme_classic(base_size=14)
np

pdf('CG_PCA.pdf',height=4,width=4.5)
np
dev.off()

##### HZ #####
#signal experiment 
xp <- '.hz'

#select data, round to 1%, and change 0% 5mC to 1% 5mC for CCA
pcd <- class %>% select(contains(xp)) %>% select(-contains(c('fdrs','DMP','Mean','divergence','Variance','Divergence'))) %>% na.omit()
pcdd <- t(pcd)
pcr <- round(pcdd,2)
rownames(pcr) <- gsub('_BL_CHK_M.hz','',rownames(pcr))

#grab ID
variables <- read.table('All-Metadata.txt',header=T)
variables <- subset(variables,Experiment=='HZ' & Retained==1) %>% select(-ID)
variables <- variables %>% rename('ID'='Short_ID','Year'='Capture_Year')
variables$Year <- as.factor(variables$Year)
keep <- as.data.frame(row.names(pcr))
names(keep) <- 'ID'
variables <- merge(variables,keep,by='ID')
#create distance vector
dd <- prcomp(variables %>% select(Latitude,Longitude),center=TRUE,scale=TRUE)
summary(dd)
variables$Distance <- dd$x[,1] 

#pre-process data 
rda_recipe =  # which consists of the formula (outcome ~ predictors)
  recipe(~ Temperature + Chr18_Hybrid_Index + Distance, data = variables) %>%
  # and some pre-processing steps
  step_novel(all_nominal(), -all_outcomes()) %>% # converts all nominal variables to factors and takes care of other issues related to categorical variables
  step_dummy(all_nominal(), -all_outcomes()) %>% # converts our factor column ocean_proximity into numeric binary (0 and 1) variables
  step_zv(all_numeric(), -all_outcomes()) %>% # removes any numeric variables that have zero variance
  step_scale(all_numeric(), -all_outcomes()) %>% # removes any numeric variables that have zero variance
  step_corr(all_predictors(), threshold = 0.99, method = "spearman") # will remove predictor variables that have large correlations with other predictor variables
proc = rda_recipe %>% # apply the recipe to the training data
  prep(variables) %>% juice()
proc
#collinearity
hzcol = variables %>% select(Temperature,Chr18_Hybrid_Index,Distance) %>% ggpairs()

pdf('HZ_Collinearity.pdf', height=6, width=7)
print(hzcol)
dev.off() 

#identify which distance is superior
rankindex(proc , pcr, indices = c("euc", "man", "gow", "bra", "kul"),stepacross= FALSE, method = "spearman")
# euc           man           gow           bra           kul 
# 0.0670912927  0.0311082977 -0.0655515988 -0.0005964364 -0.0165826116

#constrained ordination with step selection 
hznull <- dbrda(pcr ~ 1, proc, dist="euc",scaling=TRUE)  # Model with intercept only
hzdb = dbrda(pcr ~ ., proc, dist="euc",scaling=TRUE) # Model with all explanatory variables
#hzdb = dbrda(pcr ~ ., varhz %>% select(c(Temperature)), dist="euc",scaling=TRUE) # Model with all explanatory variables

## With scope present, the default direction is "both"
#ordistep(hznull, scope = formula(hzdb), perm.max = 200)

#by terms
tmh <- as.data.frame(anova(hzdb, by="terms", permu=10000)) # test for sign. environ. variables
tmh 
# Df  Variance         F Pr(>F)
# Temperature         1  157.0519 1.0864880  0.008
# Chr18_Hybrid_Index  1  145.8967 1.0093162  0.309
# Distance            1  143.9088 0.9955636  0.518
# Residual           15 2168.2507        NA     NA

# adjusted R^2
R2adjh <- round(RsquareAdj(hzdb)$adj.r.squared,3)
pvalh <- anova(hzdb) # overall test of the significant of the analysis
plabh <- round(pvalh[1,4],3)
plabh
# [1] 0.062
R2adjh
# [1] 0.005

#plot raw
p3 <- ggord(hzdb,variables$Taxon,ellipse=F,veccol="blue",labcol="blue",addcol="black",addsize=0.0000001,addpch=27,coord_fix=FALSE)+xlab('dbRDA1')+ylab('dbRDA2')
p3$layers[[1]] <- NULL
new_dat <- data.frame(p3$data$one,p3$data$two,variables$Temperature,variables$Species,variables$Chr18_Hybrid_Index,variables$Distance)
names(new_dat) <- gsub('.*\\.','',names(new_dat))
new_dat

# add points layer manually, redefine some aesthetics
hzp <- p3 + 
  geom_point(data = new_dat, show.legend = F,aes(x=one,
                                                 y=two, 
                                                 shape=Species,
                                                 fill = Temperature,
                                                 col=Temperature,
                                                 size=10)) +
  scale_fill_gradient(low='blue',high='yellow')+
  scale_color_gradient(low='blue',high='yellow')+
  scale_shape_manual(values=c(21,24,23))+
  geom_richtext(aes(x=Inf,y=Inf,label=paste("adj. R<sup>2</sup>: ",round(R2adjh,3))),vjust=1,hjust=1,size=5,label.colour=NA)+
  geom_richtext(aes(x=Inf,y=Inf,label=paste("p = ",plabh)),vjust=2,hjust=1,size=5,label.colour=NA)+
  theme_classic(base_size=12)+
  expand_limits(x = c(min(new_dat[,1])+min(new_dat[,1])*1.5,
                      max(new_dat[,1])+max(new_dat[,1])*1.5),
                y = c(min(new_dat[,2])+min(new_dat[,2])*0.5,
                      max(new_dat[,2])+max(new_dat[,2])*1.5))

hzp

pdf('RDA_HZ-Small-Temperature_2023SEPT28.pdf', height=4, width=4.5)
hzp
dev.off()  

pdf('RDA_HZ-Legend-Temperature.pdf', height=6, width=8)
hzp
dev.off() 

#chr18 hybrid index
# add points layer manually, redefine some aesthetics
hzp <- p3 + 
  geom_point(data = new_dat, show.legend = T,aes(x=one,
                                                 y=two, 
                                                 shape=Species,
                                                 fill = Chr18_Hybrid_Index,
                                                 col=Chr18_Hybrid_Index,
                                                 size=7)) +
  scale_fill_gradient(low='grey20',high='grey80')+
  scale_color_gradient(low='grey20',high='grey80')+
  scale_shape_manual(values=c(24,21,22))+
  geom_richtext(aes(x=Inf,y=Inf,label=paste("adj. R<sup>2</sup>: ",round(R2adjh,3))),vjust=1,hjust=1,size=5,label.colour=NA)+
  geom_richtext(aes(x=Inf,y=Inf,label=paste("p = ",plabh)),vjust=2,hjust=1,size=5,label.colour=NA)+
  theme_classic(base_size=12)+
  expand_limits(x = c(min(new_dat[,1])+min(new_dat[,1])*1.5,
                      max(new_dat[,1])+max(new_dat[,1])*1.5),
                y = c(min(new_dat[,2])+min(new_dat[,2])*0.5,
                      max(new_dat[,2])+max(new_dat[,2])*1.5))

hzp

pdf('RDA_HZ-Small-Chr18HI.pdf', height=3, width=3.25)
hzp
dev.off()  

pdf('RDA_HZ-Legend-Chr18HI.pdf', height=6, width=8)
hzp
dev.off() 

#also perform an unconstrained nmds
set.seed(123)
nmds = metaMDS(pcr, distance = "euc", trymax = 100)
nmds
#plot(nmds)
data.scores <- as.data.frame(scores(nmds))
data.scores$ID <- rownames(data.scores)
data.scores <- merge(data.scores,variables,by='ID')
loci <- as.data.frame(nmds$species)

np <- 
  ggplot() +
  geom_point(data = data.scores, show.legend = F,aes(x=NMDS1,
                                                     y=NMDS2, 
                                                     shape=Species,
                                                     fill = Temperature,
                                                     col=Temperature,
                                                     size=7)) +
  scale_fill_gradient(low='blue',high='yellow')+
  scale_color_gradient(low='blue',high='yellow')+
  scale_shape_manual(values=c(24,21,22))+
  geom_text(aes(x=Inf,y=Inf,label=paste0('Stress: ',round(nmds$stress,3))),vjust=2.5,hjust=1,size=5)+
  theme_classic(base_size=14)
np

pdf('NMDS_HZ.pdf', height=5, width=6)
np
dev.off()  

np <- 
  ggplot() +
  geom_point(data = data.scores, show.legend = F,aes(x=NMDS1,
                                                     y=NMDS2, 
                                                     shape=Species,
                                                     fill = Chr18_Hybrid_Index,
                                                     col=Chr18_Hybrid_Index,
                                                     size=7)) +
  scale_fill_gradient(low='grey20',high='grey80')+
  scale_color_gradient(low='grey20',high='grey80')+
  scale_shape_manual(values=c(24,21,22))+
  geom_text(aes(x=Inf,y=Inf,label=paste0('Stress: ',round(nmds$stress,3))),vjust=2.5,hjust=1,size=5)+
  theme_classic(base_size=14)
np

pdf('NMDS_HZ-Chr18HI.pdf', height=5, width=6)
np
dev.off()  

#and traditional pca
pcd$var <- apply(pcd, 1, var, na.rm=TRUE)
pcv <- subset(pcd,var > 0)
pcz <- pcv %>% select(-c(var)) %>% t(.)
pca <- prcomp(pcz,scale=TRUE)
var_explained <- pca$sdev^2/sum(pca$sdev^2)
vec <- data.frame(pca$x)
vec$ID <- row.names(vec)
vec$ID <- gsub('_BL_CHK_M.hz','',vec$ID)
vec <- merge(vec,variables,by='ID')

np <- 
  ggplot() +
  geom_point(data = vec, show.legend = F,aes(x=PC1,
                                             y=PC2, 
                                             shape=Species,
                                             fill = Temperature,
                                             col=Temperature,
                                             size=7)) +
  scale_fill_gradient(low='blue',high='yellow')+
  scale_color_gradient(low='blue',high='yellow')+
  scale_shape_manual(values=c(24,21,22))+
  xlab(paste0('PC1 : ',round(var_explained[1],3)*100,'% Explained'))+
  ylab(paste0('PC2 : ',round(var_explained[2],3)*100,'% Explained'))+
  theme_classic(base_size=14)
np

pdf('HZ_PCA.pdf',height=4,width=4.5)
np
dev.off()

#and for hybrid index 
np <- 
  ggplot() +
  geom_point(data = vec, show.legend = F,aes(x=PC1,
                                             y=PC2, 
                                             shape=Species,
                                             fill = Chr18_Hybrid_Index,
                                             col=Chr18_Hybrid_Index,
                                             size=7)) +
  scale_fill_gradient(low='grey20',high='grey80')+
  scale_color_gradient(low='grey20',high='grey80')+
  scale_shape_manual(values=c(24,21,22))+
  xlab(paste0('PC1 : ',round(var_explained[1],3)*100,'% Explained'))+
  ylab(paste0('PC2 : ',round(var_explained[2],3)*100,'% Explained'))+
  theme_classic(base_size=14)
np

pdf('HZ_PCA-Chr18HI.pdf',height=5,width=6)
np
dev.off()

#perform anosim:
ord <- as.data.frame(row.names(pcr))
names(ord) <- 'ID'
vars <- left_join(data.frame(name=ord),data.scores,by="ID")
ano1  <- anosim(pcr, vars$Chr18_Hybrid_Index, distance = "euc", permutations = 9999)
ano1
ano2  <- anosim(pcr, vars$Temperature, distance = "euc", permutations = 9999)
ano2
ano1hz <- data.frame(ano1$statistic,ano1$signif)
ano1hz$variable <- 'Chr18_Hybrid_Index'
names(ano1hz) <- c('statistic','pval','variable')
ano2hz <- data.frame(ano2$statistic,ano2$signif)
ano2hz$variable <- 'Temperature'
names(ano2hz) <- c('statistic','pval','variable')
anohz <- rbind(ano1hz,ano2hz)
anohz$experiment <- 'HZ'

#merge all the anosim
ano <- rbind(anowg,anocg,anohz)
ano
write.table(ano,file='ANOSIM_results.txt',quote=F,sep='\t',row.names=F)
write.table(variables,file='HZ-Metadata-Distance.txt',quote=F,sep='\t',row.names=F)
```

## Ordination: Year

```R
#### PCA, NMDS, db-RDA
.libPaths('~/mambaforge/envs/randomforest/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/crow_scratch/2021-07_HybridZone/Nov/Classification_Sept')
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

#import 5mC data
class <- read.table('Full.DSS.BP-10x.Classified-Annotated-YEAR_23AUG30.txt',header=TRUE)
class <- subset(class, Region != 'Missing' & Region != '.')

##### HZ #####
#signal experiment 
xp <- '.hz'

#select data, round to 1%, and change 0% 5mC to 1% 5mC for CCA
pcd <- class %>% select(contains(xp)) %>% select(-contains(c('fdrs','DMP','Mean','divergence','Variance','Divergence'))) %>% na.omit()
pcdd <- t(pcd)
pcr <- round(pcdd,2)
rownames(pcr) <- gsub('_BL_CHK_M.hz','',rownames(pcr))

#grab ID
variables <- read.table('All-Metadata.txt',header=T)
variables <- subset(variables,Experiment=='HZ' & Retained==1) %>% select(-ID)
variables <- variables %>% rename('ID'='Short_ID','Year'='Capture_Year')
variables$Year <- as.factor(variables$Year)
keep <- as.data.frame(row.names(pcr))
names(keep) <- 'ID'
variables <- merge(variables,keep,by='ID')
#create distance vector
dd <- prcomp(variables %>% select(Latitude,Longitude),center=TRUE,scale=TRUE)
summary(dd)
variables$Distance <- dd$x[,1] 

#pre-process data 
rda_recipe =  # which consists of the formula (outcome ~ predictors)
  recipe(~ Temperature + Chr18_Hybrid_Index + Distance + Year, data = variables) %>%
  # and some pre-processing steps
  step_novel(all_nominal(), -all_outcomes()) %>% # converts all nominal variables to factors and takes care of other issues related to categorical variables
  step_dummy(all_nominal(), -all_outcomes()) %>% # converts our factor column ocean_proximity into numeric binary (0 and 1) variables
  step_zv(all_numeric(), -all_outcomes()) %>% # removes any numeric variables that have zero variance
  step_scale(all_numeric(), -all_outcomes()) %>% # removes any numeric variables that have zero variance
  step_corr(all_predictors(), threshold = 0.99, method = "spearman") # will remove predictor variables that have large correlations with other predictor variables
proc = rda_recipe %>% # apply the recipe to the training data
  prep(variables) %>% juice()
proc
#collinearity
hzcol = variables %>% select(Temperature,Chr18_Hybrid_Index,Distance,Year) %>% ggpairs()

pdf('HZ_Collinearity_Year.pdf', height=6, width=7)
print(hzcol)
dev.off() 

#identify which distance is superior
rankindex(proc , pcr, indices = c("euc", "man", "gow", "bra", "kul"),stepacross= FALSE, method = "spearman")
# euc           man           gow           bra           kul 
# 0.0670912927  0.0311082977 -0.0655515988 -0.0005964364 -0.0165826116
# euc        man        gow        bra        kul   ###YEAR
# 0.01183512 0.02663682 0.13380936 0.02782850 0.01965000  ###YEAR

#constrained ordination with step selection 
hznull <- dbrda(pcr ~ 1, proc, dist="gow",scaling=TRUE)  # Model with intercept only
hzdb = dbrda(pcr ~ ., proc, dist="gow",scaling=TRUE) # Model with all explanatory variables
#hzdb = dbrda(pcr ~ ., varhz %>% select(c(Temperature)), dist="euc",scaling=TRUE) # Model with all explanatory variables

## With scope present, the default direction is "both"
#ordistep(hznull, scope = formula(hzdb), perm.max = 200)

#by terms
tmh <- as.data.frame(anova(hzdb, by="terms", permu=10000)) # test for sign. environ. variables
tmh 
# Df  Variance         F Pr(>F)
# Temperature         1  157.0519 1.0864880  0.008
# Chr18_Hybrid_Index  1  145.8967 1.0093162  0.309
# Distance            1  143.9088 0.9955636  0.518
# Residual           15 2168.2507        NA     NA

# Df   SumOfSqs         F Pr(>F)
# Temperature         1 0.04173676 1.1241054  0.002
# Chr18_Hybrid_Index  1 0.03641941 0.9808921  0.701
# Distance            1 0.03808451 1.0257385  0.200
# Year_X2013          1 0.03906159 1.0520545  0.118
# Year_X2014          1 0.03794267 1.0219184  0.309
# Residual           13 0.48267528        NA     NA

# adjusted R^2
R2adjh <- round(RsquareAdj(hzdb)$adj.r.squared,3)
pvalh <- anova(hzdb) # overall test of the significant of the analysis
plabh <- round(pvalh[1,4],3)
plabh
# [1] 0.062
# [1] 0.018 #YEAR
R2adjh
# [1] 0.005
# [1] 0.011 #YEAR

#plot raw
p3 <- ggord(hzdb,variables$Taxon,ellipse=F,veccol="blue",labcol="blue",addcol="black",addsize=0.0000001,addpch=27,coord_fix=FALSE)+xlab('dbRDA1')+ylab('dbRDA2')
p3$layers[[1]] <- NULL
new_dat <- data.frame(p3$data$one,p3$data$two,variables$Temperature,variables$Species,variables$Chr18_Hybrid_Index,variables$Distance,variables$Year)
names(new_dat) <- gsub('.*\\.','',names(new_dat))
new_dat

# add points layer manually, redefine some aesthetics
hzp <- p3 + 
  geom_point(data = new_dat, show.legend = T,aes(x=one,
                                                 y=two, 
                                                 shape=Species,
                                                 fill = Temperature,
                                                 col = Temperature,
                                                 size=Year)) +
  scale_fill_gradient(low='blue',high='yellow')+
  scale_color_gradient('Temperature',low='blue',high='yellow')+
  scale_shape_manual(values=c(24,21,22))+
  geom_richtext(aes(x=Inf,y=Inf,label=paste("adj. R<sup>2</sup>: ",round(R2adjh,3))),vjust=1,hjust=1,size=5,label.colour=NA)+
  geom_richtext(aes(x=Inf,y=Inf,label=paste("p = ",plabh)),vjust=2,hjust=1,size=5,label.colour=NA)+
  theme_classic(base_size=12)+
  expand_limits(x = c(min(new_dat[,1])+min(new_dat[,1])*2,
                      max(new_dat[,1])+max(new_dat[,1])*3),
                y = c(min(new_dat[,2])+min(new_dat[,2])*0.5,
                      max(new_dat[,2])+max(new_dat[,2])*0.1))

hzp

pdf('RDA_HZ-Small-Temperature-Year.pdf', height=6, width=8)
hzp
dev.off()  

```


