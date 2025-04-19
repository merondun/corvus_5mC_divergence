#### dbRDA
.libPaths('~/mambaforge/envs/randomforest/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_202504/')
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
class <- read.table('Full-HZ-CG.DSS.BP-10x.Classified-Annotated-SNPless_20250410.txt',header=TRUE)
class <- subset(class, Region != 'Missing' & Region != '.')

##### CG #####
# signal experiment 
xp <- '.CG'

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
# 0.4387298 0.4184629 0.4076826 0.4024229 0.4072586 

#constrained ordination with step selection 
cgnull <- dbrda(pcr ~ 1, proc, dist="euc",scaling=TRUE)  # Model with intercept only
cgdb = dbrda(pcr ~ ., proc, dist="euc",scaling=TRUE) # Model with all explanatory variables
cgdb = dbrda(pcr ~ ., proc %>% select(-Taxon_H.C.), dist="euc",scaling=TRUE) # Model with all explanatory variables

## With scope present, the default direction is "both"
#ordistep(cgnull, scope = formula(cgdb), perm.max = 200)

#by terms
tmc <- as.data.frame(anova(cgdb, by="terms", permu=10000)) # test for sign. environ. variables
tmc
# Df   Variance         F Pr(>F)
# Taxon_H.C.  1  158.08589 2.0093595  0.001
# Age_CHK     1  211.91170 2.6935154  0.001
# Age_YRL     1   56.42415 0.7171823  0.987
# Sex_Male    1  174.29972 2.2154463  0.001
# Residual   19 1494.82056        NA     NA

# adjusted R^2
R2adjc <- round(RsquareAdj(cgdb)$adj.r.squared,3)
pvalc <- anova(cgdb) # overall test of the significant of the analysis
plabc <- round(pvalc[1,4],3)
plabc
# [1] 0.001
R2adjc
# [1] 0.136

#### WITHOUT TAXON
# plabc
# [1] 0.001
# R2adjc
# [1] 0.093


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

pdf('~/symlinks/hzone/Methylation_Analyses/figures/20250415_RDA_CG-Small.pdf', height=4, width=4.5)
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

pdf('~/symlinks/hzone/Methylation_Analyses/figures/20250415__RDA_CG-Legend.pdf', height=6, width=8)
plot.new()
legend(x = "top",inset=c(0,-.1),xpd=TRUE,pt.cex=1.5,legend=allc$Legend,pt.bg=allc$Color,pch=allc$Shape,ncol=4,cex=.8)
dev.off()  

##### HZ #####
# signal experiment 
xp <- '.HZ'

#select data, round to 1%, and change 0% 5mC to 1% 5mC for CCA
pcd <- class %>% select(contains(xp)) %>% select(-contains(c('fdrs','DMP','Mean','divergence','Variance','Divergence'))) %>% na.omit()
pcdd <- t(pcd)
pcr <- round(pcdd,2)
rownames(pcr) <- gsub('.HZ','',rownames(pcr))

#grab ID
variables <- read.table('HZ_Metadata.txt',header=T)
variables <- variables %>% select(-ID) %>% dplyr::rename(ID = variable)
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
# 0.1258837 0.1608645 0.2069293 0.1883800 0.1900351 

#constrained ordination with step selection 
hznull <- dbrda(pcr ~ 1, proc, dist="euc",scaling=TRUE)  # Model with intercept only
hzdb = dbrda(pcr ~ ., proc, dist="euc",scaling=TRUE) # Model with all explanatory variables
hzdb = dbrda(pcr ~ ., proc %>% select(-Chr18_Hybrid_Index), dist="euc",scaling=TRUE) # Model with all explanatory variables

## With scope present, the default direction is "both"
#ordistep(hznull, scope = formula(hzdb), perm.max = 200)

#by terms
tmh <- as.data.frame(anova(hzdb, by="terms", permu=10000)) # test for sign. environ. variables
tmh 
# Df  Variance        F Pr(>F)
# Chr18_Hybrid_Index  1  174.3159 1.188712  0.038
# Year_X2013          1  221.2921 1.509056  0.001
# Year_X2014          1  146.7541 1.000760  0.369
# Residual           18 2639.5683       NA     NA

# adjusted R^2
R2adjh <- round(RsquareAdj(hzdb)$adj.r.squared,3)
pvalh <- anova(hzdb) # overall test of the significant of the analysis
plabh <- round(pvalh[1,4],3)
plabh
# [1] 0.003
R2adjh
# [1] 0.032


#### WITHOUT HYBRID INDEX 
# > plabh
# [1] 0.002
# > # [1] 0.003
#   > R2adjh
# [1] 0.027

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
  scale_shape_manual(values=c(24,21,23))+
  geom_richtext(aes(x=Inf,y=Inf,label=paste("adj. R<sup>2</sup>: ",round(R2adjh,3))),vjust=1,hjust=1,size=5,label.colour=NA)+
  geom_richtext(aes(x=Inf,y=Inf,label=paste("p = ",plabh)),vjust=2,hjust=1,size=5,label.colour=NA)+
  theme_classic(base_size=12)+
  expand_limits(x = c(min(new_dat[,1])+min(new_dat[,1])*1.5,
                      max(new_dat[,1])+max(new_dat[,1])*1.5),
                y = c(min(new_dat[,2])+min(new_dat[,2])*0.5,
                      max(new_dat[,2])+max(new_dat[,2])*1.5))+
  guides(fill=guide_legend(override.aes=list(shape=21)))

hzp

pdf('~/symlinks/hzone/Methylation_Analyses/figures/20250415_RDA_HZ.pdf', height=4, width=4.5)
hzp
dev.off()  

pdf('~/symlinks/hzone/Methylation_Analyses/figures/20250415_RDA_HZ-Legend.pdf', height=6, width=8)
hzp
dev.off() 

