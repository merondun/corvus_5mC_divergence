## Obtain DSS DMPs

Please see zenodo repo for the cov.gz and DSS output files.

### ComGar.RRBS

Using the cov.gz files, call DMPs 

```r
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/crow_scratch/2021-07_HybridZone/Nov/CG')
source('~/modules/R_Functions.R')
library(DSS)
library(dplyr)
library(tidyr)
library(reshape2)
library(maditr)
library(matrixStats)
library(tidyverse)
library(viridis)
library(ggplot2)
library(ggpubr)

#grab all files
files <- list.files(path=".", pattern="*.cov.gz", full.names=TRUE, recursive=FALSE)
allnames <- NULL

#loop through them..
for (samp in files) {
  
  cat('Reading in file: ',samp,'\n')
  #read data
  ztab <- gzfile(samp,'rt');
  tab <- read.table(ztab,header=F);
  #only keep positions with >= 10x
  tab <- subset(tab, (V5 + V6) >= 10)
  tab$site <- paste0(tab$V1,'_',tab$V2)
  tab$V4 <- tab$V4/100
  
  #extract name which will be our new variable
  name <- gsub('./','',samp)
  name <- gsub('.CpG_5mC.cov.gz','',name)
  
  #grab only site / countM / countU / %5mC
  tab2 <- tab[,c(7,5,6,4)]
  names(tab2) <- c('site',paste0('Cs_',name),paste0('Ts_',name),name);
  
  cat('Saved to variable: ',name,'\n')
  #assign it to a new variable
  assign(paste0(name),tab2)
  
  #add this sample to our total list
  allnames <- unique(c(allnames,name))
}

#get names
vars <- mget(allnames)
# merge the points, keeping all
masterCG <- Reduce(function(x,y) merge(x = x, y = y, by=c('site'),all=TRUE),vars)

write.table(masterCG,file='Step1.txt',quote=F,sep='\t',row.names=FALSE)
#filter according to missingness, if there's 3 more more missing individuals, drop the site
covdat <- masterCG[rowSums(is.na(masterCG[grepl('^Cs_', names(masterCG))])) <= 3, ]

#and filter for maximum coverage:
covs <- covdat[grepl('^site|^Cs_|^Ts_',names(covdat))]
covs$total <- rowSums(covs[grepl('Ts_|Cs_',names(covs))],na.rm=TRUE)
thresh <- quantile(covs$total,0.95)
keep <- subset(covs,total < thresh)
keep <- keep[grepl('^site$',names(keep))]

#calculate variable divergence
p5mC <- masterCG[grepl('^D_|^H_|^S_',names(masterCG))]
masterCG$Species <- (rowMeans(p5mC[grepl('^D_', names(p5mC))],na.rm=TRUE) - rowMeans(p5mC[grepl('^S_', names(p5mC))],na.rm=TRUE))
masterCG$StageCHK <- (rowMeans(p5mC[grepl('_CHK_', names(p5mC))],na.rm=TRUE) - rowMeans(p5mC[grepl('_ADL_', names(p5mC))],na.rm=TRUE))
masterCG$StageYRL <- (rowMeans(p5mC[grepl('_YRL_', names(p5mC))],na.rm=TRUE) - rowMeans(p5mC[grepl('_ADL_', names(p5mC))],na.rm=TRUE))
masterCG$Sex <- (rowMeans(p5mC[grepl('_M$', names(p5mC))],na.rm=TRUE) - rowMeans(p5mC[grepl('_F$', names(p5mC))],na.rm=TRUE))

divdats <- masterCG[grepl('^site$|Species$|Stage|Sex',names(masterCG))]

#filter those
maxcov <- merge(keep,divdats,by='site')
maxcov2 <- merge(maxcov,covdat,by='site')
masterCGf <- maxcov2

#remove sites with low variance, and sites undiverged
#calculate site-level range of 5mC values, and site-level variance
masterCGf$Divergence <- rowMaxs(as.matrix(masterCGf[grepl('^S_|^D_|^H_',names(masterCGf))]),na.rm=TRUE) - rowMins(as.matrix(masterCGf[grepl('^S_|^D_|^H_',names(masterCGf))]),na.rm=TRUE)
masterCGf$Variance <- apply(masterCGf[grepl('^S_|^D_|^H_',names(masterCGf))], 1, var, na.rm=TRUE)

#assess all the filters
# merged / individual missing data / CpG number filter / maximum coverage
nrow(masterCG)
nrow(covdat)
nrow(maxcov2)

#reform back to the main frame
masterCG = masterCGf

#Our data is in a poor format for looping. To fix this we can transform it long-form
Cs <- masterCG[grepl('^site|^Cs_',names(masterCG))]
Csm <- melt(Cs,id.vars=c('site'))
Csm$variable <- gsub('Cs_','',Csm$variable)
names(Csm)[names(Csm) == 'value'] <- 'Cs'

#and thymines
Ts <- masterCG[grepl('^site|^Ts_',names(masterCG))]
Tsm <- melt(Ts,id.vars=c('site'))
Tsm$variable <- gsub('Ts_','',Tsm$variable)
names(Tsm)[names(Tsm) == 'value'] <- 'Ts'

#merge
CT <- merge(Csm,Tsm,by=c('site','variable'))

#short-form for PCAs, etc, still has the divergence data
write.table(masterCG,'CG_PCA_Input-10x.txt',quote=F,row.names=FALSE)

#now ready for model
write.table(CT ,file="CG_GLMM-10x.txt",row.names=FALSE,quote=F,col.names=TRUE)
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.

source('~/modules/R_Functions.R')
#DSS DIVERGENCE
alldats <- read.table('CG_GLMM-10x.txt',header=TRUE)
alldats$TC <- alldats$Cs+alldats$Ts
mdats <- alldats[,c('site','TC','Cs','variable')]
mdats <- mdats %>% separate(site,c('chr','pos'))
names(mdats) <- c('chr','pos','N','X','ID')
mdats$pos <- as.integer(mdats$pos)
mdats <- mdats %>% arrange(chr,pos)
nrow(mdats)/24
                   
#plot
cp = mdats %>% ggplot(aes(x=variable,y=TC,fill=variable))+
  geom_boxplot(show.legend = F)+scale_fill_manual(values=viridis(24,option='turbo'))+
  theme_classic()+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+ylab('Coverage')+xlab('Sample')
png('CG_Coverage.png',units='in',height=4,width=8,res=300)
ggarrange(cp,nrow=1)
dev.off()

allnames <- NULL
#loop through them..
for (samp in unique(mdats$ID)) {
  
  var <- mdats[grepl(samp,mdats$ID),][,1:4]
  
  cat('Saved to variable: ',samp,'\n')
  #assign it to a new variable
  assign(paste0(samp),var)
  allnames <- unique(c(allnames,samp))
  
}

#get names
vars <- mget(allnames)

#create object
BSobj = makeBSseqData(vars,c(allnames))
BSobj

#import design
design <- read.table('CG_Metadat-DSS.txt',header=TRUE)
X = model.matrix(~Species+Stage+Sex, design)
X
DMLfit = DMLfit.multiFactor(BSobj, design=design, formula=~Species+Sex+Stage)
vars <- colnames(DMLfit$X)[-1]
dmpdat <- NULL
dmrdat <- NULL
for (vari in vars) {
  dats <- DMLtest.multiFactor(DMLfit,coef=vari)
  dats$var <- vari
  dmpdat <- rbind(dmpdat,dats)
}

#Reformat nicely
dmp <- dmpdat
dmp <- rename.col(dmp,'pos','start')
d1 <- dmp[,c('chr','start','fdrs','var')]
d2 <- dcast(dmp, chr + start ~ var, value.var=c('fdrs'))
d2$site <- paste0(d2$chr,'_',d2$start)
d2 <- as.data.frame(d2)
dp <- d2[!grepl('chr|start|site',names(d2))]
names(dp) <- paste0('fdrs_',names(dp))
d3 <- cbind(d2[grepl('chr|start|site',names(d2))],dp)

#Merge with PCA % Dat
pca <- read.table('CG_PCA_Input-10x.txt',header=TRUE)
#add a 'divergence_' prefix
pca2 <- pca[!grepl('^Cs_|^Ts_|site|Div|Var|^D_|^S_|^H_',names(pca))]
names(pca2) <- paste0('divergence_',names(pca2))
pca3 <- cbind(pca[grepl('site|Div|Var|^D_|^S_|^H_',names(pca))],pca2)

#merge fdr and pca data..
cgs <- merge(pca3,d3,by=c('site'))
cgm <- cgs[grepl('^D_|^S_|^H_|Diver|Var|fdrs|div',names(cgs))]
names(cgm) <- paste0(names(cgm),'.cg')
cgm1 <- cbind(cgs[,c('site','chr','start')],cgm)
write.table(cgm1,'CG.DSS.BP-10x.txt',quote=F,row.names=FALSE,sep='\t')
```

### ComGar.WGBS

```r
source('~/modules/R_Functions.R')
library(DSS)
library(dplyr)
library(tidyr)
library(reshape2)
library(maditr)
library(matrixStats)
library(tidyverse)
library(viridis)
library(ggplot2)
library(ggpubr)

#grab all files
files <- list.files(path=".", pattern="*.cov.gz", full.names=TRUE, recursive=FALSE)
allnames <- NULL

#loop through them..
for (samp in files) {

  cat('Reading in file: ',samp,'\n')
  #read data
  ztab <- gzfile(samp,'rt');
  tab <- read.table(ztab,header=F);
  #only keep positions with >= 5x
  tab <- subset(tab, (V5 + V6) >= 5)
  tab$site <- paste0(tab$V1,'_',tab$V2)
  tab$V4 <- tab$V4/100

  #extract name which will be our new variable
  name <- gsub('./','',samp)
  name <- gsub('.CpG_5mC.cov.gz','',name)

  #grab only site / countM / countU / %5mC
  tab2 <- tab[,c(7,5,6,4)]
  names(tab2) <- c('site',paste0('Cs_',name),paste0('Ts_',name),name);

  cat('Saved to variable: ',name,'\n')
  #assign it to a new variable
  assign(paste0(name),tab2)

  #add this sample to our total list
  allnames <- unique(c(allnames,name))
}

#get names
vars <- mget(allnames)
# merge the points, keeping all
masterWGBS <- Reduce(function(x,y) merge(x = x, y = y, by=c('site'),all=TRUE),vars)

write.table(masterWGBS,file='Step1.txt',quote=F,row.names=FALSE,sep='\t')

#filter according to missingness, if there's any missing individuals, drop the site, since we only have a single replicate
covdat <- masterWGBS[rowSums(is.na(masterWGBS[grepl('^Cs_', names(masterWGBS))])) == 0, ]

#and filter for maximum coverage:
covs <- covdat[grepl('^site|^Cs_|^Ts_',names(covdat))]
covs$total <- rowSums(covs[grepl('Ts_|Cs_',names(covs))],na.rm=TRUE)
thresh <- quantile(covs$total,0.95)
keep <- subset(covs,total < thresh)
keep <- keep[grepl('^site$',names(keep))]

#calculate variable divergence
p5mC <- masterWGBS[grepl('^D_|^H_|^S_',names(masterWGBS))]
masterWGBS$Species <- (rowMeans(p5mC[grepl('^D_', names(p5mC))],na.rm=TRUE) - rowMeans(p5mC[grepl('^S_', names(p5mC))],na.rm=TRUE))
masterWGBS$TissueM <- (rowMeans(p5mC[grepl('_BL_', names(p5mC))],na.rm=TRUE) - rowMeans(p5mC[grepl('_M_', names(p5mC))],na.rm=TRUE))
masterWGBS$TissueL <- (rowMeans(p5mC[grepl('_BL_', names(p5mC))],na.rm=TRUE) - rowMeans(p5mC[grepl('_L_', names(p5mC))],na.rm=TRUE))
masterWGBS$Sex <- (rowMeans(p5mC[grepl('_M$', names(p5mC))],na.rm=TRUE) - rowMeans(p5mC[grepl('_F$', names(p5mC))],na.rm=TRUE))

divdats <- masterWGBS[grepl('^site$|Species$|Tissue|Sex',names(masterWGBS))]

#filter those
maxcov <- merge(keep,divdats,by='site')
maxcov2 <- merge(maxcov,covdat,by='site')
masterWGBSf <- maxcov2

#remove sites with low variance, and sites undiverged
#calculate site-level range of 5mC values, and site-level variance
masterWGBSf$Divergence <- rowMaxs(as.matrix(masterWGBSf[grepl('^S_|^D_|^H_',names(masterWGBSf))]),na.rm=TRUE) - rowMins(as.matrix(masterWGBSf[grepl('^S_|^D_|^H_',names(masterWGBSf))]),na.rm=TRUE)
masterWGBSf$Variance <- apply(masterWGBSf[grepl('^S_|^D_|^H_',names(masterWGBSf))], 1, var, na.rm=TRUE)

#assess all the filters
# merged / individual missing data / CpG number filter / maximum coverage
nrow(masterWGBS)
nrow(covdat)
nrow(maxcov2)

#moving forward, using the final frame
masterWGBS = masterWGBSf

#Our data is in a poor format for looping. To fix this we can transform it long-form
Cs <- masterWGBS[grepl('^site|^Cs_',names(masterWGBS))]
Csm <- melt(Cs,id.vars=c('site'))
Csm$variable <- gsub('Cs_','',Csm$variable)
names(Csm)[names(Csm) == 'value'] <- 'Cs'

#and thymines
Ts <- masterWGBS[grepl('^site|^Ts_',names(masterWGBS))]
Tsm <- melt(Ts,id.vars=c('site'))
Tsm$variable <- gsub('Ts_','',Tsm$variable)
names(Tsm)[names(Tsm) == 'value'] <- 'Ts'

#merge
CT <- merge(Csm,Tsm,by=c('site','variable'))

#short-form for PCAs, etc, still has the divergence data
write.table(masterWGBS,'WGBS_PCA_Input-5x.txt',quote=F,row.names=FALSE)

#now ready for model
write.table(CT ,file="WGBS_GLMM-5x.txt",row.names=FALSE,quote=F,col.names=TRUE)
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.

source('~/modules/R_Functions.R')
alldats <- read.table('WGBS_GLMM-5x.txt',header=TRUE)
alldats$TC <- alldats$Cs+alldats$Ts
mdats <- alldats[,c('site','TC','Cs','variable')]
mdats <- mdats %>% separate(site,c('chr','pos'))
mdats <- mdats[!grepl('garbage',names(mdats))]
names(mdats) <- c('chr','pos','N','X','ID')
mdats$pos <- as.integer(mdats$pos)
mdats <- mdats %>% arrange(chr,pos)
nrow(mdats)/12
rm(alldats)
                     
#plot
cp = mdats %>% ggplot(aes(x=variable,y=TC,fill=variable))+
  geom_boxplot(show.legend = F)+scale_fill_manual(values=viridis(24,option='turbo'))+
  theme_classic()+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+ylab('Coverage')+xlab('Sample')
png('WGBS_Coverage.png',units='in',height=4,width=8,res=300)
ggarrange(cp,nrow=1)
dev.off()
                     
allnames <- NULL
#loop through them..
for (samp in unique(mdats$ID)) {

   var <- mdats[grepl(samp,mdats$ID),][,1:4]

  cat('Saved to variable: ',samp,'\n')
  #assign it to a new variable
  assign(paste0(samp),var)
  allnames <- unique(c(allnames,samp))

}

#get names
vars <- mget(allnames)

#create object
BSobj = makeBSseqData(vars,c(allnames))
BSobj

#import design
design <- read.table('WGBS_Metadat-DSS.txt',header=TRUE)
X = model.matrix(~Species+Tissue+Sex, design)
X
DMLfit = DMLfit.multiFactor(BSobj, design=design, formula=~Species+Tissue+Sex)
vars <- colnames(DMLfit$X)[-1]
dmpdat <- NULL
dmrdat <- NULL
for (vari in vars) {
  dats <- DMLtest.multiFactor(DMLfit,coef=vari)
  dats$var <- vari
  dmpdat <- rbind(dmpdat,dats)
}

#Reformat nicely
dmp <- dmpdat
dmp <- rename.col(dmp,'pos','start')
dmp$var <- gsub('SpeciesS','Species',dmp$var)
dmp$var <- gsub('SexM','Sex',dmp$var)
d1 <- dmp[,c('chr','start','fdrs','var')]
d2 <- dcast(dmp, chr + start ~ var, value.var=c('fdrs'))
d2$site <- paste0(d2$chr,'_',d2$start)
d2 <- as.data.frame(d2)
dp <- d2[!grepl('chr|start|site',names(d2))]
names(dp) <- paste0('fdrs_',names(dp))
d3 <- cbind(d2[grepl('chr|start|site',names(d2))],dp)

#Merge with PCA % Dat
pca <- read.table('WGBS_PCA_Input-5x.txt',header=TRUE)
#add a 'divergence_' prefix
pca2 <- pca[!grepl('^Cs_|^Ts_|site|Div|Var|^D_|^S_|^H_',names(pca))]
names(pca2) <- paste0('divergence_',names(pca2))
pca3 <- cbind(pca[grepl('site|Div|Var|^D_|^S_|^H_',names(pca))],pca2)

#merge fdr and pca data..
wgbss <- merge(pca3,d3,by=c('site'))
wgbsm <- wgbss[grepl('^D_|^S_|^H_|Diver|Var|fdrs|div',names(wgbss))]
names(wgbsm) <- paste0(names(wgbsm),'.wgbs')
wgbsm1 <- cbind(wgbss[,c('site','chr','start')],wgbsm)
write.table(wgbsm1,'WGBS.DSS.BP-5x.txt',quote=F,row.names=FALSE,sep='\t')

```

### HybZon.RRBS

Whole script:

```R
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/crow_scratch/2021-07_HybridZone/Nov/HZ')
.libPaths('~/mambaforge/envs/r/lib/R/library')
source('~/modules/R_Functions.R')
library(DSS)
library(dplyr)
library(tidyr)
library(reshape2)
library(maditr)
library(matrixStats)
library(tidyverse)
library(viridis)
library(ggplot2)
library(ggpubr)

#grab all files
files <- list.files(path=".", pattern="*.cov.gz", full.names=TRUE, recursive=FALSE)
allnames <- NULL

#loop through them..
for (samp in files) {
  
  cat('Reading in file: ',samp,'\n')
  #read data
  ztab <- gzfile(samp,'rt');
  tab <- read.table(ztab,header=F);
  #only keep positions with >= 10x
  tab <- subset(tab, (V5 + V6) >= 10)
  tab$site <- paste0(tab$V1,'_',tab$V2)
  tab$V4 <- tab$V4/100
  
  #extract name which will be our new variable
  name <- gsub('./','',samp)
  name <- gsub('.CpG_5mC.cov.gz','',name)
  
  #grab only site / countM / countU / %5mC
  tab2 <- tab[,c(7,5,6,4)]
  names(tab2) <- c('site',paste0('Cs_',name),paste0('Ts_',name),name);
  
  cat('Saved to variable: ',name,'\n')
  #assign it to a new variable
  assign(paste0(name),tab2)
  
  #add this sample to our total list
  allnames <- unique(c(allnames,name))
}

#get names
vars <- mget(allnames)
# merge the points, keeping all
masterHZ <- Reduce(function(x,y) merge(x = x, y = y, by=c('site'),all=TRUE),vars)

write.table(masterHZ,file='Step1.txt',sep='\t',quote=F,row.names=FALSE)

#filter according to missingness
covdat <- masterHZ[rowSums(is.na(masterHZ[grepl('^Cs_S|^Cs_D', names(masterHZ))])) <= 1, ]
covdat2 <- covdat[rowSums(is.na(covdat[grepl('^Cs_H', names(covdat))])) <= 2, ]
nrow(covdat)
nrow(covdat2)
covdat <- covdat2

#and filter for maximum coverage:
covs <- covdat[grepl('^site|^Cs_|^Ts_',names(covdat))]
covs$total <- rowSums(covs[grepl('Ts_|Cs_',names(covs))],na.rm=TRUE)
thresh <- quantile(covs$total,0.95)
keep <- subset(covs,total < thresh)
keep <- keep[grepl('^site$',names(keep))]

#calculate variable divergence
p5mC <- masterHZ[grepl('^D_|^H_|^S_',names(masterHZ))]
masterHZ$Species <- (rowMeans(p5mC[grepl('^D_', names(p5mC))],na.rm=TRUE) - rowMeans(p5mC[grepl('^S_', names(p5mC))],na.rm=TRUE))

divdats <- masterHZ[grepl('^site$|Species$',names(masterHZ))]

#filter those
maxcov <- merge(keep,divdats,by='site')
maxcov2 <- merge(maxcov,covdat,by='site')
masterHZf <- maxcov2

#remove sites with low variance, and sites undiverged
#calculate site-level range of 5mC values, and site-level variance
masterHZf$Divergence <- rowMaxs(as.matrix(masterHZf[grepl('^S_|^D_|^H_',names(masterHZf))]),na.rm=TRUE) - rowMins(as.matrix(masterHZf[grepl('^S_|^D_|^H_',names(masterHZf))]),na.rm=TRUE)
masterHZf$Variance <- apply(masterHZf[grepl('^S_|^D_|^H_',names(masterHZf))], 1, var, na.rm=TRUE)

#assess all the filters
# merged / individual missing data / CpG number filter / maximum coverage
nrow(masterHZ)
nrow(covdat)
nrow(maxcov2)
masterHZ <- masterHZf[!grepl('^minsite|^Sites_',names(masterHZf))]

#Our data is in a poor format for looping. To fix this we can transform it long-form
Cs <- masterHZ[grepl('^site|^Cs_',names(masterHZ))]
Csm <- melt(Cs,id.vars=c('site'))
Csm$variable <- gsub('Cs_','',Csm$variable)
names(Csm)[names(Csm) == 'value'] <- 'Cs'

#and thymines
Ts <- masterHZ[grepl('^site|^Ts_',names(masterHZ))]
Tsm <- melt(Ts,id.vars=c('site'))
Tsm$variable <- gsub('Ts_','',Tsm$variable)
names(Tsm)[names(Tsm) == 'value'] <- 'Ts'

#merge
CT <- merge(Csm,Tsm,by=c('site','variable'))

#short-form for PCAs, etc, still has the divergence data
write.table(masterHZ,'HZ_PCA_Input-10x.txt',quote=F,row.names=FALSE,sep='\t')

#now ready for model
write.table(CT ,file="HZ_GLMM-10x.txt",row.names=FALSE,quote=F,col.names=TRUE,sep='\t')
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.

source('~/modules/R_Functions.R')
#DSS DIVERGENCE
alldats <- read.table('HZ_GLMM-10x.txt',header=TRUE)
alldats$TC <- alldats$Cs+alldats$Ts
mdats <- alldats[,c('site','TC','Cs','variable')]
mdats <- mdats %>% separate(site,c('chr','pos'))
names(mdats) <- c('chr','pos','N','X','ID')
mdats$pos <- as.integer(mdats$pos)
mdats <- mdats %>% arrange(chr,pos)
nrow(mdats)/19
mdats$ID <- gsub('_BL_CHK_M','',mdats$ID)
mdats <- mdats %>% mutate(variable = ifelse(grepl('D_Ko',ID) == TRUE,ID,
                                            ifelse(grepl('S_Up',ID) == TRUE,ID,
                                                   gsub('D_','H_',ID))))

#plot
cp = mdats %>% ggplot(aes(x=variable,y=N,fill=variable))+
  geom_boxplot(show.legend = F)+scale_fill_manual(values=viridis(24,option='turbo'))+
  theme_classic()+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+ylab('Coverage')+xlab('Sample')
png('HZ_Coverage.png',units='in',height=4,width=8,res=300)
ggarrange(cp,nrow=1)
dev.off()

allnames <- NULL
#loop through them..
for (samp in unique(mdats$ID)) {
  
  var <- mdats[grepl(samp,mdats$ID),][,1:4]
  
  cat('Saved to variable: ',samp,'\n')
  #assign it to a new variable
  assign(paste0(samp),var)
  allnames <- unique(c(allnames,samp))
  
}

#get names
vars <- mget(allnames)

#create object
BSobj = makeBSseqData(vars,c(allnames))
BSobj

#import design
design <- read.table('HZ-Metadata-Distance.txt',header=TRUE)
design$Temperature <- as.numeric(design$Temperature)
design$Chr18_Hybrid_Index <- as.numeric(design$Chr18_Hybrid_Index)
design$Distance <- as.numeric(design$Distance)
design$Distance <- as.numeric(design$Distance)
design$Year <- as.factor(design$Year)

DMLfit = DMLfit.multiFactor(BSobj, design=design, formula=~Temperature+Chr18_Hybrid_Index+Distance+Year)
vars <- colnames(DMLfit$X)[-1]
dmpdat <- NULL
dmrdat <- NULL
for (vari in vars) {
  dats <- DMLtest.multiFactor(DMLfit,coef=vari)
  dats$var <- vari
  dmpdat <- rbind(dmpdat,dats)
}

#years
write.table(dmpdat,file='HZ_DMRs-10x_YEAR.txt',quote=F,sep='\t',row.names=F)
dmpdat = read_tsv('HZ_DMRs-10x_YEAR.txt')
yrs = dmpdat %>% select(chr,pos,fdrs,var) %>% pivot_wider(names_from=var,values_from=fdrs)

#Reformat nicely
dmp <- dmpdat
dmp <- rename.col(dmp,'pos','start')
d1 <- dmp[,c('chr','start','fdrs','var')]
d2 <- dcast(dmp, chr + start ~ var, value.var=c('fdrs'))
d2$site <- paste0(d2$chr,'_',d2$start)
d2 <- as.data.frame(d2)
dp <- d2[!grepl('chr|start|site',names(d2))]
names(dp) <- paste0('fdrs_',names(dp))
d3 <- cbind(d2[grepl('chr|start|site',names(d2))],dp)

#Merge with PCA % Dat
pca <- read.table('HZ_PCA_Input-10x.txt',header=TRUE)
#add a 'divergence_' prefix
pca2 <- pca[!grepl('^Cs_|^Ts_|site|Div|Var|^D_|^S_|^H_',names(pca))]
names(pca2) <- paste0('divergence_',names(pca2))
pca3 <- cbind(pca[grepl('site|Div|Var|^D_|^S_|^H_',names(pca))],pca2)

#merge fdr and pca data..
hzs <- merge(pca3,d3,by=c('site'))
hzm <- hzs[grepl('^D_|^S_|^H_|Diver|Var|fdrs|div',names(hzs))]
names(hzm) <- paste0(names(hzm),'.hz')
hzm1 <- cbind(hzs[,c('site','chr','start')],hzm)
write.table(hzm1,'HZ.DSS.BP-10x_YEAR.txt',quote=F,row.names=FALSE,sep='\t')
```



