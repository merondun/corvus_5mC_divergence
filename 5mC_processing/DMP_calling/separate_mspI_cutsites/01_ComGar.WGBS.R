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
masterWGBS <- Reduce(function(x,y) merge(x = x, y = y, by=c('site'),all=TRUE),vars)

write.table(masterWGBS,file='Step1.txt',quote=F,row.names=FALSE,sep='\t')

#filter according to missingness, if there's any missing individuals, drop the site, since we only have a single replicate
covdat <- masterWGBS[rowSums(is.na(masterWGBS[grepl('^Cs_', names(masterWGBS))])) == 0, ]

#and filter for maximum coverage:
covs <- covdat[grepl('^site|^Cs_|^Ts_',names(covdat))]
covs$total <- rowSums(covs[grepl('Ts_|Cs_',names(covs))],na.rm=TRUE)
thresh <- quantile(covs$total,0.99)
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
write.table(masterWGBS,'WGBS_PCA_Input-10x.txt',quote=F,row.names=FALSE)

#now ready for model
write.table(CT ,file="WGBS_GLMM-10x.txt",row.names=FALSE,quote=F,col.names=TRUE)
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.

source('~/modules/R_Functions.R')
alldats <- read.table('WGBS_GLMM-10x.txt',header=TRUE)
alldats$TC <- alldats$Cs+alldats$Ts
mdats <- alldats[,c('site','TC','Cs','variable')]
mdats <- mdats %>% separate(site,c('chr','pos'))
mdats <- mdats[!grepl('garbage',names(mdats))]
names(mdats) <- c('chr','pos','N','X','ID')
mdats$pos <- as.integer(mdats$pos)
mdats <- mdats %>% arrange(chr,pos)
nrow(mdats)/12
rm(alldats)

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
design <- design[match(allnames, design$variable), ]
rownames(design) <- design$variable

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
  dmr <- callDMR(dats,minCG=3,p.threshold=0.01,minlen=50)
  dmr$var <- vari
  dmrdat <- rbind(dmrdat,dmr)
}

#Save DMR
write.table(dmrdat,file='WGBS_DMRs-10x.txt',quote=F,sep='\t',row.names=F)

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

# Save test-statistics and p-values for RFs
dmp <- dmpdat
source('~/modules/R_Functions.R')
dmp <- rename.col(dmp,'pos','start')
d4 <- dcast(dmp, chr + start ~ var, value.var=c('stat'))
d4$site <- paste0(d4$chr,'_',d4$start)
d4 <- as.data.frame(d4)
dp1 <- d4[!grepl('chr|start|site',names(d4))]
names(dp1) <- paste0('stat_',names(dp1))
d5 <- cbind(d4[grepl('chr|start|site',names(d4))],dp1)
write.table(d5,'WGBS.DSS-10x.stat',quote=F,row.names=FALSE,sep='\t')

#Merge with PCA % Dat
pca <- read.table('WGBS_PCA_Input-10x.txt',header=TRUE)
#add a 'divergence_' prefix
pca2 <- pca[!grepl('^Cs_|^Ts_|site|Div|Var|^D_|^S_|^H_',names(pca))]
names(pca2) <- paste0('divergence_',names(pca2))
pca3 <- cbind(pca[grepl('site|Div|Var|^D_|^S_|^H_',names(pca))],pca2)

#merge fdr and pca data..
wgbss <- merge(pca3,d3,by=c('site'))
wgbsm <- wgbss[grepl('^D_|^S_|^H_|Diver|Var|fdrs|div',names(wgbss))]
names(wgbsm) <- paste0(names(wgbsm),'.wgbs')
wgbsm1 <- cbind(wgbss[,c('site','chr','start')],wgbsm)
write.table(wgbsm1,'WGBS.DSS.BP-10x.txt',quote=F,row.names=FALSE,sep='\t')
