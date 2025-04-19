setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/HZ_Full')
.libPaths('~/mambaforge/envs/r/lib/R/library')
### HZ Extract
library(DSS)
library(dplyr)
library(tidyr)
library(reshape2)
library(matrixStats)
library(ggplot2)

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
  name <- gsub('_BL_CHK_M.CpG_5mC.cov.gz','',name)
  
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
masterHZ <- read.table('Step1.txt',header=TRUE)

#filter according to missingness
covdat <- masterHZ
covdat2 <- masterHZ[rowSums(is.na(covdat[grepl('^Cs_H', names(covdat))])) <= 2, ]
nrow(covdat)
nrow(covdat2)
covdat <- covdat2

#and filter for maximum coverage:
# Add total coverage per site
covdat$total <- rowSums(covdat[grep("^Cs_|^Ts_", names(covdat))], na.rm = TRUE)
thresh <- quantile(covdat$total, 0.99)
keep <- covdat %>% filter(total < thresh) %>% select(-total)

#remove sites with low variance, and sites undiverged
#calculate site-level range of 5mC values, and site-level variance
mcdat <- keep[grepl('^H_',names(keep))] 
keep$Divergence <- rowMaxs(as.matrix(mcdat),na.rm=TRUE, useNames = FALSE) - rowMins(as.matrix(mcdat),na.rm=TRUE, ,useNames = FALSE)
keep$Variance <- rowVars(as.matrix(mcdat), na.rm = TRUE, useNames = FALSE)

#assess all the filters
# merged / individual missing data / CpG number filter / maximum coverage
nrow(masterHZ)
nrow(covdat)
nrow(keep)

#Our data is in a poor format for looping. To fix this we can transform it long-form
.libPaths('~/mambaforge/envs/R/lib/R/library')
library(data.table)
setDT(keep)
Cs <- melt(keep, id.vars = "site", measure.vars = patterns("^Cs_"), variable.name = "variable", value.name = "Cs")
Cs$variable <- gsub("^Cs_", "", Cs$variable)
Ts <- melt(keep, id.vars = "site", measure.vars = patterns("^Ts_"), variable.name = "variable", value.name = "Ts")
Ts$variable <- gsub("^Ts_", "", Ts$variable)
CT <- merge(Cs, Ts, by = c("site", "variable"))

#short-form for PCAs, etc, still has the divergence data
fwrite(keep, 'HZ_PCA_Input-10x.txt', quote = FALSE, row.names = FALSE,col.names=TRUE)

#now ready for model
fwrite(CT, 'HZ_GLMM-10x.txt', quote = FALSE, row.names = FALSE,col.names=TRUE)

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.

#DSS DIVERGENCE
alldats <- read.csv('HZ_GLMM-Parentals-10x.txt',header=TRUE)
alldats$TC <- alldats$Cs+alldats$Ts
mdats <- alldats[,c('site','TC','Cs','variable')]
mdats <- mdats %>% separate(site,c('chr','pos'))
names(mdats) <- c('chr','pos','N','X','ID')
mdats$pos <- as.integer(mdats$pos)
mdats <- mdats %>% arrange(chr,pos)
mdats$ID <- gsub('_BL_CHK_M','',mdats$ID)

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
design <- read.table('HZ_Metadat-Parentals.txt',header=TRUE)
design$Year <- as.factor(design$Year)
#design$Temperature <- as.numeric(design$Temperature)
design$Chr18_Hybrid_Index <- as.numeric(design$Chr18_Hybrid_Index)
design <- design[match(allnames, design$variable), ]
rownames(design) <- design$variable

DMLfit = DMLfit.multiFactor(BSobj, design=design, formula=~Chr18_Hybrid_Index+Year)
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
write.table(dmrdat,file='HZ_DMRs-Parentals-YearChr18HI_10x.txt',quote=F,sep='\t',row.names=F)

#Reformat nicely
dmp <- dmpdat
source('~/modules/R_Functions.R')
dmp <- rename.col(dmp,'pos','start')
d1 <- dmp[,c('chr','start','fdrs','pvals','stat','var')]
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
write.table(d5,'HZ.DSS.BP-Parentals-YearChr18HI_10x.stat',quote=F,row.names=FALSE,sep='\t')

#Merge with PCA % Dat
pca <- read.csv('HZ_PCA_Input-Parentals-10x.txt',header=TRUE)
#add a 'divergence_' prefix
pca2 <- pca[!grepl('^Cs_|^Ts_',names(pca))]

#merge fdr and pca data..
hzs <- merge(pca2,d3,by=c('site'))
hzm <- hzs[grepl('^D_|^S_|^H_|Diver|Var|fdrs|div',names(hzs))]
names(hzm) <- paste0(names(hzm),'.hz')
hzm1 <- cbind(hzs[,c('site','chr','start')],hzm)
write.table(hzm1,'HZ.DSS.BP-Parentals-YearChr18HI_10x.txt',quote=F,row.names=FALSE,sep='\t')

# Sanity checks
hzm1 %>% filter(chr == 'chr18' & start %in% c(10000910, 9048799, 9317485,9317486)) 

