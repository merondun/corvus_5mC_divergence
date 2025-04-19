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
thresh <- quantile(covs$total,0.99)
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
  dmr <- callDMR(dats,minCG=3,p.threshold=0.01,minlen=50)
  dmr$var <- vari
  dmrdat <- rbind(dmrdat,dmr)
}

#Save DMR
write.table(dmrdat,file='CG_DMRs-10x.txt',quote=F,sep='\t',row.names=F)

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
