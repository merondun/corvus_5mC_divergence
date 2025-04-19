#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_202504')
.libPaths('~/mambaforge/envs/r/lib/R/library')
### Overlap CpGs
library(DSS)
library(dplyr)
library(tidyr)
library(reshape2)
library(matrixStats)
library(ggplot2)
library(data.table)

xp = args[1]

#grab all files
files <- list.files(path=xp, pattern="*.cov.gz", full.names=TRUE, recursive=FALSE)
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
  
  #extract name which will be our new variable. Lots of file suffix across xps ...
  name <- gsub('..*/','',samp)
  name <- gsub('.CpG_5mC.cov.gz','',name)
  name <- gsub('.SE_trimmed_bismark_bt2.bismark.cov.gz','',name)
  name <- gsub('.R1_val_1_bismark_bt2_pe.bismark.cov.gz','',name)
  
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
master <- Reduce(function(x,y) merge(x = x, y = y, by=c('site'),all=TRUE),vars)

write.table(master,file=paste0(xp,'/',xp,'_Step1.txt'),sep='\t',quote=F,row.names=FALSE)
master <- read.table(paste0(xp,'/',xp,'_Step1.txt'),header=TRUE)

if (xp == "CG") {
  covdat <- master %>% filter(rowSums(is.na(select(., starts_with("Cs_")))) <= 3)
  covdat <- covdat %>%
    mutate(
      divergence_Species   = rowMeans(select(covdat, matches("^D_")), na.rm = TRUE) - rowMeans(select(covdat, matches("^S_")), na.rm = TRUE),
      divergence_StageCHK  = rowMeans(select(covdat, matches("_CHK_")), na.rm = TRUE) - rowMeans(select(covdat, matches("_ADL_")), na.rm = TRUE),
      divergence_StageYRL  = rowMeans(select(covdat, matches("_YRL_")), na.rm = TRUE) - rowMeans(select(covdat, matches("_ADL_")), na.rm = TRUE),
      divergence_Sex       = rowMeans(select(covdat, matches("_M$")), na.rm = TRUE) - rowMeans(select(covdat, matches("_F$")), na.rm = TRUE)
    )
} else if (xp == "WGBS") {
  covdat <- master %>% filter(rowSums(is.na(select(., starts_with("Cs_")))) == 0)
  covdat <- covdat %>%
    mutate(
      divergence_Species   = rowMeans(select(covdat, matches("^D_")), na.rm = TRUE) - rowMeans(select(covdat, matches("^S_")), na.rm = TRUE),
      divergence_TissueL  = rowMeans(select(covdat, matches("_BL_")), na.rm = TRUE) - rowMeans(select(covdat, matches("_M_")), na.rm = TRUE),
      divergence_TissueM  = rowMeans(select(covdat, matches("_BL_")), na.rm = TRUE) - rowMeans(select(covdat, matches("_L_")), na.rm = TRUE),
      divergence_Sex       = rowMeans(select(covdat, matches("_M$")), na.rm = TRUE) - rowMeans(select(covdat, matches("_F$")), na.rm = TRUE)
    )
} else if (xp == "HZ") {
  covdat <- master %>% filter(rowSums(is.na(select(., starts_with("Cs_")))) <= 3)
}

covdat <- covdat %>%
  mutate(Divergence = apply(select(covdat, matches("^D_|^S_|^H_")), 1, max, na.rm = TRUE) - 
           apply(select(covdat, matches("^D_|^S_|^H_")), 1, min, na.rm = TRUE))

# Total coverage and filtering
setDT(covdat)
covdat[, total := rowSums(.SD, na.rm = TRUE), .SDcols = patterns("^Cs_|^Ts_")]
thresh <- quantile(covdat$total, 0.99)
keep <- covdat[total < thresh][, total := NULL]

# Long-form with Cs and Ts
setDT(keep)
Cs <- melt(keep, id.vars = "site", measure.vars = patterns("^Cs_"), variable.name = "variable", value.name = "Cs")
Cs$variable <- gsub("^Cs_", "", Cs$variable)
Ts <- melt(keep, id.vars = "site", measure.vars = patterns("^Ts_"), variable.name = "variable", value.name = "Ts")
Ts$variable <- gsub("^Ts_", "", Ts$variable)
CT <- merge(Cs, Ts, by = c("site", "variable"))

#short-form for PCAs, etc, still has the divergence data
fwrite(keep, paste0(xp,'_PCA_Input-10x.txt'), quote = FALSE, row.names = FALSE,col.names=TRUE)

#now ready for model
fwrite(CT, paste0(xp,'_GLMM-10x.txt'), quote = FALSE, row.names = FALSE,col.names=TRUE)

# DSS DMP + DMR Detection
alldats <- read.csv(paste0(xp,'_GLMM-10x.txt'),header=TRUE)
alldats$TC <- alldats$Cs+alldats$Ts
mdats <- alldats[,c('site','TC','Cs','variable')]
mdats <- separate(mdats,site,c('chr','pos'))
names(mdats) <- c('chr','pos','N','X','ID')
mdats$pos <- as.integer(mdats$pos)
mdats <- mdats %>% arrange(chr,pos)

allnames <- NULL
#loop through them..
for (samp in unique(mdats$ID)) {
  
  var <- mdats[grepl(samp,mdats$ID),][,1:4]
  
  cat('Saved to variable: ',samp,'\n')
  #assign it to a new variable
  assign(paste0(samp),var)
  allnames <- unique(c(allnames,samp))
  
}

# Get names of samples, convert to BSobj
vars <- mget(allnames)
BSobj = makeBSseqData(vars,c(allnames))
BSobj

# Import design
design <- read.table(paste0(xp,'_Metadata.txt'),header=TRUE)

if (xp == "CG") {
  frm <- "~Species+Stage+Sex"
} else if (xp == "WGBS") {
  frm <- "~Species+Tissue+Sex"
} else if (xp == "HZ") {
  frm <- "~Chr18_Hybrid_Index+Year"
  design$Year <- as.factor(design$Year)
}

design <- design[match(allnames, design$variable), ]
rownames(design) <- design$variable

DMLfit = DMLfit.multiFactor(BSobj, design=design, formula=as.formula(frm))
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
write.table(dmrdat,file=paste0(xp,'_DMRs-10x.txt'),quote=F,sep='\t',row.names=F)
write.table(dmpdat,file=paste0(xp,'/','Step2.txt'),quote=F,sep='\t',row.names=F)
dmpdat <- read.table(paste0(xp,'/','Step2.txt'),header=TRUE)

# rename + reshape fdrs
dmp <- dmpdat %>% rename(start = pos)
d2 <- dcast(dmp[,c('chr','start','fdrs','var')], chr + start ~ var, value.var = 'fdrs') %>%
  mutate(site = paste0(chr, '_', start)) %>% as.data.frame()
dp <- d2[!grepl('chr|start|site', names(d2))]
names(dp) <- paste0('fdrs_', names(dp))
d3 <- cbind(d2[grepl('chr|start|site', names(d2))], dp)

# stat vals
d4 <- dcast(dmp[,c('chr','start','stat','var')], chr + start ~ var, value.var = 'stat') %>%
  mutate(site = paste0(chr, '_', start)) %>% as.data.frame()
dp1 <- d4[!grepl('chr|start|site', names(d4))]
names(dp1) <- paste0('stat_', names(dp1))
d5 <- cbind(d4[grepl('chr|start|site', names(d4))], dp1)
write.table(d5, paste0(xp, '_DSS.10x.stat'), quote = FALSE, row.names = FALSE, sep = '\t')

# merge with PCA
pca <- read.csv(paste0(xp, '_PCA_Input-10x.txt'))
pca3 <- pca[!grepl('^Cs_|^Ts_', names(pca))]
hzs <- left_join(d3,pca3)
hzm <- hzs[grepl('^D_|^S_|^H_|fdrs|div|Div', names(hzs))]
names(hzm) <- paste0(names(hzm),'.',xp)
hzm1 <- cbind(hzs[,c('site','chr','start')], hzm)
write.table(hzm1, paste0(xp, '_DSS.10x.txt'), quote = FALSE, row.names = FALSE, sep = '\t')

