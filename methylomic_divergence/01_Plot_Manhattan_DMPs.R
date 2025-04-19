#### Manhattan Plot, chr18 zoom in, and chi squared 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_202504')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(karyoploteR)
library(viridis)
library(LICORS)
library(ggplot2)
library(scales)
library(RColorBrewer)

#set up genome
genome <- read.table('genome5.7.bed',header=FALSE) %>% filter(str_detect(V1,'scaff|MT',negate=T)) %>% arrange(desc(V2))
names(genome) <- c('chr','end'); genome$start <- 1
genome$chr <- gsub('chr','',genome$chr)
#re-order...
genome = genome %>% mutate(chord = gsub('A','',ifelse(chr == 'Z',30,chr))) %>% arrange(as.numeric(chord))
crowG <- makeGRangesFromDataFrame(genome,seqnames.field = 'chr',start.field = 'start',end.field = 'end')
#for chr 18
crow18 <- makeGRangesFromDataFrame(subset(genome,chr=='18'),seqnames.field = 'chr',start.field = 'start',end.field = 'end')

##### Import 5mC data ####
class <- read.table('Full-HZ-CG.DSS.BP-10x.Classified-Annotated-SNPless_20250410.txt',header=TRUE)

# Positions significant in both ComGar and HybZon
class %>% select(chr,start,site,Class,DMP.HZ,DMP.CG,SupportCount,Region) %>% filter(Class == 'Taxon' & SupportCount > 1)
# chr    start           site Class DMP.HZ DMP.CG SupportCount     Region
# 1 chr18 10000909 chr18_10000910 Taxon  Taxon  Taxon            2 Intergenic
# 2 chr18  9005914  chr18_9005915 Taxon  Taxon  Taxon            2 Intergenic

#Split by FDRS, rename, create some extra columns for plotting 
dats = class %>% select(chr,start,Region,Class,contains('fdrs')) %>%  
  dplyr::rename(CG.Taxon = fdrs_Species.CG, CG.Sex = fdrs_Sex.CG, CG.AgeCHK = fdrs_StageCHK.CG, CG.AgeYRL = fdrs_StageYRL.CG,                
                HZ.HybridIndex = fdrs_Chr18_Hybrid_Index.HZ) %>%   
  pivot_longer(!c(chr,start,Region,Class),names_to = 'variable', values_to = 'fdrs') %>%   
  na.omit() %>%  
  mutate(site = paste0(chr,'_',start))

#create back-up, and assign color only with 'correct' consensus classifications 
dats2 = dats %>% 
  separate(variable,into=c('experiment','variable')) %>% 
  mutate(chr = gsub('chr','',chr))
#add colors 
cols = brewer.pal(7,name = 'Paired'); show_col(cols)
cols = viridis(7,option='turbo'); show_col(cols)
cols = c('#88CCEE', '#117733', '#332288', '#DDCC77','#CC6677', '#AA4499', '#DDDDDD'); show_col(cols)
dats3 <- dats2 %>% mutate(Color1 = ifelse(grepl('HybridIndex',variable) == TRUE,'#332288',
                                          ifelse(grepl('Taxon',variable) == TRUE,'#332288',
                                                 ifelse(grepl('Sex',variable) == TRUE,'#AA4499',
                                                        ifelse(grepl('Age',variable) == TRUE,'#117733',
                                                               ifelse(grepl('Tissue',variable) == TRUE,'#DDCC77','No'))))))

#substitute HI for taxon
dats4 = dats3 %>% 
  mutate(variable = gsub('HybridIndex','Taxon',variable),
         Color = as.character(ifelse(str_detect(variable,Class) & fdrs < 0.01,Color1,
                                     ifelse(fdrs < 0.01,'orange','grey80'))),
         Shape = ifelse(Region == 'Promoter',4,
                        ifelse(Region == 'CDS',1,
                               ifelse(Region == 'Intron',3,
                                      ifelse(Region == 'Intergenic',2,
                                             ifelse(Region == 'Repeat',5,8))))))

dats5 = dats4 %>% distinct
write.table(dats5,'5mC_10x_Manhattan_Input-20250411.txt',quote=F,sep='\t',row.names=F)

dats5 = read.table('5mC_10x_Manhattan_Input-20250411.txt',header=T,sep = "\t",comment.char = '')
dats5 = dats5 %>% mutate(experiment = gsub('CG','ComGar.RRBS',experiment),
                         experiment = gsub('HZ','HybZon.RRBS',experiment))

#import gene track
genes = read.table('GCF_000738735.5_ASM73873v5_genomic.genes.bed',header=F)
names(genes) = c('chr','start','end','strand','label')
genes$chr = gsub('chr','',genes$chr)
#correct strand, and also add a midpoint for label
genes = genes %>% mutate(sstart = ifelse(strand == '+',start,end),
                         send = ifelse(strand == '+',end,start),
                         midpoint = (round(((end-start)/2),0)+start))
#rename some of the LOCS that are known
genes = genes %>% mutate(label = gsub('LOC104696000','ABCA9',label),
                         label = gsub('LOC104696005','RGS9',label),
                         label = gsub('LOC104696026','ABCA9_x2',label),
                         label = gsub('LOC120410983','RGS9_x2',label))
melgene = genes %>% filter(str_detect(label,'AXIN2') | str_detect(label,'PRKCA') | str_detect(label,'CACNG'))
othergene = genes %>% filter(!str_detect(label,'AXIN2') & !str_detect(label,'PRKCA') & !str_detect(label,'CACNG'))

##### Import genetic data ####
divs <- read.table('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_Sept/GCF_000738735.5_ASM73873v5_genomic.K28Corvus-FullPopGen.5000bp.bed',header=T)
divs <- divs %>% mutate_at(grep("Theta|FuLi|Tajima|Dxy|FST|GC",colnames(.)),funs(as.numeric))
divs$chr <- gsub('chr','',divs$chr)

gb1 = genome[seq(1, nrow(genome), 2), ][[1]]
gb2 = genome[seq(2, nrow(genome), 2), ][[1]]

#subset each raw dataset into specific colors
div_gb1 <- divs %>% filter(chr %in% gb1)
div_gb1$Color <- "black"
div_gb2 <- divs %>% filter(chr %in% gb2)
div_gb2$Color <- "grey60"
divs_use <- rbind(div_gb1,div_gb2) %>% arrange(chr,Genstart)
divs_use = divs_use %>% mutate(Shape = ifelse(Region == 'Promoter',4,
                                              ifelse(Region == 'CDS',1,
                                                     ifelse(Region == 'Intron',3,
                                                            ifelse(Region == 'Intergenic',2,
                                                                   ifelse(Region == 'Repeat',5,8))))),
                               HapD = 1-HapD)

##### Plot All Variable manhattan plots ####
subs = c('Chr18','GW')
subs = 'GW'
subs = 'Chr18'

#for taxon only
tax = dats5 %>% filter(variable == 'Taxon')
for (sub in subs){
  
  cat('Plotting : ',sub,'\n')
  if (sub == 'Chr18') {ht = 5; wd = 6} else {ht = 7; wd = 7}
  png(paste0('~/symlinks/hzone/Methylation_Analyses/figures/20250414_',sub,'-Manhattan_Taxon.png'),units='in',res=600,height=ht,width=wd,bg='transparent')
  #pdf(paste0('20250411_',sub,'-Manhattan_Taxon_LABELS_.pdf'),height=ht,width=wd)
  
  #Plot Base, add labels.plotter=NULL for no chr names
  pp <- getDefaultPlotParams(plot.type=4)
  if (sub == 'Chr18') { 
    pp$leftmargin <- 0.08
    cs = 8.07e6 - 250000
    ce = 10.07e6 + 250000
    zoom18 = data.frame('18',cs,ce)
    kp <- plotKaryotype(plot.type=4, genome = crow18,plot.params = pp,zoom = zoom18,labels.plotter=NULL)
    #kp <- plotKaryotype(plot.type=4, genome = crow18,plot.params = pp,labels.plotter=NULL) #for entire chr 
    kpAddBaseNumbers(kp,tick.dist = 500000,minor.tick.dist = 100000,cex=0.65)
    tracks = 5
  } else { 
    pp$leftmargin <- 0.08
    kp <- plotKaryotype(plot.type=4, genome = crowG,plot.params = pp,labels.plotter=NULL)
    #for staggered chr names
    evens <- ((1:length(genome$chr))%%2)==0
    chr.names <- kp$chromosomes
    even.names <- chr.names
    even.names[!evens] <- ""
    odd.names <- chr.names
    odd.names[evens] <- ""
    kpAddChromosomeNames(kp, chr.names = odd.names,cex=0.7)
    kpAddChromosomeNames(kp, chr.names = even.names, yoffset = -7,cex=0.7)
    kpAddBaseNumbers(kp,tick.dist = 25000000,cex=0.5)
    tracks = 5
  }
  
  run <- 'tax'
  vars <- c('ComGar.RRBS','HybZon.RRBS')
  counter <- 0
  
  for (varz in vars) {
    counter <- counter+1
    cat('Current Autotrack: ',counter,' for ',varz,'\n')
    at <- autotrack(current.track = counter, total.tracks = tracks);
    at$r1 <- at$r1-0.05
    #grab data 
    data <- get(run)
    if (sub == 'Chr18') { data <- subset(data,chr=='18' & start > cs-5000 & start < ce+5000) } else {}
    dat = data %>% filter(experiment == varz)
    dat = dat %>% mutate(FDR = -log10(fdrs))
    dt = dat   
    
    ## If you want to apply a threshold to pvalues, e.g. because huge outliers, use this 
    #if (sub == 'Chr18') { thresh = 0.999 } else { thresh = 0.9999}
    #maxd = quantile(dt$FDR,thresh)
    #dt = dt %>% mutate(FDR = LICORS::threshold(FDR,min=0,max=maxd)) %>% dplyr::rename(val = FDR)
    
    ## Otherwise, just use this 
    maxd <- max(dt$FDR)
    dt <- dt %>% dplyr::rename(val = FDR)
    datNS = subset(dt,Color == 'orange' | Color == 'grey80')
    datS = subset(dt,Color != 'grey80' & Color != 'orange')
    
    #add pvalues, first non-significant
    kpPoints(kp, chr = datNS$chr, x=datNS$start, y=datNS$val,
             #col=transparent(as.character(datNS$Color1),amount=0.85),
             col=transparent('#DDDDDD',amount=0.5),
             pch=datNS$Shape,cex=0.7,
             ymin=0,ymax=maxd,r0=at$r0, r1=at$r1);
    #then significant
    kpPoints(kp, chr = datS$chr, x=datS$start, y=datS$val,
             col=as.character(datS$Color),pch=datS$Shape,cex=0.9,
             ymin=0,ymax=maxd,r0=at$r0, r1=at$r1);
    kpAxis(kp, ymin =0, ymax=maxd,cex=0.7,r0=at$r0, r1=at$r1,col="black",numticks = 2);
    kpAbline(kp, h=0,r0=at$r0, r1=at$r1, col="black", lwd=1,lty=1)
    kpAbline(kp, h=-log10(0.01),ymin=0,ymax=maxd,r0=at$r0, r1=at$r1, col="black", lwd=1,lty=2)
  }
  
  if (sub == 'Chr18') {
    #add genes
    counter = counter+1
    at <- autotrack(current.track = counter, total.tracks = tracks);
    kpArrows(kp, chr=genes$chr, x0=genes$sstart, x1=genes$send, y0=0.15, y1=0.15,length=0.1,lwd=1,r0 = at$r0,r1=at$r1)
    kpText(kp, chr=melgene$chr, x=melgene$midpoint, labels = melgene$label, col='darkmagenta',y=0.6,srt=45,cex=.5,r0 = at$r0,r1=at$r1)
    kpText(kp, chr=othergene$chr, x=othergene$midpoint, labels = othergene$label, y=0.6,srt=45,cex=.3,r0 = at$r0,r1=at$r1)
  } else {}
  
  #add genetic data
  genvars = c('FST','Dxy')
  for (var2 in genvars) {
    counter <- counter+1
    cat('Current Autotrack: ',counter,' for ',var2,'\n')
    at <- autotrack(current.track = counter, total.tracks = tracks);
    at$r1 <- at$r1-0.03
    data <- divs_use
    if (sub == 'Chr18') { data <- subset(data,chr=='18' & start > cs-5000 & start < ce+5000) } else {}
    dat <- data %>% select(c(chr,Region,Genstart,Genend,Color,var2,Shape))
    names(dat) = c('chr','Region','start','end','Color','variable','Shape')
    if (var2 == 'Dxy' & sub == 'Chr18') { maxd = 0.01} else if (var2 == 'Dxy' & sub == 'GW') { maxd = 0.02} else if (var2 == 'FST') { maxd = 1} else {maxd = max(dat$variable,na.rm=TRUE)}
    dat = dat %>% arrange(chr,start)
    # #add genetic values
    kpPoints(kp,chr=dat$chr,x=dat$start,y=dat$variable,ymin=min(dat$variable,na.rm=TRUE),
            col=dat$Color,ymax=maxd,r0=at$r0, r1=at$r1,clipping=TRUE)
    # grey = dat %>% filter(Color == 'grey60')
    # black = dat %>% filter(Color == 'grey10')
    # kpLines(kp,chr=grey$chr,x=grey$start,y=grey$variable,ymin=min(dat$variable,na.rm=TRUE),
    #         col='grey10',ymax=maxd,r0=at$r0, r1=at$r1,clipping=TRUE)
    # kpLines(kp,chr=black$chr,x=black$start,y=black$variable,ymin=min(dat$variable,na.rm=TRUE),
    #         col='grey60',ymax=maxd,r0=at$r0, r1=at$r1,clipping=TRUE)
    kpAxis(kp, ymin =min(dat$variable,na.rm=TRUE), ymax=maxd,cex=0.7,r0=at$r0, r1=at$r1,col="black",numticks = 2)
  }
  
  
  #plot chromosome 18 boundary
  if (sub == 'Chr18') {
    kpRect(kp,chr='18',border='black',x0=cs+150000,x1=ce-150000,y0=0,lty=2,y1=1,r0=0,r1=1) 
  } else { 
    kpRect(kp,chr='18',border=transparent('black',amount=0.4),col=NA,x0=0,x1=12661267,y0=0,y1=1,r0=0,r1=at$r1,lty=2) 
  }
  
  dev.off()
  
}

# Plot other region
c='4'
cs=2.5e6
ce=4e6
png(paste0('~/symlinks/hzone/Methylation_Analyses/figures/20250411_',c,'-',cs,'-',ce,'-Manhattan_Taxon.png'),units='in',res=600,height=ht,width=wd,bg='transparent')
kp <- plotKaryotype(plot.type=4, genome = crowG,plot.params = pp,zoom = data.frame(c,cs,ce),labels.plotter=NULL)
kpAddBaseNumbers(kp,tick.dist = 500000,minor.tick.dist = 100000,cex=0.65)
run <- 'tax'
vars <- c('ComGar.RRBS','HybZon.RRBS')
counter <- 0
tracks=4

for (varz in vars) {
  counter <- counter+1
  cat('Current Autotrack: ',counter,' for ',varz,'\n')
  at <- autotrack(current.track = counter, total.tracks = tracks);
  at$r1 <- at$r1-0.05
  #grab data
  data <- get(run)
  dat = data %>% filter(experiment == varz & chr == c & start > cs & start < ce)
  dat = dat %>% mutate(FDR = -log10(fdrs))
  dt = dat   
  
  ## Otherwise, just use this 
  maxd <- max(dt$FDR)
  dt <- dt %>% dplyr::rename(val = FDR)
  datNS = subset(dt,Color == 'orange' | Color == 'grey80')
  datS = subset(dt,Color != 'grey80' & Color != 'orange')
  
  #add pvalues, first non-significant
  kpPoints(kp, chr = datNS$chr, x=datNS$start, y=datNS$val,
           col=transparent('#DDDDDD',amount=0.5),
           pch=datNS$Shape,cex=0.7,
           ymin=0,ymax=maxd,r0=at$r0, r1=at$r1);
  #then significant
  kpPoints(kp, chr = datS$chr, x=datS$start, y=datS$val,
           col=as.character(datS$Color),pch=datS$Shape,cex=0.9,
           ymin=0,ymax=maxd,r0=at$r0, r1=at$r1);
  kpAxis(kp, ymin =0, ymax=maxd,cex=0.7,r0=at$r0, r1=at$r1,col="black",numticks = 2);
  kpAbline(kp, h=0,r0=at$r0, r1=at$r1, col="black", lwd=1,lty=1)
  kpAbline(kp, h=-log10(0.01),ymin=0,ymax=maxd,r0=at$r0, r1=at$r1, col="black", lwd=1,lty=2)
}

#add genetic data
genvars = c('FST','Dxy')
for (var2 in genvars) {
  counter <- counter+1
  cat('Current Autotrack: ',counter,' for ',var2,'\n')
  at <- autotrack(current.track = counter, total.tracks = tracks);
  at$r1 <- at$r1-0.03
  data <- divs_use
  dat <- data %>% select(c(chr,Region,Genstart,Genend,Color,var2,Shape))
  names(dat) = c('chr','Region','start','end','Color','variable','Shape')
  if (var2 == 'Dxy' & sub == 'Chr18') { maxd = 0.01} else if (var2 == 'Dxy' & sub == 'GW') { maxd = 0.02} else if (var2 == 'FST') { maxd = 1} else {maxd = max(dat$variable,na.rm=TRUE)}
  dat = dat %>% arrange(chr,start)
  # #add genetic values
  grey = dat %>% filter(Color == 'grey60')
  black = dat %>% filter(Color == 'grey10')
  kpLines(kp,chr=grey$chr,x=grey$start,y=grey$variable,ymin=min(dat$variable,na.rm=TRUE),
          col='grey10',ymax=maxd,r0=at$r0, r1=at$r1,clipping=TRUE)
  kpLines(kp,chr=black$chr,x=black$start,y=black$variable,ymin=min(dat$variable,na.rm=TRUE),
          col='grey60',ymax=maxd,r0=at$r0, r1=at$r1,clipping=TRUE)
  kpAxis(kp, ymin =min(dat$variable,na.rm=TRUE), ymax=maxd,cex=0.7,r0=at$r0, r1=at$r1,col="black",numticks = 2)
}
dev.off()

