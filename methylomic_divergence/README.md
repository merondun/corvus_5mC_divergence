# Manhattan Plots

## Taxon 

Manhattan plot, showing the FDR-corrected pvalues for 5mC and pop gen data:

```R
#### Manhattan Plot, chr18 zoom in, and chi squared 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_Sept')
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
class <- read.table('Full.DSS.BP-10x.Classified-Annotated_20250107.txt',header=TRUE)

#Split by FDRS, rename, create some extra columns for plotting 
dats = class %>% select(chr,start,Region,Class,contains('fdrs')) %>%  
  dplyr::rename(WGBS.Taxon = fdrs_Species.wgbs, WGBS.Sex = fdrs_Sex.wgbs, WGBS.TissueL = fdrs_TissueL.wgbs, WGBS.TissueM = fdrs_TissueM.wgbs,                
                CG.Taxon = fdrs_Species.cg, CG.Sex = fdrs_Sex.cg, CG.AgeCHK = fdrs_StageCHK.cg, CG.AgeYRL = fdrs_StageYRL.cg,                
                HZ.HybridIndex = fdrs_Chr18_Hybrid_Index.hz) %>%   
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
write.table(dats5,'5mC_10x_Manhattan_Input-20250107.txt',quote=F,sep='\t',row.names=F)

dats3 = read.table('5mC_10x_Manhattan_Input-20250107.txt',header=T,sep = "\t",comment.char = '')
dats3 = dats3 %>% mutate(experiment = gsub('CG','ComGar.RRBS',experiment),
                         experiment = gsub('WGBS','ComGar.WGBS',experiment),
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
div_gb1$Color <- "grey60"
div_gb2 <- divs %>% filter(chr %in% gb2)
div_gb2$Color <- "grey10"
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

#for taxon only
tax = dats3 %>% filter(variable == 'Taxon')
for (sub in subs){
  
  cat('Plotting : ',sub,'\n')
  if (sub == 'Chr18') {ht = 5; wd = 6} else {ht = 7; wd = 7}
  png(paste0('20250107_',sub,'-Manhattan_Taxon.png'),units='in',res=600,height=ht,width=wd,bg='transparent')
  #pdf(paste0('20250107_',sub,'-Manhattan_Taxon_LABELS_.pdf'),height=ht,width=wd)
  
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
    tracks = 6
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
    tracks = 6
  }
  
  run <- 'tax'
  vars <- c('ComGar.WGBS','ComGar.RRBS','HybZon.RRBS')
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
    kpText(kp, chr=melgene$chr, x=melgene$midpoint, labels = melgene$label, col='darkmagenta',y=0.6,srt=45,cex=.35,r0 = at$r0,r1=at$r1)
    kpText(kp, chr=othergene$chr, x=othergene$midpoint, labels = othergene$label, y=0.6,srt=45,cex=.35,r0 = at$r0,r1=at$r1)
  } else {}
  
  #add genetic data
  if (gens == 'DxyFST') { genvars = c('FST','Dxy') } else { genvars = c('FST','Dxy','FuLiD','TajimaD','HapD') }
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
    grey = dat %>% filter(Color == 'grey60')
    black = dat %>% filter(Color == 'grey10')
    kpLines(kp,chr=grey$chr,x=grey$start,y=grey$variable,ymin=min(dat$variable,na.rm=TRUE),
            col='grey10',ymax=maxd,r0=at$r0, r1=at$r1,clipping=TRUE)
    kpLines(kp,chr=black$chr,x=black$start,y=black$variable,ymin=min(dat$variable,na.rm=TRUE),
            col='grey60',ymax=maxd,r0=at$r0, r1=at$r1,clipping=TRUE)
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

#chr 7 region
dats3 %>% filter( chr == '7' & start > 33600035 & start < 33956161)
```

### Chr7 Inspections

```R
.libPaths('~/mambaforge/envs/r/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/crow_scratch/2021-07_HybridZone/Nov/Classification_Sept')
library(tidyverse)
options(scipen=999)
all = read_tsv('Full.DSS.BP-10x.Classified-Annotated_23FEB09.txt')

#chr 7 region
all %>% filter( chr == 'chr7' & start > 33600035 & start < 33956161) %>% select(chr,start,Class,DMP.hz,DMP.cg,DMP.wgbs,FST,Dxy) %>% 
  ggplot(aes(y=Dxy))

chr7 = all %>% filter(chr == 'chr7')
chr7 = chr7 %>% mutate(region = ifelse(start > 33600000 & start < 33957000, 'Region','Background'))
targs = chr7 %>% filter(start > 33e6 & start < 34e6) %>% select(chr,start,Class,Region,DMP.hz,DMP.cg,DMP.wgbs) %>% filter(DMP.cg == 'Taxon') %>% data.frame
write.table(targs %>% mutate(end=start) %>% select(chr,start,end),file='Chr7_Targets.bed',quote=F,sep='\t',row.names=F,col.names=F)
chr7 %>% ggplot(aes(x=region,y=Dxy,fill=region))+
  geom_boxplot()+
  theme_bw()

genes = read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/crow_scratch/2021-07_HybridZone/Nov/Classification_Sept/tmp.2',col_names = F)
names(genes) = c('chr','start','end','gene','score','strand')
genes %>% filter( chr == 'chr7' & start >= 33600035 & start <= 33956161) 

chr7 %>% filter(start >= 33943990 & start <= 33951404 ) %>% select(chr,start,Class,DMP.hz,DMP.cg,DMP.wgbs) %>% data.frame

```

## All Variables

Plot FDR-p-values for all variables in 5mC experiments: 

```R
#### Manhattan Plot for all variables, chr18 zoom in, and chi squared 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_Sept')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(karyoploteR)
library(viridis)
library(LICORS)
library(ggplot2)
library(scales)
library(RColorBrewer)

# Set up genome
genome <- read.table('genome5.7.bed',header=FALSE) %>% filter(str_detect(V1,'scaff|MT',negate=T)) %>% arrange(desc(V2))
names(genome) <- c('chr','end'); genome$start <- 1
genome$chr <- gsub('chr','',genome$chr)
#re-order...
genome = genome %>% mutate(chord = gsub('A','',ifelse(chr == 'Z',30,chr))) %>% arrange(as.numeric(chord))
crowG <- makeGRangesFromDataFrame(genome,seqnames.field = 'chr',start.field = 'start',end.field = 'end')
#for chr 18
crow18 <- makeGRangesFromDataFrame(subset(genome,chr=='18'),seqnames.field = 'chr',start.field = 'start',end.field = 'end')

dats3 = read.table('5mC_10x_Manhattan_Input-20250107.txt',header=T,sep = "\t",comment.char = '')
dats3 = dats3 %>% mutate(experiment = gsub('CG','ComGar.RRBS',experiment),
                         experiment = gsub('WGBS','ComGar.WGBS',experiment),
                         experiment = gsub('HZ','HybZon.RRBS',experiment))

##### Plot All Variable manhattan plots ####
subs = c('Chr18','GW')
subs = 'GW'
subs = 'Chr18'

# For All Variables 
for (sub in subs){
  
  cat('Plotting : ',sub,'\n')
  if (sub == 'Chr18') {ht = 5; wd = 6} else {ht = 7; wd = 7}
  png(paste0('20250107_',sub,'-Manhattan_AllVariables.png'),units='in',res=600,height=ht,width=wd,bg='transparent')

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
    tracks = 9
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
    tracks = 9
  }
  
  run <- 'dats3'
  vars <- c('ComGar.WGBS','ComGar.RRBS','HybZon.RRBS')
  counter <- 0
  
  for (varz in vars) {

    #grab data for each experimnt
    data <- get(run)
    dat = data %>% filter(experiment == varz)
    dat = dat %>% mutate(FDR = -log10(fdrs))
    
    # Grab data for each variable 
    for (cov in unique(dat$variable)) {
      counter <- counter+1
      cat('Current Autotrack: ',counter,' for ',varz,' and variable: ',cov,'\n')
      at <- autotrack(current.track = counter, total.tracks = tracks);
      at$r1 <- at$r1-0.02
      dt = dat %>% filter(variable == cov)
    
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
      kpAxis(kp, ymin =0, ymax=maxd,cex=0.5,r0=at$r0, r1=at$r1,col="black",numticks = 2);
      kpAbline(kp, h=0,r0=at$r0, r1=at$r1, col="black", lwd=1,lty=1)
      kpAbline(kp, h=-log10(0.01),ymin=0,ymax=maxd,r0=at$r0, r1=at$r1, col="black", lwd=1,lty=2)
      kpAddLabels(kp,cex=0.5,labels =cov,r0=at$r0+0.02, r1=at$r1,col="black",srt=0,label.margin = 0.02)
    }
    
  }
  
  #plot chromosome 18 boundary
  if (sub == 'Chr18') {
    kpRect(kp,chr='18',border='black',x0=cs+150000,x1=ce-150000,y0=0,lty=2,y1=1,r0=0,r1=1) 
  } else { 
    kpRect(kp,chr='18',border=transparent('black',amount=0.4),col=NA,x0=0,x1=12661267,y0=0,y1=1,r0=0,r1=at$r1,lty=2) 
  }
  
  dev.off()
  
}

```

## Volcano Plots

```R
#### Volcano plot of FDR ~ 5mC Divergence 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_Sept')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(karyoploteR)
library(viridis)
library(LICORS)
library(ggplot2)
library(scales)
library(RColorBrewer)

class = read.table('Full.DSS.BP-10x.Classified-Annotated_20250107.txt',header=T,sep = "\t",comment.char = '')#volcano
v1 = class %>% select(site,Class,contains(c('divergence'),ignore.case = FALSE))  %>% pivot_longer(!c(site,Class),values_to = 'Divergence') %>% na.omit %>% mutate(name = gsub('divergence_','',name))
v1 = v1 %>% separate(name,into=c('Variable','Experiment'))
v2 = class %>% select(site,Class,contains(c('fdrs'),ignore.case = FALSE))  %>% pivot_longer(!c(site,Class),values_to = 'FDRS') %>% na.omit %>% mutate(name = gsub('fdrs_','',name))
v2 = v2 %>% separate(name,into=c('Variable','Experiment'))
v3 = left_join(v1,v2)
v4 <- v3 %>% filter(!grepl('Hybrid|hz',Experiment)) %>% na.omit
v4 = v4 %>% mutate(Experiment = gsub('^cg','ComGar.RRBS',Experiment),
                   Experiment = gsub('^wgbs','ComGar.WGBS',Experiment),
                   Variable = gsub('Species','Taxon',Variable),
                   Variable = gsub('StageCHK','StageChick',Variable),
                   Variable = gsub('StageYRL','StageYearling',Variable),
                   Variable = gsub('TissueL','TissueLiver',Variable),
                   Variable = gsub('TissueM','TissueSpleen',Variable))

# Checkpoint
write.table(v4,file='20250107_Volcano_Input.txt',quote=F,sep='\t',row.names=F)
v4 = read.table('20250107_Volcano_Input.txt',header=T,sep = "\t",comment.char = '')

# Plot
png('20250107_VolcanoPlots-Known.png',units='in',height=5,width=7,res=600)
vp = v4 %>% 
  filter(Class != 'Unknown') %>% 
  ggplot(aes(y=-log10(FDRS),x=Divergence,col=Variable))+
  geom_point()+
  geom_hline(yintercept=-log10(0.01),lty=2,col='blue')+
  scale_color_viridis(discrete=TRUE,option='turbo')+
  facet_wrap(Variable+Experiment~.,scales='free')+
  theme_bw()
vp
dev.off()

```



## Hybrid Index to 5mC

```R
#### Plot 5mC ~ Hybrid Index 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_Sept')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(ggplot2)
library(scales)
library(RColorBrewer)

# Import 5mC data
class <- read.table('Full.DSS.BP-10x.Classified-Annotated_20250107.txt',header=TRUE)

#Plot Sites of Interest 
md = read.table('All-Metadata.txt',header=T)
hz = class %>% filter(Class == 'Taxon' & SupportCount >= 2 & DMP.hz == 'Taxon') %>% select(site,contains('.hz'))

#first hybrid index
hi = hz %>% 
  select(-c(DMP.hz,Divergence.hz,Variance.hz,contains('divergence'),contains('fdrs'))) %>%
  pivot_longer(!site,names_to = 'Short_ID',values_to = 'Methylation') %>% 
  mutate(Short_ID = gsub('.hz','',Short_ID))
him = merge(hi,md %>% select(Short_ID,Chr18_Hybrid_Index,Species),by='Short_ID')
him = him %>% mutate(HI = ifelse(Species == 'C.corone','C.c. corone',
                                 ifelse(Species == 'C.cornix','C.c. cornix',
                                        ifelse(Chr18_Hybrid_Index < 0.25,'C.c. corone-like',
                                               ifelse(Chr18_Hybrid_Index > 0.75,'C.c. cornix-like','Hybrid')))))

him$HI = factor(him$HI,levels=c('C.c. corone','C.c. corone-like','Hybrid','C.c. cornix-like','C.c. cornix'))

hip = him %>% dplyr::rename(Taxon = Species) %>% 
  ggplot(aes(x=Chr18_Hybrid_Index,y=Methylation,fill=HI,shape=HI))+
  geom_jitter(width=0.01,size=3,show.legend = F)+
  scale_fill_grey(start=0,end=1)+
  scale_shape_manual(values=c(24,25,23,22,21))+ylim(c(-.01,1.01))+
  xlab('Genetic Hybrid Index')+ylab('Methylation (%)')+
  facet_wrap(site~.,scales='free')+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
hip

ggsave('20250107_Hybrid5mC_Distributions.pdf',hip,dpi=300,height=2,width=6)

```

### Variation Across Hybrid Groups

```bash
.libPaths('~/miniconda3/envs/randomforest/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_Sept/')
library(tidyverse)
library(data.table)
library(viridis)

##### Import 5mC data ####
class <- read.table('Full.DSS.BP-10x.Classified-Annotated_20250107.txt',header=TRUE)

#Grab important columns
hz = class %>% dplyr::select(chr,start,Class,Focal_Region,contains(c('.hz'))) %>% dplyr::select(-contains(c('fdrs','DMP','Div','Var','div')))

#only keep sites which are taxon/tissue/conserved 
hzm = hz %>% filter(grepl('Taxon|Tissue',Class)) 
hzm =  hzm[rowSums(is.na(hzm[grepl('^S_|^H_|^D_', names(hzm))])) <= 3, ] #remove sites with 3 or more individuals missing data 

#count sites
hzm %>% dplyr::count(Class,Focal_Region)
hzm2 = hzm %>% 
  pivot_longer(!c(chr,start,Class,Focal_Region))
hzm2 = hzm2 %>% mutate(name = gsub('.hz','',name)) %>% na.omit
mdat = read.table('All-Metadata.txt',header=TRUE)
hzm3 = left_join(hzm2,mdat %>% select(Short_ID,Chr18_Hybrid_Index,Species) %>% dplyr::rename(name = Short_ID)) %>% na.omit
hzm3 = hzm3 %>% mutate(HI = ifelse(Species == 'C.corone','C.c. corone',
                                   ifelse(Species == 'C.cornix','C.c. cornix',
                                          ifelse(Chr18_Hybrid_Index < 0.25,'C.c. corone-like',
                                                 ifelse(Chr18_Hybrid_Index > 0.75,'C.c. cornix-like','Hybrid')))))

fullboots = NULL
options(dplyr.summarise.inform = FALSE)
set.seed(123)
B <- 999 # number of bootstrap replicates

for (grp in unique(hzm3$HI)) {
  for (vr in unique(hzm3$Class)) {
    cat('Working on Experiment & Variable: ',grp,' ',vr,'\n')
    bd = hzm3 %>% filter(HI == grp & Class == vr)
    n = bd %>% filter(Focal_Region == 'Focal') %>% nrow
    bdat=NULL
    ## Run the bootstrap procedure:
    for (b in 1:B) {
      cat('Iteration: ',b,'\n')
      d1 = bd %>% group_by(Focal_Region) %>% slice_sample(n=n,replace=TRUE) %>% summarize(mean = mean(value,na.rm=TRUE))
      bdat = rbind(bdat,d1 )
    }
    bdat$HI = grp
    bdat$Variable = vr
    names(bdat) = c('Focal_Region','Boots','HI','Variable')
    fullboots = rbind(fullboots,bdat)
  }
}

write.table(fullboots,file='HypoBoots_bg_20250107.txt',quote=F,sep='\t',row.names=F)
bdatplot = read.table('HypoBoots_bg_20250107.txt',header=TRUE,sep='\t')

#Plot 
bdatplot$HI = factor(bdatplot$HI,levels=c('C.c. corone','C.c. corone-like','Hybrid','C.c. cornix-like','C.c. cornix'))
#calculate CIs
cis = bdatplot %>% 
  group_by(HI,Variable,Focal_Region) %>% 
  summarize(LowerBound = quantile(Boots,0.025), #can change based on your significance threshold 
            UpperBound = quantile(Boots,0.975))
cis = cis %>% group_by(HI,Variable) %>% 
  mutate(Signif = ifelse(LowerBound[Focal_Region == 'Focal'] > UpperBound[Focal_Region == 'Undiverged'],'*',
                         ifelse(UpperBound[Focal_Region == 'Focal'] < LowerBound[Focal_Region == 'Undiverged'],'*','n.s.')))
bp = bdatplot %>% 
  ggplot(aes(x=HI,fill=Focal_Region,y=Boots))+
  geom_boxplot(alpha=0.5)+
  #geom_point(data=pvals,aes(x=HI,group=HI,y=Observed),pch=21,size=5,fill='darkviolet',position=position_dodge(width=1))+
  geom_text(data=cis,aes(x=HI,group=HI,y=Inf,label=Signif,size=Signif),vjust=1.5,position=position_dodge(width=1))+
  scale_size_manual(values=c(5,4))+
  theme_bw(base_size=8)+
  geom_hline(yintercept=0,lty=2)+
  scale_fill_manual(values=c('darkorchid2','grey60'))+
  ylim(c(0,1))+
  facet_grid(Variable ~., scales='free')+
  xlab('')+ylab('Relative Methylation Status')
bp

pdf('20250107_Methylation_levels_Hybrids_Focal_Bootstraps.pdf',height=3,width=5)
bp
dev.off()


#HI, plot raw 
hzm3$HI = factor(hzm3$HI,levels=c('C.c. corone','C.c. corone-like','Hybrid','C.c. cornix-like','C.c. cornix'))
samples <- hzm3 %>% select(name,Chr18_Hybrid_Index) %>% distinct %>% arrange(Chr18_Hybrid_Index)
hzm3$name <- factor(hzm3$name,levels=samples$name)

mcp = hzm3 %>% 
  filter(Class == 'Taxon' & Focal_Region == 'Focal') %>% 
  ggplot(aes(x=name,y=value,fill=HI))+
  geom_violin(alpha=0.5)+
  geom_jitter()+
  facet_grid(Class ~ HI,scales='free')+
  scale_fill_grey()+xlab('')+ylab('Methylation Level (%) In Focal Region CpGs')+
  theme_bw(base_size=8)+
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
mcp
pdf('20250107_5mC_DistributionsRaw.pdf',height=2.5,width=6)
mcp 
dev.off()

```


