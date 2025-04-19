# **DNA methylation reflects cis-genetic differentiation across the European crow hybrid zone**

```
title: "DNA methylation reflects cis-genetic differentiation across the European crow hybrid zone"
author: "Justin Merondun"
date: "`r Sys.Date()`"
output:
  html_document:
    theme: journal
    toc: true
    toc_float:
      collapsed: true
```

This study comprises 58 unique methylation libraries spanning a common garden experiment (ComGar) for RRBS and WGBS data, a hybrid zone experiment (HybZon), and 28 illumina resequencing libraries. Sample details can be found in Supplementary Table 1. 

## File Paths

Paths to relevant files.

Using the newest crow genome, with RefSeq annotation: [here](https://www.ncbi.nlm.nih.gov/assembly/GCF_000738735.5).

```bash
#genome-based files
genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Genomes/Corvus_ASM73873v5
genome=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.fa
cgi=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.CGI.bed
repeats=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.REPEATS.bed
cutsites=${genomedir}/Corvus_cornix__S_Up_H32_v5.6_polished.CUT-SITES.bed
chrsbed=${genomedir}/Corvus_cornix__S_Up_H32_v5.6_polished.CHRS.bed
chrs=${genomedir}/Corvus_cornix__S_Up_H32_v5.6_polished.genome
```

I first renamed chromosomes from the NC_ format of NCBI to chr1 with this:

```bash
library(phylotools)
#Just have a .list file with the new chromosomes **IN THE SAME ORDER AS THE ORIGINAL FILE** 
old_name <- get.fasta.name('GCF_000738735.5_ASM73873v5_genomic.fa')
new_name <- read.table('newnames.list',header=FALSE)
ref2 <- data.frame(old_name,new_name)
rename.fasta(infile = "GCF_000738735.5_ASM73873v5_genomic.fa", ref_table = ref2, outfile = "GCF_000738735.5_ASM73873v5_genomic.fasta")
```

# Extract Genome Features

Create tracks for analyses: CpG island, genomic context track, cut-site track, and repeats. 

## CpG Islands

For our analyses, let's start by making a CpG island track, so we can at least look at CpG islands on top of our gff gene track. I'll do this with HMMs using [makeCGI](http://www.haowulab.org/software/makeCGI/index.html)

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=48:00:00

SCRATCH=/tmp/$SLURM_JOB_ID

Rscript makeCGI.R

```

Rscript:

```r
library(makeCGI);
.CGIoptions=CGIoptions();
.CGIoptions$rawdat.type="txt";
.CGIoptions$species="GCF_000738735.5_ASM73873v5_genomic";
.CGIoptions;
makeCGI(.CGIoptions);
5000
```

Now output a bed file following strict CpG island definitions, the standard definition is:

* length at least 200
* Percent GC > 50%
* obs/expected CG ratio > 60%

```bash
awk '($4 > 200) && ($7 > 0.5) && ($8 > 0.6){print $1, $2, $3}' CGI-GCF_000738735.5_ASM73873v5_genomic.txt | sed '1d' | tr ' ' \\t  > ../../GCF_000738735.5_ASM73873v5_genomic.CGI.bed
```

Which gives us about 30K CpG islands:

```bash
cat GCF_000738735.5_ASM73873v5_genomic.CGI.bed | wc -l
30459
```

Then, grab the promoter regions (bearing in mind strand specificity) using this script, using the .gff and .fa files: [script source]([GitHub - RimGubaev/extract_promoters: Bash script for promoter sequences extraction](https://github.com/RimGubaev/extract_promoters)) 

```bash
bedtools intersect -a GCF_000738735.5_ASM73873v5_genomic.CGI.bed -b GCF_000738735.5_ASM73873v5_genomic.promoters.bed -wa | bedtools sort -i - |  bedtools merge -i - > GCF_000738735.5_ASM73873v5_genomic.CGI-Promoters.bed
```

## Repeat Track

Using repeatmasker, I will subset our 5mC into repeat and non-repeat. 

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=48:00:00

SCRATCH=/tmp/$SLURM_JOB_ID

#genome-based files
genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Genomes/Corvus_ASM73873v5
genome=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.fa

RepeatMasker ${genome} -species chicken -nolow -pa 10 -gff -dir ${genomedir};
```

And we can use the raw gff file output for masking, or convert it to a standard bed with the following:

```bash
grep -v '#' GCF_000738735.5_ASM73873v5_genomic.fa.out.gff | awk '{OFS="\t"}{print $1, $4, $5, $7}' > GCF_000738735.5_ASM73873v5_genomic.REPEATS.bed
```

## Annotation Feature Track

Now, to get an intergenic, CDS, intron, promoter, or repeat annotation. Start here:

```bash
#grab gene coordinates
awk '{OFS="\t"}{print $1, $4, $5, $3, $7}' GCF_000738735.5_ASM73873v5_genomic.CHR.gff | grep -v '#' | egrep 'gene' | bedtools sort -i - | bedtools merge -i - > GENES.bed
#grab cds
awk '{OFS="\t"}{print $1, $4, $5, $3, $7}' GCF_000738735.5_ASM73873v5_genomic.CHR.gff | grep -v '#' | egrep 'CDS' | bedtools sort -i - | bedtools merge > CDS.bed

chrs=GCF_000738735.5_ASM73873v5_genomic.CHRS.bed
repeats=GCF_000738735.5_ASM73873v5_genomic.REPEATS.bed
genes=GCF_000738735.5_ASM73873v5_genomic.GENES.bed
promoters=GCF_000738735.5_ASM73873v5_genomic.CGI-Promoters.bed
cds=GCF_000738735.5_ASM73873v5_genomic.CDS.bed
intron=GCF_000738735.5_ASM73873v5_genomic.INTRONS.bed

#subtract overlapping, starting with largest. Therefore, priority given to promoters > repeats > cds > introns > intergenic
#intergenic
bedtools subtract -a $chrs -b $intron | bedtools subtract -a - -b $cds | bedtools subtract -a - -b $repeats | bedtools subtract -a - -b $promoters | grep -v 'scaff' | grep -v 'chrMT' > genomic_featuresCDS/INTERGENIC.bed
#introns
bedtools subtract -a $intron -b $cds | bedtools subtract -a - -b $repeats | bedtools subtract -a - -b $promoters | grep -v 'scaff' | grep -v 'chrMT' > genomic_featuresCDS/INTRON.bed
#cds
bedtools subtract -a $cds -b $repeats | bedtools subtract -a - -b $promoters | grep -v 'scaff' | grep -v 'chrMT' > genomic_featuresCDS/CDS.bed
#repeats
bedtools subtract -a $repeats -b $promoters | grep -v 'scaff' | grep -v 'chrMT' > genomic_featuresCDS/REPEATS.bed
#promoters
grep -v 'scaff' $promoters | grep -v 'chrMT' > genomic_featuresCDS/PROMOTER.bed

#add identifier
awk '{OFS="\t"}{print $0,"Intergenic"}' genomic_featuresCDS/INTERGENIC.bed > genomic_featuresCDS/INTERGENIC.txt
awk '{OFS="\t"}{print $0,"Intron"}' genomic_featuresCDS/INTRON.bed > genomic_featuresCDS/INTRON.txt
awk '{OFS="\t"}{print $0,"CDS"}' genomic_featuresCDS/CDS.bed > genomic_featuresCDS/CDS.txt
awk '{OFS="\t"}{print $1, $2, $3, "Repeat"}' genomic_featuresCDS/REPEATS.bed > genomic_featuresCDS/REPEATS.txt
awk '{OFS="\t"}{print $0,"Promoter"}' genomic_featuresCDS/PROMOTER.bed > genomic_featuresCDS/PROMOTER.txt

#sort
cat genomic_featuresCDS/*.txt | bedtools sort -i - > genomic_featuresCDS/features.txt
bedtools merge -i genomic_featuresCDS/features.txt -d -1 -c 4 -o distinct > genomic_featuresCDS/distinct.txt

#we have overlapping ranges, subtract 1 from the end to ensure overlap won't duplicate
awk '{OFS="\t"}{print $1, $2, $3-1, $4}' distinct.txt > Sub.bed

#now remove any zero length
awk '{OFS="\t"}{print $1, $2, $3, $4, $3-$2}' Sub.bed | awk '$5 > 0' | awk '{OFS="\t"}{print $1, $2, $3, $4}'  > finfeatures.bed

#now check that there are NO commas in the distinct.txt, which indicates we have overlapping conflicts
grep ',' finfeatures.bed
#should be empty

#final file:
cp finfeatures.bed ../GCF_000738735.5_ASM73873v5_genomic.GenomicFeaturesCDS.bed
```

Plot:

```r
setwd('~/files/class')
source('~/R_Functions.R')
library(dplyr)
library(karyoploteR)
library(tidyr)
library(maditr)
library(RColorBrewer)
library(viridis)
library(yarrr)
library(LICORS)
library(ggplot2)
library(scales)

features <- read.table('GCF_000738735.5_ASM73873v5_genomic.GenomicFeaturesCDS.bed',header=FALSE)
features$length <- features$V3-features$V2
names(features) <- c('chr','start','end','feature','length')
col1 <- as.data.frame(unique(features$feature))
col1$Color <- viridis(5)
names(col1) <- c('feature','Color')
feat <- merge(features,col1,by='feature')

#Set up genome
genome <- read.table('GCF_000738735.5_ASM73873v5_genomic.genome',header=FALSE)
genome$start <- 1
genome <- genome[!grepl('scaff|MT',genome$V1),]
genome$lab <- gsub('chr','',genome$V1)
genome <- genome %>% arrange(as.integer(lab))
crowG <- makeGRangesFromDataFrame(genome,seqnames.field = 'V1',start.field = 'start',end.field = 'V2')

#plot
png('Features.png',height=8,width=15,units='in',res=600)
kp <- plotKaryotype(plot.type=4, genome = crowG[grepl('chr18',crowG@seqnames)],labels.plotter=NULL)
#kpAddBaseNumbers(kp,tick.dist = 25000000, tick.len = 4, tick.col="gray", cex=0.5,)
#kpAddChromosomeNames(kp,srt=0, cex=.75,yoffset = -1)

kpRect(kp,chr=feat$chr,x0=feat$start,x1=feat$end,y0=0,y1=1,border=NA,col=feat$Color,r0=0.05,r1=0.25)

#legend
legend(x = "top",inset=c(0,-.1),xpd=TRUE,pch=15,legend=col1$feature,col=col1$Color,ncol=4,cex=0.6)
dev.off()

sum <- feat %>% group_by(feature) %>% summarize(Total.BP=sum(length),
                                                MeanLength.BP=mean(length),
                                                MinLength.BP=min(length),
                                                MaxLength.BP=max(length),
                                                Count = n())
write.table(sum,'Feature_Counts.txt',quote=F,sep='\t',row.names=FALSE)
```

## GC-Content

Calculate in 500bp :

```bash
seqkit sliding -s 500 -W 500 GCF_000738735.5_ASM73873v5_genomic.fa | seqkit fx2tab -n -g | tr ':' '\t' | tr '-' '\t' | sed 's/_sliding//g' > GCF_000738735.5_ASM73873v5_genomic.GC-Content.500bp.bed
```

## Cut-site Track

Let's find all reasonable cut-sites.

simRAD script:

```r
#simRAD
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SimRAD")

library(SimRAD)

#this will only import the chromosome-specific contig 
crow <- ref.DNAseq("SEQ.fasta",subselect.contigs=FALSE)

#I only have a single-cutter enzyme... 
#Define the restriction enzyme 1 recognition pattern:
#MspI :  C'CGG
cs_5p1 <- "C"
cs_3p1 <- "CGG"

#Simulation of the digestion just like before, but by adding the new recognition site
crow.dig <- insilico.digest(crow, cs_5p1, cs_3p1, verbose=T)

#will need to change type="AA" if you have 2 enzymes 
#selecting fragment with ends corresponding to each enzyme: since engle end, only enzyme A to enzyme A
crow.sel <- adapt.select(crow.dig, type="AA", cs_5p1, cs_3p1)

#modify insert size 
#selecting the fraction of fragments between 40 and 350bp:
nar.crow <- size.select(crow.sel,  min.size = 40, max.size = 350, graph=F, verbose=T)

#this will have the output of the explicit genomic ranges 
write.table(nar.crow@ranges,file="RANGES.txt")
```

Now, loop through each chromosome and run the above script, then merge all the intervals. The genome bed file simply looks like this:

```bash
head Corvus_cornix__S_Up_H32_v5.6_polished.CHRS.bed
chr2    0       154820475
chr1A   0       74117537
chr27   0       4995057
chr18   0       12661267
chr10   0       20564847
chr28   0       5328479
chr9    0       25711917
```

Loop: 

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=48:00:00

SCRATCH=/tmp/$SLURM_JOB_ID

#genome-based files
genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Genomes/Corvus_ASM73873v5
genome=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.fa

#you need R ... and SimRAD installed 

for i in $(awk '{print $1}' GCF_000738735.5_ASM73873v5_genomic.CHRS.bed); do

#simRAD can only hold on ~250 mb scaffolds at once, subselect scaffold
grep -w ${i} GCF_000738735.5_ASM73873v5_genomic.CHRS.bed > GRAB.bed
#use seqtk to grab this scaffold from the genome 
seqtk subseq $genome GRAB.bed > SEQ.fasta

#run simRAD on this scaffold
Rscript simRAD.R

#reformat the output, so that we have $chr $start $end 
awk -v c=${i} '{OFS="\t"}{print c, $3, $4}' RANGES.txt | sed '1d' > ${i}.ranges.bed

done

#merge all the chromosomes back together 
cat *.ranges.bed | bedtools sort -i - > GCF_000738735.5_ASM73873v5_genomic.CUT-SITES.bed
```

This gives us 306088 fragments, of which 299587 are on named chromosomes. 

We can count the 'CpG' motifs with:

```bash
seqkit subseq --bed GCF_000738735.5_ASM73873v5_genomic.CUT-SITES.bed $genome > frag_seqs.fa
grep -v "^>" frag_seqs.fa  | tr -d "\n" | tr -c "[GCgc]" "\n" | grep -v '^$' | grep -o -i 'CG' | wc -l
```

Which leaves us with 1,601,786 'CG' motifs, meaning we could have a maximum of 3,203,572 methylated/non-methylated positions, since each CpG contains 2 potential C positions. 

This is compared to the entire genome, which contains 9,764,344 CG motifs, leaving 19,528,688 potential positions. 

# Resequencing Reanalysis

## Pre-processing

Take our fastqs (139 libraries for 30 samples), remove adapters, and map them to the reference. Note -- one was identified as female after mapping, so the lowest coverage C. cornix was also dropped, bringing us to n=28 illumina resequencing. Mean cov: 15.13x, range (8.03-30.85). 

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=48:00:00

SCRATCH=/tmp/$SLURM_JOB_ID

genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Genomes/Corvus_ASM73873v5
genome=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.fa
trimdir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Crow_Population_Resequencing_ASM73873v5/trimmed
bamdir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Crow_Population_Resequencing_ASM73873v5/bams

#positional arguments
RUN=$1
SAMPLE=$2
R1=$3
R2=$4

#adapter trim
bbduk.sh t=8 -Xmx30g in=${R1} in2=${R2} out=$trimdir/${RUN}.trim.fastq.gz minlen=25 qtrim=rl trimq=2 ktrim=r k=23 mink=11 ref=~/modules/bbmap/adapters.fa hdist=1 tpe tbo

#assign read groups
ID=${RUN}
SM=${SAMPLE}
LB=${SM}
PU=${RUN}

#alignment
bwa mem -M -p -t 10 -R "@RG\tID:${ID}\tSM:${SM}\tPL:ILLUMINA\tLB:${LB}\tPU:${PU}" ${genome} $trimdir/${RUN}.trim.fastq.gz | samtools sort -@10 -o $SCRATCH/${RUN}.raw.bam -;

samtools view -b -f 1 -F 1536 $SCRATCH/${RUN}.raw.bam > ${bamdir}/${RUN}.bam
samtools index -b ${bamdir}/${RUN}.bam
```

And then merge the bams:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=48:00:00

SCRATCH=/tmp/$SLURM_JOB_ID

genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Genomes/Corvus_ASM73873v5
genome=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.fa
crowdir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Crow_Population_Resequencing_ASM73873v5
trimdir=$crowdir/trimmed
bamdir=$crowdir/bams

#positional argument at sample level
RUN=$1

mkdir ${crowdir}/stats/alignment
mkdir ${crowdir}/sample_bams
mkdir ${crowdir}/stats/coverage

#first, add all the library files into a file that ends with .lib
samtools merge - ${bamdir}/${RUN}*.bam | samtools sort -@10 -o $crowdir/merged_bams/${RUN}.bam
samtools index $crowdir/merged_bams/${RUN}.bam

#Mark AND remove duplicates (just to be safe in case ANGSD is used later)
gatk MarkDuplicatesSpark -I $crowdir/merged_bams/${RUN}.bam -M ${crowdir}/stats/alignment/${RUN}.duplicates.txt -O $SCRATCH/${RUN}.dedup.bam --create-output-bam-index true --remove-all-duplicates true

#Remove reads over scaffold end
gatk CleanSam -I $SCRATCH/${RUN}.dedup.bam -O ${crowdir}/sample_bams/${RUN}.bam -R ${genome} --CREATE_INDEX true

#Summary Stats & Coverage
gatk CollectAlignmentSummaryMetrics -R ${genome} -I ${crowdir}/sample_bams/${RUN}.bam -O ${crowdir}/stats/alignment/${RUN}.alignments.txt --METRIC_ACCUMULATION_LEVEL=SAMPLE --METRIC_ACCUMULATION_LEVEL=READ_GROUP

#calculate coverage
mosdepth --threads 10 --use-median --by 100000 --fast-mode --no-per-base ${crowdir}/stats/coverage/${RUN} ${crowdir}/sample_bams/${RUN}.bam
```

And submit at sample level:

```bash
for sample in $(cat Samples.list); do sbatch -J ${sample} Merge_Clean.sh ${sample}; done 
```

**Examine coverage for resequencing data:**

Take the mosdepth outputs, merge em:

```bash
for i in $(ls *.bed.gz| sed 's/\..*//g'); do zcat ${i}*bed.gz | awk -v s=${i} '{OFS="\t"}{print $0, s}' > ${i}.bed; done
```

And in R:

```R
#### Manhattan Plot, chr18 zoom in, and chi squared 
setwd('F:/Research/scratch/crow_hybrid_paper/resequencing/2022_07/')
library(dplyr)
library(karyoploteR)
library(tidyverse)
library(viridis)
library(LICORS)
library(ggplot2)
library(scales)
library(reshape2)

#set up genome
genome <- read.table('F:/Research/scratch/crow_hybrid_paper/genome5.7.bed',header=FALSE) %>% filter(str_detect(V1,'scaff|MT',negate=T)) %>% arrange(desc(V2))
names(genome) <- c('chr','end'); genome$start <- 1
genome$chr <- gsub('chr','',genome$chr)
crowG <- makeGRangesFromDataFrame(genome,seqnames.field = 'chr',start.field = 'start',end.field = 'end')
#for chr 18
crow18 <- makeGRangesFromDataFrame(subset(genome,chr=='18'),seqnames.field = 'chr',start.field = 'start',end.field = 'end')

#import
c <- read.table('Coverage-Resequencing.bed',header=F)
names(c) <- c('chr','start','end','Coverage','ID')

m <- read.table('F:/Research/scratch/crow_hybrid_paper/All-Metadata.txt',header=T)
m <- subset(m,Experiment=='WGS')
cm <- merge(c,m,by='ID')
cm <- cm %>% mutate(Discrete = cut(Coverage, breaks = c(0,5,20,5000), right = F, labels = c('< 5x','5 - 20x','> 20x')),
                    DiscreteColor = cut(Coverage, breaks = c(0,5,20,5000), right = F, labels = c('black','yellow','red')))
cm$chr <- gsub('chr','',cm$chr)
leg <- cm %>% select(Discrete,DiscreteColor) %>% unique()
leg$Discrete <- factor(leg$Discrete,levels=c('< 5x','5 - 20x','> 20x'))

##### Plot Coverage WGS ####
png('WGS_Coverage.png',units='in',res=600,height=7,width=10,bg='white')

#Plot Base, add labels.plotter=NULL for no chr names
pp <- getDefaultPlotParams(plot.type=4)
pp$leftmargin <- 0.08
kp <- plotKaryotype(plot.type=4, genome = crowG,plot.params = pp,labels.plotter=NULL)

legend(x = "top",inset=c(0,-.1),xpd=TRUE,pch=22,legend=levels(leg$Discrete),col=as.character(levels(leg$DiscreteColor)),fill=as.character(levels(leg$DiscreteColor)),ncol=3,cex=0.7)

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

run <- 'cm'
vars <- unique(cm$ID)
tracks <- length(vars)

counter <- 0
for (varz in vars) {
  counter <- counter+1
  cat('Current Autotrack: ',counter,' for ',varz,'\n')
  at <- autotrack(current.track = counter, total.tracks = tracks);
  at$r1 <- at$r1-0.01
  data <- get(run)
  dat <- subset(data,ID==varz)
  
  #add pvalues
  kpRect(kp, chr = dat$chr, x0=dat$start,x1=dat$end,y0=0,y1=1,ymin=0,ymax=1,
         border=as.character(dat$DiscreteColor),col=NA,
         r0=at$r0, r1=at$r1);
  kpAddLabels(kp,cex=0.7,labels = varz,r0=at$r0, r1=at$r1,col="black",srt=0,label.margin = 0.01)
  
}

dev.off()
```

## Variant Calling

Split job by chromosome:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --cpus-per-task=10
#SBATCH --partition=biohpc_gen_production
#SBATCH --time=200:00:00
TMPDIR=$SCRATCH_DSS

genome=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Genomes/Corvus_ASM73873v5/GCF_000738735.5_ASM73873v5_genomic.fa
popfile=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Crow_Population_Resequencing_ASM73873v5/SPECIES_CROW_K28.txt
CHR=$1

mkdir raw
mkdir filtered
mkdir allsites
mkdir phased
mkdir divergence
mkdir lostruct_in

#bcftools 1.16

grep -w ${CHR} Genome.bed > intervals/${CHR}.bed

#variant calling
echo "CALLING VCF FOR ${CHR}"
bcftools mpileup -Ou -I -C 50 --threads 10 -f ${genome} -b K28.bamlist -R intervals/${CHR}.bed -a "AD,DP" | \
        bcftools call -a GQ --ploidy 2 --threads 10 -m -Ob -o raw/${CHR}.bcf
bcftools index --threads 10 raw/${CHR}.bcf

#variant filtering
echo "FILTERING AND PREPPING VCF FOR ${CHR}"
AVGDP=$(bcftools query -f '%DP\n' raw/${CHR}.bcf | datamash mean 1 | datamash round 1)
#after you run it and create QC stats, re-run this with the sample files input for the final samples you want to run
bcftools view --types snps --min-alleles 2 --max-alleles 2 --threads 10 -e "QUAL < 20 || INFO/DP>2*$AVGDP || INFO/DP < 28" -Oz -o filtered/${CHR}.IF.vcf.gz raw/${CHR}.bcf
bcftools index --threads 10 filtered/${CHR}.IF.vcf.gz

MINDP=3
#set bad genotypes to missing (require at least 2x for a genotype)
bcftools filter --threads 10 -e "FORMAT/DP < ${MINDP}" --set-GTs . filtered/${CHR}.IF.vcf.gz -Ou | bcftools view --threads 10 -U -i 'TYPE=="snp"' -Oz -o filtered/${CHR}.IF-GF.vcf.gz
bcftools view --threads 10 -i 'F_MISSING<0.1' -Oz -o filtered/${CHR}.IF-GF-MM1.vcf.gz filtered/${CHR}.IF-GF.vcf.gz
bcftools index --threads 5 filtered/${CHR}.IF-GF.vcf.gz
bcftools index --threads 5 filtered/${CHR}.IF-GF-MM1.vcf.gz

#phase vcf with beagle
java -jar ~/modules/beagle.28Jun21.220.jar gt=filtered/${CHR}.IF-GF-MM1.vcf.gz out=phased/${CHR}.IF-GF-MM1.phased nthreads=10 window=40 overlap=2 impute=true
bcftools index --threads 10 phased/${CHR}.IF-GF-MM1.phased.vcf.gz

#calculate pop gen metrics
grep -w "${CHR}" Genome.region > intervals/${CHR}.region
~/modules/Theta_D_H.Est/Theta_D_H.Est --gzvcf phased/${CHR}.IF-GF-MM1.phased.vcf.gz --region intervals/${CHR}.region --window_shift 5000@5000 --out divergence/${CHR}_popgen
~/modules/Theta_D_H.Est/Theta_D_H.Est --gzvcf phased/${CHR}.IF-GF-MM1.phased.vcf.gz --region intervals/${CHR}.region --out divergence/${CHR}_popgen-WC

#separate invariant sites
#echo "FILTERING AND MERGING INVARIANT SITES FOR ${CHR}"
bcftools view --max-ac 0 --threads 10 -Oz -o filtered/${CHR}.1N.vcf.gz raw/${CHR}.bcf
bcftools index --threads 10 filtered/${CHR}.1N.vcf.gz

#re-merge the invariant and filtered SNPs
bcftools concat --threads 10 -Ou filtered/${CHR}.1N.vcf.gz filtered/${CHR}.IF-GF.vcf.gz | bcftools sort -Oz -o allsites/${CHR}.AllSites.vcf.gz
bcftools index --threads 10 allsites/${CHR}.AllSites.vcf.gz

#create simon's divergence input file from the filtered vcf and the raw vcf
echo "CALCULATING DIVERGENCE FOR ${CHR}"
python ~/modules/genomics_general/VCF_processing/parseVCF.py --ploidyMismatchToMissing --ploidy 2 --skipIndels -i allsites/${CHR}.AllSites.vcf.gz | bgzip > allsites/${CHR}.geno.gz

#calculate DXY and FST, per window and chromosome-wide
popgenWindows.py -w 5000 -m 1000 -g allsites/${CHR}.geno.gz -o divergence/${CHR}.DXY.csv.gz -f phased --windType coordinate --ploidy 2 -T 10 -p CORNIX -p CORONE --popsFile ${popfile}

echo "${CHR} with window size ${win})"
win=$(grep -w "${CHR}" Genome.region | awk '{print $4}')
#also calculate dxy/fst on whole chromosome
popgenWindows.py -w ${win} -m 1000 -g allsites/${CHR}.geno.gz -o divergence/${CHR}.DXY-WC.csv.gz -f phased --windType coordinate --ploidy 2 -T 10 -p CORNIX -p CORONE --popsFile ${popfile}

```

Summarize SNP quality:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=36:00:00

CHR=$1

#INFO stats
echo -e "CHROM\tPOS\tQUAL\tDP\tMQBZ\tRPBZ\tSCBZ\tBQBZ" > ${CHR}.snp.qc
bcftools query -i 'TYPE="SNP"' -f '%CHROM\t%POS\t%QUAL\t%DP\t%MQBZ\t%RPBZ\t%SCBZ\t%BQBZ\n' ../${CHR}.IF.vcf.gz >> ${CHR}.snp.qc

Rscript QC.R ${CHR}

#FORMAT stats
bcftools stats -S- ../${CHR}.IF.vcf.gz | grep 'PSC' | tr ' ' '_' | awk -v c=${CHR} '{OFS="\t"}{print c, $3,$4,$5,$6,$7,$8,$10,$11,$14}' | sed '1,2d' > ${CHR}_sample.stats

```

With this R script ready to create histograms:

```R
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=36:00:00

CHR=$1

#INFO stats
echo -e "CHROM\tPOS\tQUAL\tDP\tMQBZ\tRPBZ\tSCBZ\tBQBZ" > ${CHR}.snp.qc
bcftools query -i 'TYPE="SNP"' -f '%CHROM\t%POS\t%QUAL\t%DP\t%MQBZ\t%RPBZ\t%SCBZ\t%BQBZ\n' ../${CHR}.IF.vcf.gz >> ${CHR}.snp.qc

Rscript QC.R ${CHR}

#FORMAT stats
bcftools stats -S- ../${CHR}.IF.vcf.gz | grep 'PSC' | tr ' ' '_' | awk -v c=${CHR} '{OFS="\t"}{print c, $3,$4,$5,$6,$7,$8,$10,$11,$14}' | sed '1,2d' > ${CHR}_sample.stats
di39dux@mpp3-login8:~/symlinks/atac/Crow_Population_Resequencing_ASM73873v5/bcftools/filtered/vcfqc$ cat QC.R
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
.libPaths('~/miniconda3/envs/mtdna/lib/R/library')
library(tidyverse)

qca = read.table(paste0(args[1],'.snp.qc'),header=T)
chrdat = NULL
qca[,2:8] <- sapply(qca[,2:8],as.numeric)
for (chr in unique(qca$CHROM)){
  qc = subset(qca,CHROM==chr)
  vals = NULL
  for (col in 3:ncol(qc)) {
    val = hist(qc[,col],plot=F)
    vald = data.frame(val$mids,val$counts)
    vald$variable = names(qc)[col]
    names(vald) = c('Break','Count','Variable')
    vals = rbind(vals,vald)
  }
  vals$Chr = chr
  chrdat = rbind(chrdat,vals)
}
write.table(chrdat,file=paste0(args[1],'.snpout.qc'),quote=F,sep='\t',row.names=F)

```

And plot:

```R
#### Sample stats 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Crow_Population_Resequencing_ASM73873v5/bcftools/filtered/vcfqc')
.libPaths('~/miniconda3/envs/r/lib/R/library')
library(tidyverse)
library(karyoploteR)
library(viridis)
library(LICORS)
library(ggplot2)
library(openxlsx)
library(scales)
library(gmodels)

#for chr length
genome = read.table('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Crow_Population_Resequencing_ASM73873v5/bcftools/Genome.region') %>% select(V1,V4)
names(genome) = c('chr','length')
samp = read.table('SAMPLE.info',header=F)
names(samp) = c('chr','ID','nREF','nALT','nHET','nTs','nTv','avgDP','singletons','missing')
samp = merge(samp,genome)
samp = samp %>% mutate(TsTv = nTs/nTv,
                       pMiss = missing/length,
                       hetratio = nHET/nALT,
                       hettotal = nHET/length)
#grab what you want to plot
sm = samp %>% select(chr,ID,nALT,nHET,avgDP,singletons,TsTv,pMiss,hetratio,hettotal) %>% pivot_longer(!c(chr,ID))
sm$value = as.numeric(sm$value)
smp = sm %>% separate(ID,into=c('Taxon','Locality','IDnum'),remove=F)
smp = smp %>% mutate(Taxon = gsub('^D','C.c.corone',Taxon),
                     Taxon = gsub('^S','C.c.cornix',Taxon))
namesold = unique(smp$name)
namesnew = c('Num_Alt_Alleles','Num_Het_Sites','Average_DP','Num_Singletons','TsTv','Percent_Missing','Het_Ratio','Het_Total')
renames = data.frame(name = namesold,vari = namesnew)
smp = left_join(smp,renames) %>% select(-name)
smps = smp %>% group_by(ID,Taxon,vari) %>%
  summarize(m = ci(value)[1],
            l = ci(value)[2],
            u = ci(value)[3])
spp = smps %>% 
  #filter(name == 'nHET' | name == 'singletons' | name =='pMiss' | name =='hetratio') %>% 
  ggplot(aes(x=ID,y=m,ymin=l,ymax=u,col=Taxon))+
  geom_point()+
  geom_errorbar()+
  facet_wrap(vari~.,scales='free',nrow = 4)+
  scale_color_manual(values=c('grey20','grey60'))+
  theme_classic()+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
pdf('Site_Summaries.pdf',height=8,width=9)
spp
dev.off()

#### inspect plots, find individuals that grossly violate distributions

#site info
sf = read.table('SITE.info',header=T)
sfs = sf %>% group_by(Break,Variable) %>%
  summarize(m = ci(Count)[1],
            l = ci(Count)[2],
            u = ci(Count)[3])
sfp = sf %>% subset(Chr=='chr1' | Chr =='chr2' | Chr =='chr3' | Chr == 'chr18') %>%
  ggplot(aes(x=Break,y=Count,col=Variable))+
  geom_line()+
  xlab('Values')+
  scale_color_manual(values=viridis(6))+
  facet_wrap(Chr~Variable,scales='free',ncol=6)+
  theme_classic()
pdf('Site_info.pdf',height=6,width=15)
sfp
dev.off()

```



## Identify Transitions

We will subset our variant-sites for C-T and A-G SNPs, that we will remove from our bisulfite-converted datasets, since they will be incorrectly inferred as methylated sites. We convert the vcf to zero-based bed format by subtracting '1' from the vcf coordinate: 

```bash
zcat *V4.vcf.gz | awk 'BEGIN{OFS="\t"}{print $1, $2-1, $2, $4, $5}' | grep -v '#' > GCF_000738735.5_ASM73873v5_genomic.AllSNPs.bed

cat GCF_000738735.5_ASM73873v5_genomic.AllSNPs.bed | wc -l
1470103

head GCF_000738735.5_ASM73873v5_genomic.AllSNPs.bed
chr6    3621    3622    T       C
chr6    6774    6775    C       T
chr6    6780    6781    C       T
chr6    9278    9279    C       T
chr6    9384    9385    T       A
chr6    9409    9410    G       A
```

Now we subset CT and AG SNPs, which could be incorrectly identified as 5mC loci in BS experiments. 

```bash
awk '(($4 == "C" && $5 == "T") || ($4 == "T" && $5 == "C") || ($4 == "G" && $5 == "A") || ($4 == "A" && $5 == "G") || $0 ~ /^#/ )' GCF_000738735.5_ASM73873v5_genomic.AllSNPs.bed | awk '{OFS="\t"}{print $1, $2, $3}' > GCF_000738735.5_ASM73873v5_genomic.TS_SNPs.bed

cat GCF_000738735.5_ASM73873v5_genomic.TS_SNPs.bed | wc -l
995822
```

Which leaves us with 995,822 transitions and 474,281 sites for masking, giving us a Ts/Tv of 2.10).

Simon Martin's estimation creates some negative FST measurements -- we will set these simply to zero. Also, any NaN FST values, which are present in Dxy, are simply zero (no diveristy, so FST errors). I set these specific cases to FST=0.

```r
library(LICORS)
library(dplyr)

#import data
vcf <- read.table("Master.gen",header=FALSE)
names(vcf) <- c('CHROM','START','END','Dxy','FST','PiAll')
vcf$FST <- threshold(vcf$FST,min=0,max=1)
vcf$Dxy <- threshold(vcf$Dxy,min=0,max=1)

vcf[is.na(vcf$Dxy),]
vcfd <- vcf[!is.na(vcf$Dxy),]
vcfd$FST[is.na(vcfd$FST)] <- 0;
options(scipen=999)
write.table(vcfd,file='Corvus_K28_Gen.bed',row.names=F,quote=F,sep='\t')
```

Leaving a total of 2,023,280 windows. 

## Population Genetics

The above variant calling script should have output many popgen and dxy files. Combine them and analyze in R below. Statistics were called either with the Shuhua group's github `https://github.com/Shuhua-Group/Theta_D_H.Est` or Simon Martin's (for FST / DXY).

Citation for hap diveristy can be found within Nei, or [here](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0251878#pone.0251878.ref009) 

Assess correlations between statistics, identify uncorrelated ones to retain (note that the tidymodels package will drop automatically anyways):

```R
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Crow_Population_Resequencing_ASM73873v5/bcftools')
.libPaths('~/miniconda3/envs/r/lib/R/library')
library(tidyverse)
.libPaths('~/miniconda3/envs/mtdna/lib/R/library')
library(LICORS)
library(GGally)
library(corrr)
library(gt)

as = read.table('divergence/Popgen.txt',header=T)
dx = read.table('divergence/DXY.txt',header=T)
d1 = dx %>% select(!c(mid))
names(d1) = c('chr','start','end','sites','pi_CORNIX','pi_CORONE','Dxy','FST')
#drop any with no pi data
d2 = d1 %>% drop_na(Dxy)
#now, FST is NA only when dxy/pi is 0 (no variance), set these to fst 0 because they are completely undiverged
#first double check
d2[is.na(d2$FST),] %>% summarize(max=max(Dxy))
d2$FST[is.na(d2$FST)] <- 0;
options(scipen=999)
#set threshold to fst between 0 and 1
d2 = d2 %>% mutate(FST = threshold(FST,min=0,max=1))

#merge
as = as %>% select(-c(regionID,Hfaywu,norm_Hfaywu))
fu = full_join(as,d2)

#quick cor
fucorL = fu %>% 
  select(singleton,ThetaPI,ThetaK,segregating,Hap_diversity,Ffuli,Dfuli,Dtajima,Dxy,FST) %>% 
  correlate(method = 'spearman')
f1 = fucorL %>% select(-term)
new <- matrix(NA, nrow = 10, ncol = 10)
new[upper.tri(new)] <- f1[upper.tri(f1)]
cordf = new %>% as.data.frame
names(cordf) = names(f1)
cordf$term = names(cordf)
cordf = cordf %>% select(term,everything())
cordf
cordf %>%
  gt() %>%
  fmt_number(columns = !(term),
             rows = everything(),
             n_sigfig = 3)

write.table(cordf,file = 'Rho_PopGen.txt',quote=F,sep='\t',row.names=F)
rplot(cordf)

#on final subset
fucors = fu %>% 
  select(Hap_diversity,Dfuli,Dtajima,Dxy,FST) %>% 
  correlate(method = 'spearman')
fucors[upper.tri(fucors )] = 0
fucors  %>%
  gt() %>% 
  fmt_number(columns = !(term),
             rows = everything(),
             n_sigfig = 3)
rplot(fucors)
write.table(fucors,file = 'Rho_PopGenUsed-LD.txt',quote=F,sep='\t',row.names=F)

#correlations 
png('Correlations-Reseq.png',units='in',res=300,bg='transparent',height=10,width=12)
g <- fu %>% 
  select(ThetaK,Dfuli,Dtajima,Dxy,FST) %>% 
  ggpairs(upper = list(continuous = "cor", corMethod = "spearman") + theme_minimal())
print(g)
dev.off()

#save final output
fs = fu %>% select(c(chr,start,end,Hap_diversity,Dfuli,Dtajima,Dxy,FST))
write.table(fs,file='Corvus28_PopGen-2023JAN04.txt',quote=F,sep='\t',row.names=F)

#calculate watterson's theta manually to ensure it's consistent (num segreg sites / harmonic mean of chrs) 
harmonicNumber = 0
numChromosomes = 56
for (i in 1:(numChromosomes - 1)) {
  harmonicNumber = harmonicNumber + 1.0/i
}
print(harmonicNumber)
```

Calculate thetaK / Da: 

```R
library(tidyverse)
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Crow_Population_Resequencing_ASM73873v5/bcftools')
d = read.table('Corvus28_AllPopGen-2022NOV18_Unlinked.txt',header=T)
names(d)
#da = dxy - (intrapop variation) / 2
d = d %>% mutate(Da = Dxy - (pi_CORNIX + pi_CORONE) / 2,
                 DaTheta = Da/ThetaK)
range(d$DaTheta,na.rm=T)
d = d %>% mutate(AvZvMT = ifelse(chr == 'chrMT','MT',
                                 ifelse(chr == 'chrZ','Z',
                                        'Autosome')))
d %>% ggplot(aes(x=DaTheta,fill=AvZvMT))+
  geom_histogram(alpha=0.5)+
  facet_grid(AvZvMT~.,scales='free')
  theme_classic()
```

### Focal Region Divergence

```R
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/crow_scratch/2021-07_HybridZone/Nov/Classification_Sept')
.libPaths('~/miniconda3/envs/r/lib/R/library')
library(tidyverse)
library(karyoploteR)
library(viridis)
library(LICORS)
library(ggplot2)
library(scales)
library(rsample)
library(rstatix)
library(ggpubr)
library(stringr)
library(gmodels)
library(ggpol)

set.seed(123)
g = read.table('GCF_000738735.5_ASM73873v5_genomic.K28Corvus-FullPopGen.5000bp.bed',header=TRUE)
g %>% mutate(thresh = 0.3) %>% filter(FST > thresh & chr == 'chr18') %>% 
  summarize(min = min(Genstart),
            max = max(Genend))  #this identifies the focal region, the min/max of chr18 which contains FST peaks > 0.3
g %>% mutate(thresh = 0.3) %>% filter(FST > thresh) %>% count(chr)
g = g %>% mutate(Focal_Region = ifelse(chr == 'chr18' & Genstart > 8.07e6 & Genstart < 10.07e6,'Focal','Undiverged'))
gm = g %>% select(-c(chr,Genstart,Genend,site,start,Region,GC,Site)) %>% unique %>% pivot_longer(!(Focal_Region))
gm %>% group_by(name) %>% 
  summarize(mean = mean(value,na.rm=TRUE),
            min = min(value,na.rm=TRUE),
            max = max(value,na.rm=TRUE))
gdat = NULL
options(dplyr.summarise.inform = FALSE)
for (i in seq(1,1000)){
  cat('Iteration: ',i,'\n')
  nsamp = gm %>% group_by(Focal_Region) %>% 
    filter(name == 'HapD' & Focal_Region == 'Focal') %>% count(Focal_Region)
  b1 = gm %>%
    group_by(Focal_Region,name) %>% 
    slice_sample(n=round(nsamp$n/2,0),replace=TRUE) %>% 
    summarize(mean = mean(value,na.rm=TRUE))
  gdat = rbind(gdat,b1)
  
}

write.table(gdat,file='BootstrappedGenetic.boots',quote=F,sep='\t',row.names=F)
gdat = read.table('BootstrappedGenetic.boots',header=TRUE)
gdats = gdat %>% ggplot(aes(y=mean,x=name,fill=Focal_Region))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA, 
                 jitter.height = 0.05, jitter.width = 0.075, errorbar.draw = TRUE) +
  scale_fill_manual(values=c('darkviolet','grey40'))+
  ylab('Bootstrapped Mean')+xlab('')+
  facet_wrap(.~name,scales='free')+
  scale_color_manual(values=c('darkviolet','grey40'))+
  theme_bw()

pdf('GeneticVariation_FocalBG.pdf',height=4,width=7)
gdats
dev.off()

#points
gs = gdat %>% 
  group_by(Focal_Region,name) %>% 
  summarise(m = ci(mean)[1], 
            lo = ci(mean)[2],
            hi = ci(mean)[3],
            sd = ci(mean)[4])

write.table(gs,file=paste0('Genetic_Bootstraped_ResultsFocal-2023JAN27.txt'),quote=F,sep='\t',row.names=F)

```



# Spatial Distribution

Plot the distributions of our crows. 

```R
# Local on lenovo
install.packages(c('cowplot', 'googleway', 'ggplot2', 'ggrepel', 'ggspatial', 'libwgeom', 'sf', 'rnaturalearth', 'rnaturalearthdata'))
setwd('G:/My Drive/Research/Crow_Hybrid_Epigenetics/spatial/2025_01')
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(viridis)
library(ggspatial)
library(tidyverse)

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
ggplot(data = world) +
  geom_sf(color = "black", fill = "palegreen")

cd <- read_tsv('../../All-Metadata.txt')
hz = cd %>% filter(Experiment == 'HZ')
sites <- st_as_sf(hz, coords = c("Longitude", "Latitude"), 
                  crs = 4326, agr = "constant")
plum =
  ggplot(data = world) +
  geom_rect(aes(xmin = min(hz$Longitude)-5, xmax = max(hz$Longitude)+5, 
                ymin = min(hz$Latitude)-1, ymax=max(hz$Latitude)+1),fill='lightblue1')+
  geom_sf(fill = "white") +
  geom_sf(data = sites, size = 6, aes(fill = Chr18_Hybrid_Index, shape=as.factor(Capture_Year)),col='black') +
  scale_fill_gradient(low='grey30',high='grey90')+
  scale_shape_manual(values=c(24,21,23))+
  annotation_scale(location = "br", width_hint = 0.3,pad_x = unit(0.3, "in")) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.3, "in"), pad_y = unit(0.25, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(min(cd$Longitude)-5, max(cd$Longitude)+5), 
           ylim = c(min(cd$Latitude)-1, max(cd$Latitude)+1), expand = FALSE)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
plum

ggsave('20250107_Crow_Distribution.pdf',plum,dpi=300,height=5,width=8)

## Focus only on hybrids
hybs <- cd %>% filter(Experiment == 'HZ' & grepl('H_',Short_ID))
hybs <- hybs %>% mutate(Latitude = jitter(Latitude,amount=0.1),
                        Longitude = jitter(Longitude,amount=0.1))
sites_hybs <- st_as_sf(hybs, coords = c("Longitude", "Latitude"), 
                  crs = 4326, agr = "constant")
plum_hyb =
  ggplot(data = world) +
  geom_sf(fill = "white") +
  geom_sf(data = sites_hybs, size = 6, aes(fill = Chr18_Hybrid_Index, shape=as.factor(Capture_Year)),col='black') +
  scale_fill_gradient(low='grey30',high='grey90')+
  scale_shape_manual(values=c(24,21,23))+
  annotation_scale(location = "br", width_hint = 0.3,pad_x = unit(0.3, "in")) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.3, "in"), pad_y = unit(0.25, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(min(hybs$Longitude)-2, max(hybs$Longitude)+2), 
           ylim = c(min(hybs$Latitude)-1, max(hybs$Latitude)+1), expand = FALSE)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
plum_hyb

ggsave('20250107_Crow_Distribution-Hybrids.pdf',plum,dpi=300,height=4,width=6)


```

Hybrid zoom:

```R
#install.packages(c('cowplot', 'googleway', 'ggplot2', 'ggrepel', 'ggspatial', 'libwgeom', 'sf', 'rnaturalearth', 'rnaturalearthdata'))
setwd('F:/Research/scratch/crow_hybrid_paper')
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(sf)
library(ggplot2)
library(viridis)
library(ggspatial)

set.seed(123)
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
ggplot(data = world) +
  geom_sf(color = "black", fill = "palegreen")

cd <- read.table('All-Metadata.txt',header=TRUE)
hz = subset(cd,Retained==1 & Experiment == 'HZ' & Species == 'Hybrid')
hz = hz %>% mutate(Longitude= jitter(Longitude,25),
                   Latitude= jitter(Latitude,25))
sites <- st_as_sf(hz, coords = c("Longitude", "Latitude"), 
                  crs = 4326, agr = "constant")
plum =
  ggplot(data = world) +
  geom_rect(aes(xmin = min(hz$Longitude)-1, xmax = max(hz$Longitude)+1, 
                ymin = min(hz$Latitude)-1, ymax=max(hz$Latitude)+1),fill='lightblue1')+
  geom_sf(fill = "white") +
  geom_sf(data = sites, size = 8, aes(fill = Temperature, shape=as.factor(Species)),col='black') +
  scale_fill_gradient(low='blue',high='yellow')+
  #scale_fill_gradient(low='grey10',high='grey75')+
  scale_shape_manual(values=c(23))+
  annotation_scale(location = "br", width_hint = 0.3,pad_x = unit(0.3, "in")) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.3, "in"), pad_y = unit(0.25, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(min(hz$Longitude)-1, max(hz$Longitude)+1), 
           ylim = c(min(hz$Latitude)-0.5, max(hz$Latitude)+0.5), expand = FALSE)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_x_continuous(breaks = seq(from = 11.5, to = 15.5, by = 2))+
  scale_y_continuous(breaks = seq(from = 50.5, to = 53, by = 2))

plum

png('Crow_Distribution-2023JAN31_HYBRIDS.png',units='in',res=600,height=3,width=4,bg='transparent')
plum
dev.off()

```



## Resequencing Individuals

For the resequencing data, we can plot the coordinates, and also calculate an IBD matrix. 

```r
#install.packages(c('cowplot', 'googleway', 'ggplot2', 'ggrepel', 'ggspatial', 'libwgeom', 'sf', 'rnaturalearth', 'rnaturalearthdata'))
setwd('F:/Research/scratch/crow_hybrid_paper/2022_02/Spatial')
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggsci)
library(viridis)
library(ggplot2)
library(reshape2)
library(rgeos)
library(forcats)
library(ggmap)
library(tidyverse)
library(lwgeom)
library(cluster)
library(factoextra)
source('F:/Research/scratch/R_Functions.R')

#load in cuckoo data
all <- read.table('Reseq.txt',header=T)
plm <- all

#convert
sites <- st_as_sf(plm, coords = c("Longitude", "Latitude"), 
                  crs = 4326, agr = "constant")
world <- ne_countries(scale = "medium", returnclass = "sf")

#show full range, colored by host and species PRE-CLUSTERING BY DISTANCE 
plum <- ggplot(data = world) +
  geom_sf() +
  geom_sf(data = sites, size = 6, aes(col = Taxon, shape=Locality)) +
  scale_color_manual(values=c('grey10','grey60'))+
  coord_sf(xlim = c(min(plm$Longitude)-5, max(plm$Longitude)+5), 
           ylim = c(min(plm$Latitude)-5, max(plm$Latitude)+5), expand = FALSE)
plum


#try with topographic
bbox <- c(left = min(plm$Longitude)-5, bottom = min(plm$Latitude)-5, right = max(plm$Longitude)+5, top = max(plm$Latitude)+5)
shapes <- c(15,16,8,17)
basemap <- ggmap(get_stamenmap(bbox,maptype = "terrain-background", zoom = 5))+
  theme_bw() + 
  labs(x = "Longitude", y = "Latitude")

#add points ( filled in ) 
shapes <- c(21,22,24,25)
fm <- basemap + 
  geom_point(data=plm,aes(x=Longitude,y=Latitude,fill=Taxon,col='black',shape=Locality),size=6,
             position=position_jitter(width=0,height=0),show.legend = T)+
  scale_shape_manual(values=shapes)+
  scale_color_manual(values='black')+
  scale_fill_manual(values=c('grey10','grey60'))+
  guides(fill=guide_legend(override.aes=list(shape=21)))
fm

png('Reseq_Spatial.png',units='in',height=6,width=10,res=600,bg='white')
fm
dev.off()


dist = st_distance(sites$geometry,  by_element = F)
d1 <- as.data.frame(dist)
d2 <- as.data.frame(lapply(d1, function(y) gsub(" [m]", "", y)))
names(d2) <- sites$ID
row.names(d2) <- sites$ID
d3 <- as.dist(as.matrix(d2))

#perform nmds
set.seed(123)
gmt <- as.data.frame(t(d3))
#read ids
ids <- read.table('K20.id')
names(gmt) <- ids$V1
gmut <- t(gmt)

nmds = metaMDS(d3, distance = "bray", trymax = 100)
nmds
data.scores <- as.data.frame(scores(nmds))
data.scores$ID <- rownames(data.scores)
data.scores <- merge(data.scores,variables,by='ID')
data.scores <- data.scores %>% mutate(Taxon = ifelse(grepl('^D_',ID)==TRUE,'C.C.','H.C.'))
data.scores <- data.scores %>% separate(ID,into=c('Species','Locality','Number'))

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

png('Reseq_Distance_NMDS.png',units='in',height=6,width=10,res=600,bg='white')
np
dev.off()


#traditional PCA 
#perform pca..
pca <- prcomp(d3,scale=TRUE)
var_explained <- pca$sdev^2/sum(pca$sdev^2)
vec <- data.frame(pca$x)
vec$ID <- row.names(vec)
vec <- vec %>% mutate(Taxon = ifelse(grepl('^D_',ID)==TRUE,'C.C.','H.C.'))
vec <- vec %>% separate(ID,into=c('Species','Locality','Number'))

np <- 
  ggplot() +
  geom_point(data = vec, show.legend = T,aes(x=PC1,
                                                     y=PC2, 
                                                     shape = Locality,
                                                     colour = Taxon,
                                                     fill = Taxon,
                                                     size=15)) +
  scale_color_manual(values=c('grey10','grey60'))+
  scale_fill_manual(values=c('grey10','grey60'))+
  scale_shape_manual(values=c(16,15,18,17))+
  xlab(paste0('PC1 : ',round(var_explained[1],3)*100,'% Explained'))+
  ylab(paste0('PC2 : ',round(var_explained[2],3)*100,'% Explained'))+
  theme_classic(base_size=14)
np

png('Reseq_Distance_PCA.png',units='in',height=6,width=10,res=600,bg='white')
np
dev.off()
```

# Relatedness

I will call SNPs with BS-SNPer to see relatedness across all our individuals:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=24:00:00

#modify the v5.7 genome (ASM) to maximum 80 characters per line for BSSNPer
genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/crow_scratch/2021-07_HybridZone/Nov/BSSNPs
genome=${genomedir}/genome.80.fa

RUN=$1

#call BS snps
perl ~/modules/BS-Snper-master/BS-Snper.pl --fa ${genome} --input ${RUN}.sorted.bam --output ${RUN} \
        --methcg ${RUN}.CG --methchg ${RUN}.CHG --methchh ${RUN}.CHH \
        --mincover 10 --maxcover 500 --mapvalue 20 --minread2 2 --errorate 0.02 --minquali 15 --minhetfreq 0.1 --minhomfreq 0.85 \
        > ${RUN}.VCF
```

```bash
xp=$1

for i in $(ls *VCF); do bgzip -c ${i} > ${i}.gz; tabix ${i}.gz;  done
bcftools merge -Ou -o - *.VCF.gz | bcftools view --min-alleles 2 --max-alleles 2 -V indels -Oz -o $xp.vcf.gz
bcftools query -l $xp.vcf.gz | sed 's/\..*//g' > Samps.list
bcftools reheader --samples Samps.list -o ${xp}.2N.vcf.gz ${xp}.vcf.gz
bcftools view -Ou -o - --min-ac 2 ${xp}.2N.vcf.gz | \
    bcftools view --exclude "F_MISSING>0.2" -Ov -o - - | \
    sed 's/VCFv4.3/VCFv4.2/g' | \
    bgzip -c > ${xp}.2N_MAC2-MM2.vcf.gz

#relatedness using KING's statistic    
vcftools --gzvcf ${xp}.2N_MAC2-MM2.vcf.gz --relatedness2 --out ${xp}.2N_MAC2-MM2

#and also HC vs CC fst
bcftools query -l ${xp}.2N_MAC2-MM2.vcf.gz > Samps.list
grep 'D_Ko' Samps.list > CC.list
grep 'S_Up' Samps.list > HC.list
vcftools --gzvcf ${xp}.2N_MAC2-MM2.vcf.gz --weir-fst-pop CC.list --weir-fst-pop HC.list --out ${xp}.2N_MAC2-MM2

#format experiment
awk -v x=${xp} '{OFS="\t"}{print $0,x}' ${xp}.2N_MAC2-MM2.relatedness2 > ${xp}_RELATE.txt
awk -v x=${xp} '{OFS="\t"}{print $0,x}' ${xp}.2N_MAC2-MM2.weir.fst > ${xp}_FST.txt
```

From the BS-called SNPs, plot relatedness heatmaps and a manhattan plot for FST.

```r
#### Plot Relatedness & FST
setwd('~/files/class')
source('~/R_Functions.R')
library(dplyr)
library(karyoploteR)
library(tidyverse)
library(viridis)
library(LICORS)
library(ggplot2)
library(scales)
library(gmodels)
library(reshape2)
library(ggpubr)

#relate first
rel <- read.table('Master.rel',header=T)
rel <- rel %>% select(INDV1,INDV2,RELATEDNESS_PHI,Experiment)
names(rel) <- c('ID1','ID2','PHI','XP')
rel <- rel %>% mutate(XP=gsub('CG','Co.Gar.RRBS',XP),
                        XP=gsub('HZ','Hyb.Zo.RRBS',XP),
                        XP=gsub('WGBS','Co.Gar.WGBS',XP))

for (xp in unique(rel$XP)) {
  cat('Generating plot for: ',xp,'\n')
  d <- subset(rel,XP==xp)
  dp <- ggplot(d,aes(x=ID1,y=ID2,fill=PHI))+
    geom_tile()+
    scale_fill_gradient(low='yellow',high='red')+
    theme_minimal()+xlab('')+ylab('')+
    theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1))+
    ggtitle(xp)
  assign(xp,dp)
}
pdf('Relatedness-BSSNP.pdf',height=9,width=12)
ggarrange(Co.Gar.WGBS,Co.Gar.RRBS,Hyb.Zo.RRBS)
dev.off()

#plot FST
#set up genome
genome <- read.table('genome5.7.bed',header=FALSE) %>% filter(str_detect(V1,'scaff|MT',negate=T)) %>% arrange(desc(V2))
names(genome) <- c('chr','end'); genome$start <- 1
genome$chr <- gsub('chr','',genome$chr)
crowG <- makeGRangesFromDataFrame(genome,seqnames.field = 'chr',start.field = 'start',end.field = 'end')
#for chr 18
crow18 <- makeGRangesFromDataFrame(subset(genome,chr=='18'),seqnames.field = 'chr',start.field = 'start',end.field = 'end')

#import fst 
divs <- read.table('Master.fst',header=T)
names(divs) <- c('chr','start','FST','experiment')
divs$FST <- threshold(divs$FST,min=0,max=1)
divs$chr <- gsub('chr','',divs$chr)
divs <- divs %>% mutate(experiment=gsub('CG','Co.Gar.RRBS',experiment),
                        experiment=gsub('HZ','Hyb.Zo.RRBS',experiment),
                        experiment=gsub('WGBS','Co.Gar.WGBS',experiment))

gb1 = genome[seq(1, nrow(genome), 2), ][[1]]
gb2 = genome[seq(2, nrow(genome), 2), ][[1]]

#assign colors based on chromosome
div_gb1 <- divs[grepl(paste(gb1, collapse="$|"), divs$chr), ]
div_gb1$Color <- "grey60"
div_gb2 <- divs[grepl(paste(gb2, collapse="$|"), divs$chr), ]
div_gb2$Color <- "grey10"
divs_use <- rbind(div_gb1,div_gb2)
divs_use$experiment <- factor(divs_use$experiment,levels=c('Co.Gar.WGBS','Co.Gar.RRBS','Hyb.Zo.RRBS'))

##### Plot chr18 #####
png('Chr18_BS-SNPer-Manhattan.png',units='in',res=600,height=6,width=6,bg='transparent')

#Plot Base, add labels.plotter=NULL for no chr names
pp <- getDefaultPlotParams(plot.type=4)
pp$leftmargin <- 0.08
kp <- plotKaryotype(plot.type=4, genome = crow18,plot.params = pp,labels.plotter=NULL)

kpAddBaseNumbers(kp,tick.dist = 5000000,minor.tick.dist = 1000000,cex=0.65)

run <- 'divs_use'
vars <- unique(divs_use$experiment)
tracks <- length(vars)

counter <- 0
for (varz in vars) {
  counter <- counter+1
  cat('Current Autotrack: ',counter,' for ',varz,'\n')
  at <- autotrack(current.track = counter, total.tracks = tracks);
  at$r1 <- at$r1-0.05
  data <- get(run)
  #if chr 18 unhash
  data <- subset(data,chr=='18')
  dat <- data[grepl(varz,data$experiment),];

  #plot fst
  kpPoints(kp, chr = dat$chr, x=dat$start, y=dat$FST,
           col="grey20",
           pch=16,
           ymin=0,ymax=1,r0=at$r0, r1=at$r1,cex=0.75);
  kpAxis(kp, ymin =0, ymax=1,cex=0.8,r0=at$r0, r1=at$r1,col="black");
  kpAddLabels(kp,cex=0.8,labels = paste0(varz),r0=at$r0+0.2, r1=at$r1,col="black",srt=90,label.margin = 0.06)
  kpAbline(kp, h=0, r0=at$r0, r1=at$r1, col="black", lwd=1,lty=1)


}

dev.off()

##### Plot GW ####
png('WG_BS-SNPer-Manhattan.png',units='in',res=600,height=7,width=10,bg='transparent')

#Plot Base, add labels.plotter=NULL for no chr names
pp <- getDefaultPlotParams(plot.type=4)
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

run <- 'divs_use'
vars <- unique(divs_use$experiment)
tracks <- length(vars)

counter <- 0
for (varz in vars) {
  counter <- counter+1
  cat('Current Autotrack: ',counter,' for ',varz,'\n')
  at <- autotrack(current.track = counter, total.tracks = tracks);
  at$r1 <- at$r1-0.05
  data <- get(run)
  dat <- data[grepl(varz,data$experiment),];

  #plot fst
  kpPoints(kp, chr = dat$chr, x=dat$start, y=dat$FST,
           col=dat$Color,
           pch=16,
           ymin=0,ymax=1,r0=at$r0, r1=at$r1,cex=0.75);
  kpAxis(kp, ymin =0, ymax=1,cex=0.8,r0=at$r0, r1=at$r1,col="black");
  kpAddLabels(kp,cex=0.8,labels = paste0(varz),r0=at$r0+0.2, r1=at$r1,col="black",srt=90,label.margin = 0.06)
  kpAbline(kp, h=0, r0=at$r0, r1=at$r1, col="black", lwd=1,lty=1)


}

dev.off()
```

# 5mC Processing

## ComGar

Script that will take raw data and extract it, NO deduplication because it's RRBS data. 

Prepare genome for mapping:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=24:00:00

SCRATCH=/tmp/$SLURM_JOB_ID

genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Genomes/Corvus_ASM73873v5
genome=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.fa
cgi=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.CGI.bed
repeats=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.REPEATS.bed
cutsites=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.CUT-SITES.bed
ctsnp=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.TS_SNPs.bed
chrsbed=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.CHRS.bed

#working directories
rawdataWGBS=/dss/dsslegfs01/pr53da/pr53da-dss-0018/rawdata/Bisulfite-Seq/WGBS
rawdataRRBS=/dss/dsslegfs01/pr53da/pr53da-dss-0018/rawdata/Bisulfite-Seq/RRBS
basedir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/2021-07_HybridZone
outdir=${basedir}/5mC_raw

bismark_genome_preparation ${genomedir}
```

And the full pipeline:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=24:00:00

SCRATCH=/tmp/$SLURM_JOB_ID

genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Genomes/Corvus_ASM73873v5
genome=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.fa
cgi=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.CGI.bed
repeats=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.REPEATS.bed
cutsites=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.CUT-SITES.bed
ctsnp=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.TS_SNPs.bed
chrsbed=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.CHRS.bed

#working directories
rawdataWGBS=/dss/dsslegfs01/pr53da/pr53da-dss-0018/rawdata/Bisulfite-Seq/WGBS
rawdataRRBS=/dss/dsslegfs01/pr53da/pr53da-dss-0018/rawdata/Bisulfite-Seq/RRBS
basedir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/2021-07_HybridZone
outdir=${basedir}/5mC_raw/cg
workdir=${basedir}/scratch/cg

#positional arguments
RUN=$1

mkdir ${workdir}
mkdir ${outdir}

cat ${rawdataRRBS}/${RUN}*__SE__*fastq.gz > $SCRATCH/${RUN}.SE.fastq.gz

cd ${workdir}

#Trim adapters
trim_galore --fastqc -j 10 --rrbs --quality 20 --output_dir ${workdir} ${SCRATCH}/${RUN}.SE.fastq.gz

#Map reads
bismark --parallel 6 --output_dir ${workdir} --genome ${genomedir} ${workdir}/${RUN}.SE_trimmed.fq.gz 

#Overlap bam with bed
bedtools intersect -a ${workdir}/${RUN}.SE_trimmed_bismark_bt2.bam -b ${cutsites} > ${workdir}/${RUN}.cut-sites.bam

#Extract methylation
bismark_methylation_extractor --parallel 7 --ignore 2 --ignore_3prime 1 --gzip --bedgraph --buffer_size 20G --cytosine_report --genome_folder ${genomedir} --output ${outdir} ${workdir}/${RUN}.cut-sites.bam
```

And create a bash script for each sample:

```
for sample in $(cat CG-Samples.list); do sbatch -J ${sample} RRBS_Pipeline_CG.sh ${sample}; done 
```

Now, subset the sites based on 

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=2:00:00

SCRATCH=/tmp/$SLURM_JOB_ID

genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Genomes/Corvus_ASM73873v5
genome=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.fa
cgi=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.CGI.bed
repeats=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.REPEATS.bed
cutsites=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.CUT-SITES.bed
ctsnp=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.TS_SNPs.bed
chrsbed=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.CHRS.bed

#working directories
basedir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/crow_scratch/2021-07_HybridZone
outdir=${basedir}/Sept/cg
workdir=${basedir}/5mC_raw/cg_sub

#positional arguments
RUN=$1

#Remove any positions overlapping a transition
zcat $workdir/${RUN}.cut-sites.bismark.cov.gz | bedtools subtract -A -a - -b ${ctsnp} > ${workdir}/${RUN}.filt_cutsite_ctsnp.tmp

#Only keep sites on the major chromosomes
bedtools intersect -wb -a ${chrsbed} -b ${workdir}/${RUN}.filt_cutsite_ctsnp.tmp | awk '{OFS="\t"}{print $4, $5, $6, $7, $8, $9}' | sort | uniq | bgzip -c > ${outdir}/${RUN}.CpG_5mC.cov.gz
```

And count sites:

```bash
mkdir counts
echo -e "Sample\tFilter\tCount" > counts/header.txt

for i in $(cat Samples.list); do 

#make a file with the first 2 columns
echo -e "${i}\n${i}\n${i}\n${i}" > counts/${i}.sample
echo -e "Cutsite\nTsSNP\nChr\nPostCoverage" > counts/${i}.field


#now for the counts
zcat ${i}.cut-sites.bismark.cov.gz | wc -l > counts/${i}.counts
cat ${i}.filt_cutsite_ctsnp.tmp | wc -l >> counts/${i}.counts
zcat ${i}.CpG_5mC.cov.gz | wc -l >> counts/${i}.counts
zcat ${i}.CpG_5mC.cov.gz | awk '{print $5+$6}' | awk '$1 > 5' | wc -l >> counts/${i}.counts

paste counts/${i}.sample counts/${i}.field counts/${i}.counts > counts/${i}.count

done
cat counts/header.txt counts/*.count > 5mC_Counts.txt
rm -rf counts
```

## HybZon

Our PE data from the Hybrid Zone will be treated largely similar, except it is paired end, so will require trimming from both sides.

And full pipeline, we have multiple libraries for each sample, so we will merge the fastqs beforehand: 

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=24:00:00

SCRATCH=/tmp/$SLURM_JOB_ID

genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Genomes/Corvus_ASM73873v5
genome=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.fa
cgi=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.CGI.bed
repeats=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.REPEATS.bed
cutsites=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.CUT-SITES.bed
ctsnp=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.TS_SNPs.bed
chrsbed=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.CHRS.bed

#working directories
rawdataWGBS=/dss/dsslegfs01/pr53da/pr53da-dss-0018/rawdata/Bisulfite-Seq/WGBS
rawdataRRBS=/dss/dsslegfs01/pr53da/pr53da-dss-0018/rawdata/Bisulfite-Seq/RRBS
basedir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/crow_scratch/2021-07_HybridZone
outdir=${basedir}/5mC_raw/hz_sub
workdir=${basedir}/scratch/hz

#positional arguments
RUN=$1

mkdir ${workdir}
mkdir ${outdir}

cat ${rawdataRRBS}/${RUN}*__R1__*fastq.gz > $SCRATCH/${RUN}.R1.fastq.gz
cat ${rawdataRRBS}/${RUN}*__R2__*fastq.gz > $SCRATCH/${RUN}.R2.fastq.gz

cd ${SCRATCH}

#Trim adapters
trim_galore --fastqc -j 10 --rrbs --paired --quality 20 --output_dir ${workdir} ${SCRATCH}/${RUN}.R1.fastq.gz ${SCRATCH}/${RUN}.R2.fastq.gz

#Map reads
bismark --parallel 6 --output_dir ${workdir} --genome ${genomedir} -1 ${workdir}/${RUN}.R1_val_1.fq.gz -2 ${workdir}/${RUN}.R2_val_2.fq.gz

#Overlap bam with bed, since these are paired-end files, we need to use a different tool than the CG 
bedtools pairtobed -abam ${workdir}/${RUN}.R1_val_1_bismark_bt2_pe.bam -b ${cutsites}  > ${workdir}/${RUN}.cut-sites.bam

#Extract methylation
bismark_methylation_extractor --parallel 7 --ignore 2 --ignore_3prime 1 --ignore_r2 2 --ignore_3prime_r2 1 --gzip --bedgraph --buffer_size 20G --cytosine_report --genome_folder ${genomedir} --output ${outdir} ${workdir}/${RUN}.cut-sites.bam

#Remove any positions overlapping a transition
zcat $outdir/${RUN}.cut-sites.bismark.cov.gz | bedtools subtract -A -a - -b ${ctsnp} > ${outdir}/${RUN}.filt_cutsite_ctsnp.tmp

#Only keep sites on the major chromosomes
bedtools intersect -wb -a ${chrsbed} -b ${outdir}/${RUN}.filt_cutsite_ctsnp.tmp | awk '{OFS="\t"}{print $4, $5, $6, $7, $8, $9}' | sort | uniq | bgzip -c > ${outdir}/${RUN}.CpG_5mC.cov.gz
```

And run it:

```bash
for sample in $(cat HZ-samples.list); do sbatch -J ${sample} HZ_Pipeline.sh ${sample}; done 
```

And then count the sites:

```bash
mkdir counts
echo -e "Sample\tFilter\tCount" > counts/header.txt

for i in $(cat Samples.list); do 

#make a file with the first 2 columns
echo -e "${i}\n${i}\n${i}\n${i}" > counts/${i}.sample
echo -e "Cutsite\nTsSNP\nChr\nPostCoverage" > counts/${i}.field


#now for the counts
zcat ${i}.cut-sites.bismark.cov.gz | wc -l > counts/${i}.counts
cat ${i}.filt_cutsite_ctsnp.tmp | wc -l >> counts/${i}.counts
zcat ${i}.CpG_5mC.cov.gz | wc -l >> counts/${i}.counts
zcat ${i}.CpG_5mC.cov.gz | awk '{print $5+$6}' | awk '$1 > 5' | wc -l >> counts/${i}.counts

paste counts/${i}.sample counts/${i}.field counts/${i}.counts > counts/${i}.count

done
cat counts/header.txt counts/*.count > 5mC_Counts.txt
rm -rf counts
```

## ComGar.WGBS

Full pipeline:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=24:00:00

SCRATCH=tmp/$SLURM_JOB_ID
mkdir tmp
mkdir $SCRATCH

genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Genomes/Corvus_ASM73873v5
genome=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.fa
cgi=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.CGI.bed
repeats=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.REPEATS.bed
cutsites=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.CUT-SITES.bed
ctsnp=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.TS_SNPs.bed
chrsbed=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.CHRS.bed

#working directories
rawdataWGBS=/dss/dsslegfs01/pr53da/pr53da-dss-0018/rawdata/Bisulfite-Seq/WGBS
tmpdataWGBS=/dss/dsslegfs01/pr53da/pr53da-dss-0021/INBOX/2021-11-09-WGBS/tmp_rename
rawdataRRBS=/dss/dsslegfs01/pr53da/pr53da-dss-0018/rawdata/Bisulfite-Seq/RRBS
basedir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/crow_scratch/2021-07_HybridZone
outdir=${basedir}/5mC_raw/wgbs
workdir=${basedir}/scratch/wgbs

#positional arguments
RUN=$1

mkdir ${workdir}
mkdir ${outdir}

cat ${tmpdataWGBS}/${RUN}*__R1__*fastq.gz > $SCRATCH/${RUN}.R1.fastq.gz
cat ${tmpdataWGBS}/${RUN}*__R2__*fastq.gz > $SCRATCH/${RUN}.R2.fastq.gz

cd $workdir

#Trim adapters
trim_galore --fastqc -j 10 --clip_r1 9 --clip_r2 9 --three_prime_clip_R2 1 --paired --quality 20 --output_dir ${workdir} ${SCRATCH}/${RUN}.R1.fastq.gz ${SCRATCH}/${RUN}.R2.fastq.gz

#Map reads
bismark --parallel 7 --output_dir ${workdir} --genome ${genomedir} -1 ${workdir}/${RUN}.R1_val_1.fq.gz -2 ${workdir}/${RUN}.R2_val_2.fq.gz

#Deduplicate
deduplicate_bismark --paired ${workdir}/${RUN}.R1_val_1_bismark_bt2_pe.bam --output_dir ${workdir}

#Extract methylation
bismark_methylation_extractor --parallel 7 --gzip --bedgraph --ignore_3prime_r2 1 --ignore_r2 3 --buffer_size 20G --cytosine_report --genome_folder ${genomedir} --output ${workdir} ${workdir}/${RUN}.R1_val_1_bismark_bt2_pe.deduplicated.bam
```

And submit a bash script for each sample:

```bash
for sample in $(cat WGBS-samples.list); do sbatch -J ${sample} WGBS_Pipeline.sh ${sample}; done 
```

Now, subset:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=2:00:00

genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Genomes/Corvus_ASM73873v5
genome=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.fa
cgi=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.CGI.bed
repeats=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.REPEATS.bed
cutsites=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.CUT-SITES.bed
ctsnp=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.TS_SNPs.bed
chrsbed=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.CHRS.bed

#working directories
basedir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/crow_scratch/2021-07_HybridZone
outdir=${basedir}/Nov/WGBS
workdir=${basedir}/scratch/wgbs

#positional arguments
RUN=$1

#conver to temporary file 
zcat ${workdir}/${RUN}*deduplicated.bismark.cov.gz | awk '{OFS="\t"}{print $1, $2, $3, $4, $5, $6}' | sort | uniq > ${workdir}/${RUN}.tmp

#Remove any positions overlapping a transition
bedtools subtract -A -a ${workdir}/${RUN}.tmp -b ${ctsnp} > ${workdir}/${RUN}.filt_ctsnp.tmp

#Only keep sites on the major chromosomes
bedtools intersect -wb -a ${chrsbed} -b ${workdir}/${RUN}.filt_ctsnp.tmp | awk '{OFS="\t"}{print $4, $5, $6, $7, $8, $9}' | sort | uniq | bgzip -c > ${outdir}/${RUN}.CpG_5mC.cov.gz
```

And count sites:

```bash
mkdir counts
echo -e "Sample\tFilter\tCount" > counts/header.txt

for i in $(cat WGBS-samples.list); do 

#make a file with the first 2 columns
echo -e "${i}\n${i}\n${i}\n${i}" > counts/${i}.sample
echo -e "RAW\nTsSNP\nChr\nPostCoverage" > counts/${i}.field


#now for the counts
zcat ${i}.*bismark.cov.gz | wc -l > counts/${i}.counts
cat ${i}.filt_ctsnp.tmp | wc -l >> counts/${i}.counts
cat ../../5mC_raw/wgbs/${i}.CpG_5mC.cov.gz | wc -l >> counts/${i}.counts
cat ../../5mC_raw/wgbs/${i}.CpG_5mC.cov.gz | awk '{print $5+$6}' | awk '$1 > 5' | wc -l >> counts/${i}.counts

paste counts/${i}.sample counts/${i}.field counts/${i}.counts > counts/${i}.count

done
cat counts/header.txt counts/*.count > 5mC_Counts.txt
rm -rf counts
```

## Bisulfite Conversion

I will also assess bisulfite conversion efficiency by mapping our libraries to the crow mitochondrion, and our libraries with a lambda phage spike (n = 30; common garden + WGBS libraries) to the lambda phage genome. 

Common garden:

```bash
#!/bin/bash

#SBATCH -A snic2018-3-655
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 1:00:00

#genome-based files
mitogen=/crex/proj/uppstore2019047/nobackup/Genomes/Corvus_mtDNA/corvus_mitochondrion.fa
lamgen=/crex/proj/uppstore2019047/nobackup/Genomes/Lambda/LambdaPhage.fa

#Directory paths
workdir=/crex/proj/uppstore2019047/nobackup/rrbs/corvus/scratch/bs_conversion/cg
outdir=/crex/proj/uppstore2019047/nobackup/5mC_Filtered_Data/Conversion/CG
datadir=/crex/proj/uppstore2019047/nobackup/rrbs/corvus/scratch/cg

mkdir $workdir
mkdir $outdir

RUN=$1

#Map reads
bismark --parallel 6 --output_dir ${workdir} --genome ${lamgen} --prefix ${RUN}_lambda ${datadir}/${RUN}.cleaned.fastq.gz
bismark --parallel 6 --output_dir ${workdir} --genome ${mitogen} --prefix ${RUN}_mtdna ${datadir}/${RUN}.cleaned.fastq.gz
```

WGBS:

```bash
#!/bin/bash

#SBATCH -A snic2018-3-655
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 4:00:00

#genome-based files
mitogen=/crex/proj/uppstore2019047/nobackup/Genomes/Corvus_mtDNA/corvus_mitochondrion.fa
lamgen=/crex/proj/uppstore2019047/nobackup/Genomes/Lambda/LambdaPhage.fa

#Directory paths
workdir=/crex/proj/uppstore2019047/nobackup/rrbs/corvus/scratch/bs_conversion/wgbs
outdir=/crex/proj/uppstore2019047/nobackup/5mC_Filtered_Data/Conversion/WGBS
datadir=/crex/proj/uppstore2019047/nobackup/wgbs/corvus/scratch

mkdir $workdir
mkdir $outdir

RUN=$1

#Map reads
bismark --parallel 6 --output_dir ${workdir} --genome ${lamgen} --prefix ${RUN}_lambda -1 ${datadir}/${RUN}.R1.cleaned.fastq.gz -2 ${datadir}/${RUN}.R2.cleaned.fastq.gz
bismark --parallel 6 --output_dir ${workdir} --genome ${mitogen} --prefix ${RUN}_mtdna -1 ${datadir}/${RUN}.R1.cleaned.fastq.gz -2 ${datadir}/${RUN}.R2.cleaned.fastq.gz
```

And hybrid zone:

```bash
#!/bin/bash

#SBATCH -A snic2018-3-655
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 1:00:00

#genome-based files
mitogen=/crex/proj/uppstore2019047/nobackup/Genomes/Corvus_mtDNA/corvus_mitochondrion.fa

#Directory paths
workdir=/crex/proj/uppstore2019047/nobackup/rrbs/corvus/scratch/bs_conversion/hz
outdir=/crex/proj/uppstore2019047/nobackup/5mC_Filtered_Data/Conversion/HZ
datadir=/crex/proj/uppstore2019047/nobackup/rrbs/corvus/scratch/hz

mkdir $workdir
mkdir $outdir

RUN=$1

#Map reads
bismark --parallel 6 --output_dir ${workdir} --genome ${mitogen} --prefix ${RUN}_mtdna -1 ${datadir}/${RUN}.R1.cleaned.fastq.gz -2 ${datadir}/${RUN}.R2.cleaned.fastq.gz
```

And submit all 3:

```bash
for sample in $(cat /crex/proj/uppstore2019047/nobackup/rrbs/corvus/lists/CG_Samples.list); do sbatch -J ${sample}_BScg CG_conversion.sh ${sample}; done
for sample in $(cat /crex/proj/uppstore2019047/nobackup/rrbs/corvus/lists/HZ_Samples.list); do sbatch -J ${sample}_BShz HZ_conversion.sh ${sample}; done
for sample in $(cat /crex/proj/uppstore2019047/nobackup/wgbs/corvus/lists/WGBS_Samples.list); do sbatch -J ${sample}_BSwg WGBS_conversion.sh ${sample}; done
```

## M-Bias & Counts

Create m-bias plots with this:

```bash
for RUN in $(ls *M-bias.txt | sed 's/\..*//g' ); do 

sed -n '/CHG context (R1)/q;p' ${RUN}.*M-bias.txt | egrep -v 'methyl|=|context' \
| grep . | awk -v s=${RUN} '{OFS="\t"}{print s,$0,"R1","CpG"}' > ${RUN}.MBIAS1.tmp

sed -n '/CpG context (R2)/,$p' ${RUN}.*M-bias.txt | sed -n '/CHG context (R2)/q;p'  | egrep -v 'methyl|=|context' \
| grep . | awk -v s=${RUN} '{OFS="\t"}{print s,$0,"R2","CpG"}' > ${RUN}.MBIAS2.tmp

cat ${RUN}.MBIAS1.tmp ${RUN}.MBIAS2.tmp > ${RUN}.M-Bias.txt
rm ${RUN}*MBI*tmp
done
```

And within R:

```r
setwd('F:/Research/scratch/crow_hybrid_paper/2022_02/QC')

library(ggplot2)
library(gridExtra)
library(ggsci)
library(viridis)
library(dplyr)
library(ggpubr)

bias <- read.table('M_Bias.txt',header=T)
bias$Experiment <- factor(bias$Experiment,levels=c('Co.Gar.WGBS','Co.Gar.RRBS','Hyb.Zo.RRBS'))
a <- ggplot(bias,aes(x=Positions,col=Sample))+
  geom_line(aes(y=Methylated),stat="identity",show.legend=F,size=2)+ylab('Percent Methylation')+
  theme_classic(base_size=16)+facet_grid(Experiment~Read)+
  scale_color_viridis(discrete=T,option='turbo')+
  coord_cartesian(ylim=c(0,100))+
  ggtitle('Methylation Along Read Length')

a

b <- ggplot(bias,aes(x=Positions,col=Sample,y=log(Coverage)))+
  geom_line(stat="identity",show.legend=F,size=2)+ylab('log(Coverage)')+
  theme_classic(base_size=16)+facet_grid(Experiment~Read)+
  scale_color_viridis(discrete=T,option='turbo')+
  ggtitle('Coverage Along Read Length')
b

png('M-Bias.png',units='in',height=8,width=12,res = 300)
ggarrange(a,b,
          ncol=2)
dev.off()

#high HZ individuals removed
hz <- read.table('HZ_Kept.txt',header=F)
names(hz) <- 'Sample'
hzf <- subset(bias,Experiment == 'Hyb.Zo.RRBS')
hzr <- merge(hzf,hz,by='Sample')
hzb <- hzf[!grepl(paste0(hz$Sample,collapse='|'),hzf$Sample),]
hzb$Experiment <- 'Hyb.Zo.Removed'
hzr$Experiment <- 'Hyb.Zo.RRBS'
hzf <- rbind(hzb,hzr)
ab <- subset(bias,Experiment != 'Hyb.Zo.RRBS')
mb <- rbind(ab,hzf)
sampdat <- mb %>% group_by(Experiment,Sample) %>% summarize(meanM=mean(Methylated),
                                                            meanC=mean(Coverage))
samps <- ggplot(sampdat,aes(x=Sample,y=meanM,col=log(meanC)))+
  geom_point(size=5)+
  facet_grid(.~Experiment,scales='free')+
  scale_color_continuous(low='yellow',high='red')+
  theme_minimal(base_size=10)+
  theme(axis.text.x=element_text(angle=90))
samps

png('Sample_Coverage-Bias.png',units='in',height=5,width=12,res = 300)
samps
dev.off()
```

## No Cutsite or SNP Filtering

ComGarRRBS:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=48:00:00

SCRATCH=tmp/$SLURM_JOB_ID
mkdir tmp
mkdir $SCRATCH

genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Genomes/Corvus_ASM73873v5
genome=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.fa
chrsbed=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.CHRS.bed

#working directories
rawdataWGBS=/dss/dsslegfs01/pr53da/pr53da-dss-0018/rawdata/Bisulfite-Seq/WGBS
rawdataRRBS=/dss/dsslegfs01/pr53da/pr53da-dss-0018/rawdata/Bisulfite-Seq/RRBS
basedir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/bams
outdir=${basedir}/cg
workdir=${basedir}/cg

#positional arguments
RUN=$1

echo "Working on ${RUN}"

#Extract methylation
bismark_methylation_extractor --parallel 7 --ignore 2 --ignore_3prime 1 --gzip --bedgraph --buffer_size 20G --cytosine_report --genome_folder ${genomedir} --output ${outdir} ${workdir}/${RUN}.SE_trimmed_bismark_bt2.bam

```

HybZonRRBS:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=48:00:00

genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Genomes/Corvus_ASM73873v5
genome=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.fa
chrsbed=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.CHRS.bed

#working directories
rawdataWGBS=/dss/dsslegfs01/pr53da/pr53da-dss-0018/rawdata/Bisulfite-Seq/WGBS
rawdataRRBS=/dss/dsslegfs01/pr53da/pr53da-dss-0018/rawdata/Bisulfite-Seq/RRBS
basedir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/bams
outdir=${basedir}/hz
workdir=${basedir}/hz

#positional arguments
RUN=$1

echo "Working on ${RUN}"

#Extract methylation
bismark_methylation_extractor --parallel 7 --ignore 2 --ignore_3prime 1 --ignore_r2 2 --ignore_3prime_r2 1 --gzip --bedgraph --buffer_size 20G --cytosine_report --genome_folder ${genomedir} --output ${outdir} ${workdir}/${RUN}.R1_val_1_bismark_bt2_pe.bam

```

And full for ComGarWGBS:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=48:00:00

TRASH=$(realpath tmp/$SLURM_JOB_ID)
mkdir -p $TRASH

genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Genomes/Corvus_ASM73873v5
genome=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.fa
chrsbed=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.CHRS.bed

#working directories
rawdataWGBS=/dss/dsslegfs01/pr53da/pr53da-dss-0018/rawdata/Bisulfite-Seq/WGBS
rawdataRRBS=/dss/dsslegfs01/pr53da/pr53da-dss-0018/rawdata/Bisulfite-Seq/RRBS
basedir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/bams
outdir=${basedir}/wgbs
workdir=${basedir}/tmp

#positional arguments
RUN=$1

mkdir -p ${workdir}
mkdir -p ${outdir}

echo "Working on ${RUN}"

cat ${rawdataWGBS}/${RUN}*__R1__*fastq.gz > $TRASH/${RUN}.R1.fastq.gz
cat ${rawdataWGBS}/${RUN}*__R2__*fastq.gz > $TRASH/${RUN}.R2.fastq.gz

cd $workdir

#Trim adapters
trim_galore --fastqc -j 10 --clip_r1 9 --clip_r2 9 --three_prime_clip_R2 1 --paired --quality 20 --output_dir ${workdir} ${TRASH}/${RUN}.R1.fastq.gz ${TRASH}/${RUN}.R2.fastq.gz

#Map reads
bismark --parallel 7 --output_dir ${workdir} --genome ${genomedir} -1 ${workdir}/${RUN}.R1_val_1.fq.gz -2 ${workdir}/${RUN}.R2_val_2.fq.gz

#Deduplicate
deduplicate_bismark --paired ${workdir}/${RUN}.R1_val_1_bismark_bt2_pe.bam --output_dir ${outdir}

#Extract methylation
bismark_methylation_extractor --parallel 7 --gzip --bedgraph --ignore_3prime_r2 1 --ignore_r2 3 --buffer_size 20G --cytosine_report --genome_folder ${genomedir} --output ${outdir} ${outdir}/${RUN}.R1_val_1_bismark_bt2_pe.deduplicated.bam

```

Cleanup:

```bash
rm *bt2.txt.gz *bt2_pe.txt.gz *splitting* *M-bias* *bedGraph* *CpG_report* *context*
```

Now we just have `.bam` files and `.cov.gz` files for analysis. 

# 5mC DMP Analyses

## Overlap, Estimate DMPs

```R
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

```

Just submit with an ol:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_highmem
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00

Rscript 1.Merge.R $1
```

## Including ComGar.WGBS

### ComGar.RRBS

Import the bismark `cov.gz` files, filter for 10x coverage minimum and remove sites at the top 1% of the coverage distribution, calculate mean divergence between groups for factor variables, and use DSS to estimate DMPs for Age, Sex, Taxon.

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
```

### ComGar.WGBS

Import the bismark `cov.gz` files, filter for 10x coverage minimum and remove sites at the top 1% of the coverage distribution, calculate mean divergence between groups for factor variables, and use DSS to estimate DMPs for Sex, Tissue, Taxon.

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
```

### HybZon.RRBS

Import the bismark `cov.gz` files, filter for 10x coverage minimum and remove sites at the top 1% of the coverage distribution, and use DSS to estimate DMPs for Genetic Hybrid Index and Year effects. 

```R
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

```

## More DMPs on Chr18?

```R
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_202504')
.libPaths('~/mambaforge/envs/r/lib/R/library')
### Are DMPs enriched on chr18?
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

hz <- read.table('HZ_DSS.10x.txt',header=TRUE)
hz_cands <- hz %>% filter(fdrs_Chr18_Hybrid_Index.HZ < 0.01) %>% select(chr,start) %>% mutate(Experiment = 'HybZon')

cg <- read.table('CG_DSS.10x.txt',header=TRUE)
cg_cands <- cg %>% filter(fdrs_SpeciesS.CG < 0.01) %>% select(chr,start) %>% mutate(Experiment = 'ComGar')

taxcands <- rbind(hz_cands, cg_cands) %>% filter(!grepl('scaf',chr)) %>% as_tibble
write_zero <- taxcands %>% mutate(end = start, start = start -1 ) %>% select(chr,start,end,Experiment)
write.table(write_zero,file='CG-HZ_DMPsTaxon.bed',quote=F,sep='\t',row.names=F,col.names = F)

chr_ord <- taxcands %>% select(chr) %>% unique %>%  
  mutate(chrs = gsub('[chr.A]', '', chr),   # Remove 'chr', '.', and 'A'
         chrs = as.numeric(gsub('Z', '99', chrs))) %>% arrange(chrs)
taxcands$chr <- factor(taxcands$chr,levels=chr_ord$chr)

# Unadjusted 
count_plot <- taxcands %>% 
  ggplot(aes(y=chr,fill=Experiment))+
  geom_bar()+
  xlab('Count of Taxon DMPs')+ylab('')+
  scale_fill_manual(values=c('#33638d','#b8de29'))+
  theme_bw(base_size=8)
count_plot

taxcounts <- taxcands %>% 
  count(chr,Experiment)


## correct for counts by the number of DMPs assayed by chromosome
cg_counts <- cg %>% count(chr) %>% filter(!grepl('MT|scaf',chr)) %>% arrange(desc(n)) %>% mutate(Experiment = 'ComGar')
hz_counts <- hz %>% count(chr) %>% filter(!grepl('MT|scaf',chr)) %>% arrange(desc(n)) %>% mutate(Experiment = 'HybZon')
xp_counts <- rbind(cg_counts,hz_counts)
all_chrs <- unique(xp_counts$chr)

# init df
binomial_results <- data.frame()

# loop xp & chr
for (exp in unique(taxcounts$Experiment)) {
  exp_data <- taxcounts %>% 
    filter(Experiment == exp) %>% 
    complete(chr = all_chrs, Experiment, fill = list(n = 0)) %>% left_join(.,xp_counts %>% dplyr::rename(chr_total = n)) %>% 
    mutate(gw_total = sum(chr_total),
           expected_prop = chr_total / gw_total)
  exp_total_dmps <- sum(exp_data$n)
  
  for (chr in unique(exp_data$chr)) {
    chr_dmps <- exp_data %>% filter(chr == !!chr) %>% pull(n)
    chr_expected_prop <- exp_data %>% filter(chr == !!chr) %>% pull(expected_prop)
    
    # bin test
    bin_test <- binom.test(chr_dmps, exp_total_dmps, chr_expected_prop)
    
    result_row <- data.frame(
      Experiment = exp,
      Chromosome = chr,
      chr_DMPs = chr_dmps,
      Total_DMPs = exp_total_dmps,
      Expected_proportion = chr_expected_prop,
      p_value = bin_test$p.value,
      lower_CI = bin_test$conf.int[1],
      upper_CI = bin_test$conf.int[2],
      observed_proportion = bin_test$estimate
    )
    
    binomial_results <- bind_rows(binomial_results, result_row)
  }
}

binomial_results <- binomial_results %>%
  mutate(significant = (lower_CI > Expected_proportion) )
binomial_results$Chromosome <- factor(binomial_results$Chromosome,levels=chr_ord$chr)
write.table(binomial_results,file='~/symlinks/hzone/Methylation_Analyses/figures/20250418_Binomial_Results.txt',quote=F,sep='\t',row.names=F)

binom_plot <- binomial_results %>%
  ggplot(aes(y = Chromosome, color = Experiment)) +
  
  # Observed proportion and 95% CI
  geom_point(aes(x = observed_proportion), size = 2) +
  geom_errorbarh(aes(xmin = lower_CI, xmax = upper_CI), height = 0.2) +
  
  # Expected proportion
  geom_point(aes(x = Expected_proportion), shape = 4, size = 2, stroke = 0.5, color = "black") +
  
  # asterisk for significant deviations
  geom_text(data = . %>% filter(significant),
            aes(x = max(upper_CI, Expected_proportion) + 0.005, label = "*"),  # nudge slightly right
            color = "black", size = 3) +
  
  facet_wrap(~ Experiment, scales = "free_y") +
  scale_color_manual(values=c('#33638d','#b8de29'))+
  labs(
    x = "Proportion of DMPs",
    y = "Chromosome"
  ) +
  theme_bw(base_size = 8)
binom_plot

ggsave(ggarrange(count_plot,binom_plot,common.legend = TRUE),file='~/symlinks/hzone/Methylation_Analyses/figures/20250418_DMP_Counts.pdf',height=5,width=6)


```

## chrZ Saturation

```R
.libPaths('~/mambaforge/envs/tidymodels/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_202504/RFS_202504/')
library(tidyverse)

all <- read.table('../Full-HZ-CG.DSS.BP-10x.Classified-Annotated-SNPless_20250410.txt',header=T)

chr_z_sex <- all %>% filter(chr == 'chrZ' & Class == 'Sex') %>% nrow
total_sex <- all %>% filter(Class == 'Sex') %>% nrow
chr_z_CpGs <- all %>% filter(chr == 'chrZ' & DMP.CG != 'MISSING') %>% nrow
total_CpGs <- all %>% filter(DMP.CG != 'MISSING') %>% nrow

bin_test <- binom.test(x = chr_z_sex, n = total_sex, p = chr_z_CpGs / total_CpGs)

result_row <- data.frame(
  chr_Z_sex_DMPs = chr_z_sex,
  total_sex_DMPs = total_sex, 
  chr_Z_CpGs = chr_z_CpGs,
  total_CpGs = total_CpGs,
  Expected_proportion = chr_z_CpGs / total_CpGs,
  p_value = bin_test$p.value,
  lower_CI = bin_test$conf.int[1],
  upper_CI = bin_test$conf.int[2],
  observed_proportion = bin_test$estimate)
result_row  

# chr_Z_sex_DMPs total_sex_DMPs chr_Z_CpGs total_CpGs Expected_proportion      p_value   lower_CI
# probability of success            163           2768      27890    1047318          0.02662993 5.741463e-20 0.05040697
# upper_CI observed_proportion
# probability of success 0.06831449          0.05888728

```



## DMRs

```bash
 grep 'Chr18_' HZ_DMRs-10x.txt | bedtools sort -i - > taxonfocus/HZ_DMRs.taxon.bed
 grep 'Species' CG_DMRs-10x.txt | bedtools sort -i - > taxonfocus/CG_DMRs.taxon.bed
```

in R:

```R
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_202504')
.libPaths('~/mambaforge/envs/r/lib/R/library')
### DMRs
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

dmrs <- read_tsv('DMRs.txt') %>% filter(!grepl('scaf|MT',chr))
dmrs %>% count(var)
dmrs <- dmrs %>% mutate(Experiment = case_when(
  var == 'Chr18_Hybrid_Index' ~ 'HybZon', 
  grepl('Year',var) ~ 'HybZon', 
  TRUE ~ 'ComGar'
)) 

taxdmr <- dmrs %>% filter(var == 'Chr18_Hybrid_Index' | var == 'SpeciesS')
dmpplot <- taxdmr %>% ggplot(aes(y=chr,fill=Experiment))+
  geom_bar()+
  scale_fill_manual(values=brewer.pal(3,'Set2'))+
  theme_bw(base_size=8)

ggsave('~/symlinks/hzone/Methylation_Analyses/figures/20250415_DMR_Counts.pdf',dmpplot, height=3,width=4)

# Intersected closest
int <- read_tsv('taxonfocus/Intersect.txt',col_names = F)
int %>% arrange(X15)

# # A tibble: 46 × 15
# X1          X2       X3    X4    X5    X6 X7                 X8         X9    X10   X11   X12    X13 X14     X15
# <chr>    <dbl>    <dbl> <dbl> <dbl> <dbl> <chr>              <chr>   <dbl>  <dbl> <dbl> <dbl>  <dbl> <chr> <dbl>
#   1 chr11  1431927  1432029   103    18 -59.9 Chr18_Hybrid_Index chr11  1.43e6 1.43e6   412    10 -48.2  Spec…     0
# 2 chr1A   374930   375111   182    30 -88.8 Chr18_Hybrid_Index chr1A  3.75e5 3.75e5   169    16 -39.6  Spec…     0
# 3 chr28  4327237  4327301    65     7  15.2 Chr18_Hybrid_Index chr28  4.33e6 4.33e6   195    21 130.   Spec…     0
# 4 chr26  6273580  6273635    56     8 -25.4 Chr18_Hybrid_Index chr26  6.28e6 6.28e6    71     4  11.1  Spec…  1681
# 5 chr4A  3586190  3586242    53     5 -11.8 Chr18_Hybrid_Index chr4A  3.58e6 3.58e6    82     9  26.2  Spec…  2683
# 6 chr3  18729087 18729150    64     8 -20.7 Chr18_Hybrid_Index chr3   1.87e7 1.87e7   483     5  22.8  Spec…  7577
```

| HybZon                | ComGar |         |         |        |             |           |       |          |          |        |             |           |              |               |
| --------------------- | ------ | ------- | ------- | ------ | ----------- | --------- | ----- | -------- | -------- | ------ | ----------- | --------- | ------------ | ------------- |
| coordinates           | chr    | start   | end     | Length | Number CpGs | Area Stat | chr   | start    | end      | Length | Number CpGs | Area Stat | Gene         | Distance      |
| chr11:1431927-1432029 | chr11  | 1431927 | 1432029 | 103    | 18          | -59.9     | chr11 | 1.43E+06 | 1.43E+06 | 412    | 10          | -48.2     | CBFA2T3      | Within        |
| chr1A:374930-375111   | chr1A  | 374930  | 375111  | 182    | 30          | -88.8     | chr1A | 3.75E+05 | 3.75E+05 | 169    | 16          | -39.6     | DENND6B      | 2 Kb upstream |
| chr28:4327237-4327301 | chr28  | 4327237 | 4327301 | 65     | 7           | 15.2      | chr28 | 4.33E+06 | 4.33E+06 | 195    | 21          | 130       | LOC104697407 | Within        |



## Consensus Classification

Require an FDR of 1%, a primary effect divergence of 25%, non-significance for other effects. Conserved sites exhibit less than a 10% difference in methylation range. 

```R
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_202504/')
.libPaths('~/mambaforge/envs/randomforest/lib/R/library')
library(tidyverse)
library(data.table)
library(viridis)

# Thresholds: FDR corrected pvalue, divergence between groups for factors (e.g. 25%), min% for conserved
pval <- 0.01
div <- 0.25
dd <- 0.1

#### CG
datsCG <- read.table('CG_DSS.10x.txt',header=TRUE)
datsCG <- datsCG %>% dplyr::rename(fdrs_Species.CG = fdrs_SpeciesS.CG, fdrs_Sex.CG = fdrs_SexM.CG)

# With interacting effects, exclude variables significant for other variables
cg <- datsCG %>% mutate(DMP.CG = ifelse(fdrs_Species.CG < pval & abs(divergence_Species.CG) > div &
                                          fdrs_Sex.CG > pval & fdrs_StageCHK.CG > pval & fdrs_StageYRL.CG > pval, 'Taxon',  #species
                                        ifelse(fdrs_Sex.CG < pval & abs(divergence_Sex.CG) > div &
                                                 fdrs_Species.CG > pval & fdrs_StageCHK.CG > pval & fdrs_StageYRL.CG > pval, 'Sex', #sex
                                               ifelse((fdrs_StageCHK.CG < pval | fdrs_StageYRL.CG < pval) & (abs(divergence_StageCHK.CG) > div | abs(divergence_StageYRL.CG) > div) &
                                                        (fdrs_Sex.CG > pval & fdrs_Species.CG > pval),'Age',  #stage
                                                      ifelse(Divergence.CG < dd,  'Conserved','Unknown')))))

cg %>% dplyr::count(DMP.CG)


#merge with HZ data..
hz <- read.table('HZ_DSS.10x.txt',header=T)
hz = hz %>% mutate(DMP.HZ = ifelse(fdrs_Chr18_Hybrid_Index.HZ < pval & fdrs_Year2013.HZ > pval & fdrs_Year2014.HZ > pval, 'Taxon',
                                   ifelse(Divergence.HZ < dd, 'Conserved','Unknown')))
hz %>% dplyr::count(DMP.HZ)
dmr2 <- merge(cg,hz,by=c('site','chr','start'),all=TRUE)

#count NAs by experiment
nac = dmr2 %>% select(site,DMP.CG,DMP.HZ)
nac$NAs = rowSums(is.na(nac))
nac %>% dplyr::count(NAs)

#replace NA with missing
dmr3 = dmr2 %>% replace_na(list(DMP.CG ='MISSING',DMP.HZ = 'MISSING'))

# Assign sites based on conflicts between experiments
dmr4 <- dmr3 %>% 
  mutate(
    Class = case_when(
      # Taxon classification
      (DMP.CG == "Taxon" & DMP.HZ %in% c("Taxon", "Unknown", "Conserved", "MISSING")) ~ "Taxon",
      (DMP.HZ == "Taxon" & DMP.CG %in% c("Taxon", "Unknown", "Conserved", "MISSING")) ~ "Taxon",
      
      # Age classification
      (DMP.CG == "Age" & DMP.HZ %in% c("Unknown", "Conserved", "MISSING"))  ~ "Age",
      
      # Sex classification
      (DMP.CG == "Sex" & DMP.HZ %in% c("Unknown", "Conserved", "MISSING")) ~ "Sex",
      
      # Conserved classification
      (DMP.CG == "Conserved" & DMP.HZ == "Conserved") ~ "Conserved",
      
      # Conflict
      TRUE ~ "Unknown"
    ),
    # Count experiments supporting each classification
    SupportCount = case_when(
      Class == "Taxon" ~ rowSums(cbind(DMP.CG == "Taxon", DMP.HZ == "Taxon")),
      Class == "Age" ~ as.numeric(DMP.CG == "Age"),
      Class == "Sex" ~ as.numeric(DMP.CG == "Sex"),
      Class == "Conserved" ~ 2, # All three experiments must agree
      TRUE ~ 0 # For conflicts or unknown classifications
    )
  )

dmr4 <- dmr4 %>% mutate(Focal_Region = ifelse(chr == 'chr18' & start > 8.07e6 & start < 10.07e6,'Focal','Undiverged'))
dmr4 %>% count(SupportCount,Class)
dmr4 %>% filter(Class == 'Taxon' & SupportCount >= 2 & DMP.HZ == 'Taxon')
dmr4 %>% filter(Class == 'Taxon' & SupportCount >= 2 & DMP.HZ == 'Taxon') %>% select(matches('site|Focal|fdrs')) %>% select(matches('site|Focal|Species|Hybrid'))


dmr4 %>% filter(chr == 'chr18' & start %in% c(10000910, 9005914, 9005915)) %>% select(site,contains('.HZ'))
dmr5 = dmr4 
dmr5$Mean5mC.CG = rowMeans(dmr5[grepl('^D_.*CG|^S_.*CG',names(dmr5))],na.rm=TRUE)
dmr5$Mean5mC.HZ = rowMeans(dmr5[grepl('^D_.*HZ|^S_.*HZ|^H_.*HZ',names(dmr5))],na.rm=TRUE)                                                     

write.table(dmr5,'Full-HZ-CG.DSS.BP-10x.Classified-20250410.txt',quote=F,sep='\t',row.names=FALSE)

```

### Including ComGar.WGBS

```R
.libPaths('~/miniconda3/envs/randomforest/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_Sept/')
library(tidyverse)
library(data.table)
library(viridis)

# Thresholds: FDR corrected pvalue, divergence between groups for factors (e.g. 25%), min% for conserved
pval <- 0.01
div <- 0.25
dd <- 0.1

#### CG
datsCG <- read.table('CG.DSS.BP-10x.txt',header=TRUE)
datsCG <- datsCG %>% dplyr::rename(fdrs_Species.cg = fdrs_SpeciesS.cg, fdrs_Sex.cg = fdrs_SexM.cg)

# With interacting effects, exclude variables significant for other variables
CG <- datsCG %>% mutate(DMP = ifelse(fdrs_Species.cg < pval & abs(divergence_Species.cg) > div &
                                       fdrs_Sex.cg > pval & fdrs_StageCHK.cg > pval & fdrs_StageYRL.cg > pval, 'Taxon',  #species
                                     ifelse(fdrs_Sex.cg < pval & abs(divergence_Sex.cg) > div &
                                              fdrs_Species.cg > pval & fdrs_StageCHK.cg > pval & fdrs_StageYRL.cg > pval, 'Sex', #sex
                                            ifelse((fdrs_StageCHK.cg < pval | fdrs_StageYRL.cg < pval) & (abs(divergence_StageCHK.cg) > div | abs(divergence_StageYRL.cg) > div) &
                                                     (fdrs_Sex.cg > pval & fdrs_Species.cg > pval),'Age',  #stage
                                                   ifelse(Divergence.cg < dd,  'Conserved','Unknown')))))

CG %>% dplyr::count(DMP)

#### WGBS
datsWGBS <- read.table('WGBS.DSS.BP-10x.txt',header=TRUE)

#with interacting effects, exclude variables significant for other variables
WGBS <- datsWGBS %>% mutate(DMP = ifelse(fdrs_Species.wgbs < pval & abs(divergence_Species.wgbs) > div &
                                           fdrs_Sex.wgbs > pval & fdrs_TissueM.wgbs > pval & fdrs_TissueL.wgbs > pval, 'Taxon',  #species
                                         ifelse(fdrs_Sex.wgbs < pval & abs(divergence_Sex.wgbs) > div &
                                                  fdrs_Species.wgbs > pval & fdrs_TissueM.wgbs > pval & fdrs_TissueL.wgbs > pval, 'Sex', #sex
                                                ifelse((fdrs_TissueM.wgbs < pval | fdrs_TissueL.wgbs < pval) & (abs(divergence_TissueM.wgbs) > div | abs(divergence_TissueL.wgbs) > div) &
                                                         (fdrs_Sex.wgbs > pval & fdrs_Species.wgbs > pval),'Tissue',  #tissue
                                                       ifelse(Divergence.wgbs < dd,  'Conserved','Unknown')))))


WGBS %>% dplyr::count(DMP)

#merged
dmr1 <- merge(WGBS,CG,by=c('site','chr','start'),
              suffixes=c('.wgbs','.cg'),all=TRUE)

#merge with HZ data..
hz <- read.table('HZ.DSS.BP-Parentals-YearChr18HI_10x.txt',header=T)
hz = hz %>% mutate(DMP.hz = ifelse(fdrs_Chr18_Hybrid_Index.hz < pval & fdrs_Year2013.hz > pval & fdrs_Year2014.hz > pval, 'Taxon',
                                   ifelse(Divergence.hz < dd, 'Conserved','Unknown')))
hz %>% dplyr::count(DMP.hz)
dmr2 <- merge(dmr1,hz,by=c('site','chr','start'),all=TRUE)

#count NAs by experiment
nac = dmr2 %>% select(site,DMP.wgbs,DMP.cg,DMP.hz)
nac$NAs = rowSums(is.na(nac))
nac %>% dplyr::count(NAs)

#replace NA with missing
dmr3 = dmr2 %>% replace_na(list(DMP.wgbs = 'MISSING',DMP.cg ='MISSING',DMP.hz = 'MISSING'))

# Assign sites based on conflicts between experiments
dmr4 <- dmr3 %>% 
  mutate(
    Class = case_when(
      # Taxon classification
      (DMP.cg == "Taxon" & DMP.hz %in% c("Taxon", "Unknown", "Conserved", "MISSING") & DMP.wgbs %in% c("Taxon", "Unknown", "Conserved", "MISSING")) ~ "Taxon",
      (DMP.hz == "Taxon" & DMP.cg %in% c("Taxon", "Unknown", "Conserved", "MISSING") & DMP.wgbs %in% c("Taxon", "Unknown", "Conserved", "MISSING")) ~ "Taxon",
      (DMP.wgbs == "Taxon" & DMP.hz %in% c("Taxon", "Unknown", "Conserved", "MISSING") & DMP.cg %in% c("Taxon", "Unknown", "Conserved", "MISSING")) ~ "Taxon",
      
      # Tissue classification
      (DMP.wgbs == "Tissue" & DMP.hz %in% c("Unknown", "Conserved", "MISSING") & DMP.cg %in% c("Unknown", "Conserved", "MISSING")) ~ "Tissue",
      
      # Age classification
      (DMP.cg == "Age" & DMP.hz %in% c("Unknown", "Conserved", "MISSING") & DMP.wgbs %in% c("Unknown", "Conserved", "MISSING")) ~ "Age",
      
      # Sex classification
      (DMP.cg == "Sex" & DMP.hz %in% c("Unknown", "Conserved", "MISSING") & DMP.wgbs %in% c("Sex", "Unknown", "Conserved", "MISSING")) ~ "Sex",
      (DMP.wgbs == "Sex" & DMP.hz %in% c("Unknown", "Conserved", "MISSING") & DMP.cg %in% c("Sex", "Unknown", "Conserved", "MISSING")) ~ "Sex",
      
      # Conserved classification
      (DMP.cg == "Conserved" & DMP.hz == "Conserved") ~ "Conserved",
      
      # Conflict
      TRUE ~ "Unknown"
    ),
    # Count experiments supporting each classification
    SupportCount = case_when(
      Class == "Taxon" ~ rowSums(cbind(DMP.cg == "Taxon", DMP.wgbs == "Taxon", DMP.hz == "Taxon")),
      Class == "Tissue" ~ as.integer(DMP.wgbs == "Tissue"),
      Class == "Age" ~ as.integer(DMP.cg == "Age"),
      Class == "Sex" ~ rowSums(cbind(DMP.cg == "Sex", DMP.wgbs == "Sex")),
      Class == "Conserved" ~ 3, # All three experiments must agree
      TRUE ~ 0 # For conflicts or unknown classifications
    )
  )

dmr4 <- dmr4 %>% mutate(Focal_Region = ifelse(chr == 'chr18' & start > 8.07e6 & start < 10.07e6,'Focal','Undiverged'))
dmr4 %>% count(SupportCount,Class)
dmr4 %>% filter(Class == 'Taxon' & SupportCount >= 2 & DMP.hz == 'Taxon')
dmr4 %>% filter(Class == 'Taxon' & SupportCount >= 2 & DMP.hz == 'Taxon') %>% select(matches('site|Focal|fdrs')) %>% select(matches('site|Focal|Species|Hybrid'))


dmr4 %>% filter(chr == 'chr18' & start %in% c(10000910, 9005914, 9005915)) %>% select(site,contains('.hz'))
dmr5 = dmr4 
dmr5$Mean5mC.cg = rowMeans(dmr5[grepl('^D_.*cg|^S_.*cg',names(dmr5))],na.rm=TRUE)
dmr5$Mean5mC.wgbs = rowMeans(dmr5[grepl('^D_.*wgbs|^S_.*wgbs',names(dmr5))],na.rm=TRUE)
dmr5$Mean5mC.hz = rowMeans(dmr5[grepl('^D_.*hz|^S_.*hz|^H_.*hz',names(dmr5))],na.rm=TRUE)                                                     

write.table(dmr5,'Full.DSS.BP-10x.Classified-20250107.txt',quote=F,sep='\t',row.names=FALSE)

```
## Merge with Popgen Data

And now intersect with our REGION, DXY, FST, PI, data: 

```bash
rgndir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Regions
gc=$rgndir/GCF_000738735.5_ASM73873v5_genomic.GC-Content.500bp.bed
rgn=$rgndir/GCF_000738735.5_ASM73873v5_genomic.GenomicFeaturesCDS.bed
popdat=$rgndir/GCF_000738735.5_ASM73873v5_genomic.PopGen2023JAN.bed

#create sites file 
sed '1d' Full-HZ-CG.DSS.BP-10x.Classified-20250410.txt | awk '{OFS="\t"}{print $2, $3, $3, $1}' > Sites.bed

#create master overlapped genetic file 
bedtools intersect -wao -a Sites.bed -b $rgn | bedtools intersect -wao -a - -b $gc | bedtools intersect -header -wao -a - -b $popdat | awk '{OFS="\t"}{print $1, $16, $17, $4, $2, $8, $13, $18, $19, $20, $21, $22}' > All-Overlapped.bed
```

and re-merge with the full dataset

```r
library(tidyverse)
options(scipen=999)
init <- read.table('Full-HZ-CG.DSS.BP-10x.Classified-20250410.txt',header=T)
anno <- read.table('All-Overlapped.bed',header=F)
names(anno) <- c('chr','GenStart','GenEnd','site','start', 'Region','GC','HapD','FuLiD','TajimaD','Dxy','FST')
all <- full_join(init,anno)
all$Region <- gsub('\\.','Missing',all$Region)
nrow(all)
all <- all %>%
  mutate(end = start,
         start = start - 1) %>%
  relocate(chr, start, end, .before = 1)
write.table(all,'Full-HZ-CG.DSS.BP-10x.Classified-Annotated_20250410.txt',quote=F,sep='\t',row.names=F)
```

### Including ComGar.WGBS

Merge with our REGION, DXY, FST, PI, data: 

```bash
rgndir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/crow_scratch/2021-07_HybridZone/Nov/Regions
gc=$rgndir/GCF_000738735.5_ASM73873v5_genomic.GC-Content.500bp.bed
rgn=$rgndir/GCF_000738735.5_ASM73873v5_genomic.GenomicFeaturesCDS.bed
popdat=$rgndir/GCF_000738735.5_ASM73873v5_genomic.PopGen2023JAN.bed

#create sites file 
sed '1d' Full.DSS.BP-10x.Classified-20250107.txt | awk '{OFS="\t"}{print $2, $3, $3, $1}' > Sites_Year.bed

#create master overlapped genetic file 
bedtools intersect -wao -a Sites_Year.bed -b $rgn | bedtools intersect -wao -a - -b $gc | bedtools intersect -header -wao -a - -b $popdat | awk '{OFS="\t"}{print $1, $16, $17, $4, $2, $8, $13, $18, $19, $20, $21, $22}' > All-Overlapped-Year.bed
```

and re-merge with the full dataset

```r
library(tidyverse)
options(scipen=999)
init <- read.table('Full.DSS.BP-10x.Classified-23AUG30_YEAR.txt',header=T)
anno <- read.table('All-Overlapped-Year.bed',header=F)
names(anno) <- c('chr','GenStart','GenEnd','site','start', 'Region','GC','HapD','FuLiD','TajimaD','Dxy','FST')
all <- full_join(init,anno)
all$Region <- gsub('\\.','Missing',all$Region)
nrow(all)
all = all %>% dplyr::rename(Class = DMP_Strict)
write.table(all,'Full.DSS.BP-10x.Classified-Annotated-YEAR_23AUG30.txt',quote=F,sep='\t',row.names=F)
```

Afterwards, exclude any SNPs:

```bash
head -n 1 Full-HZ-CG.DSS.BP-10x.Classified-Annotated_20250410.txt > header.txt
sed '1d' Full-HZ-CG.DSS.BP-10x.Classified-Annotated_20250410.txt | bedtools subtract -a - -b GCF_000738735.5_ASM73873v5_genomic.AllSNPs.sorted.bed > no_snp.bed
cat header.txt no_snp.bed > Full-HZ-CG.DSS.BP-10x.Classified-Annotated-SNPless_20250410.txt
```



## UpSet Plot

```bash
.libPaths('~/mambaforge/envs/r/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_202504')
library(tidyverse)
library(UpSetR)
options(scipen=999)
all = read_tsv('Full-HZ-CG.DSS.BP-10x.Classified-Annotated-SNPless_20250410.txt')

ap = all %>% select(site,DMP.HZ,DMP.CG)
ap = ap %>% mutate(across(starts_with("DMP"), ~ ifelse(. == "MISSING", "Unknown", .)))
names(ap) = c('site','HybZon','ComGar')

#only include sites which at least 1 taxon occurrence 
taxon_data <- ap %>%
  filter(HybZon == "Taxon" | ComGar == "Taxon")

# binary presence/absence matrix 
binary_matrix <- taxon_data %>%
  mutate(HybZon = ifelse(HybZon == "Taxon", 1, 0),
         ComGar = ifelse(ComGar == "Taxon", 1, 0)) %>%
  select(HybZon, ComGar)

# sum up the binary flags to get counts for each combination
aggregated_counts <- binary_matrix %>%
  group_by(HybZon, ComGar) %>%
  summarise(count = n(), .groups = 'drop')

# 'exression' data format for upsetR 
upset_data <- c(
  HybZon = aggregated_counts$count[aggregated_counts$HybZon == 1 & aggregated_counts$ComGar == 0],
  ComGar = aggregated_counts$count[aggregated_counts$HybZon == 0 & aggregated_counts$ComGar == 1],
  "HybZon&ComGar" = aggregated_counts$count[aggregated_counts$HybZon == 1 & aggregated_counts$ComGar == 1]
)

up = upset(fromExpression(upset_data), order.by = "freq")
up

pdf('~/symlinks/hzone/Methylation_Analyses/figures/20250416_Upset_Taxon_Plot.pdf',height=4,width=6)
up
dev.off()

```

### Including ComGar.WGBS

```bash
.libPaths('~/mambaforge/envs/r/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_Sept')
library(tidyverse)
library(UpSetR)
options(scipen=999)
all = read_tsv('Full.DSS.BP-10x.Classified-Annotated_20250107.txt')

ap = all %>% select(site,DMP.wgbs,DMP.hz,DMP.cg)
ap = ap %>% mutate(across(starts_with("DMP"), ~ ifelse(. == "MISSING", "Unknown", .)))
names(ap) = c('site','ComGar.WGBS','HybZon.RRBS','ComGar.RRBS')

#only include sites which at least 1 taxon occurrence 
taxon_data <- ap %>%
  filter(ComGar.WGBS == "Taxon" | HybZon.RRBS == "Taxon" | ComGar.RRBS == "Taxon")

# binary presence/absence matrix 
binary_matrix <- taxon_data %>%
  mutate(ComGar_WGBS = ifelse(ComGar.WGBS == "Taxon", 1, 0),
         HybZon_RRBS = ifelse(HybZon.RRBS == "Taxon", 1, 0),
         ComGar_RRBS = ifelse(ComGar.RRBS == "Taxon", 1, 0)) %>%
  select(ComGar_WGBS, HybZon_RRBS, ComGar_RRBS)

# sum up the binary flags to get counts for each combination
aggregated_counts <- binary_matrix %>%
  group_by(ComGar_WGBS, HybZon_RRBS, ComGar_RRBS) %>%
  summarise(count = n(), .groups = 'drop')

# 'exression' data format for upsetR 
upset_data <- c(
  ComGar_WGBS = aggregated_counts$count[aggregated_counts$ComGar_WGBS == 1 & aggregated_counts$HybZon_RRBS == 0 & aggregated_counts$ComGar_RRBS == 0],
  HybZon_RRBS = aggregated_counts$count[aggregated_counts$ComGar_WGBS == 0 & aggregated_counts$HybZon_RRBS == 1 & aggregated_counts$ComGar_RRBS == 0],
  ComGar_RRBS = aggregated_counts$count[aggregated_counts$ComGar_WGBS == 0 & aggregated_counts$HybZon_RRBS == 0 & aggregated_counts$ComGar_RRBS == 1],
  "ComGar_WGBS&HybZon_RRBS" = aggregated_counts$count[aggregated_counts$ComGar_WGBS == 1 & aggregated_counts$HybZon_RRBS == 1 & aggregated_counts$ComGar_RRBS == 0],
  "ComGar_WGBS&ComGar_RRBS" = aggregated_counts$count[aggregated_counts$ComGar_WGBS == 1 & aggregated_counts$HybZon_RRBS == 0 & aggregated_counts$ComGar_RRBS == 1],
  "HybZon_RRBS&ComGar_RRBS" = aggregated_counts$count[aggregated_counts$ComGar_WGBS == 0 & aggregated_counts$HybZon_RRBS == 1 & aggregated_counts$ComGar_RRBS == 1],
  "All" = aggregated_counts$count[aggregated_counts$ComGar_WGBS == 1 & aggregated_counts$HybZon_RRBS == 1 & aggregated_counts$ComGar_RRBS == 1]
)

up = upset(fromExpression(upset_data), order.by = "freq")
up

pdf('20250107_Upset_Taxon_Plot.pdf',height=4,width=7)
up
dev.off()

```

And for all variables, streamlined for tons of predictors:

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
library(UpSetR)
options(scipen=999)

##### Import 5mC data ####
class <- read.table('Full.DSS.BP-10x.Classified-Annotated_20250107.txt',header=TRUE)

# Summaries
class %>% count(DMP.hz,Focal_Region) %>% filter(DMP.hz =='Taxon')
class %>% count(DMP.cg,Focal_Region) %>% filter(DMP.cg =='Taxon')
class %>% count(DMP.wgbs,Focal_Region) %>% filter(DMP.wgbs =='Taxon')

#Sex
class %>% count(Class,chr) %>% filter(Class =='Sex') %>% summarize(sum=sum(n))

#Overall
class %>% count(Class)

# Grab dmps 
ap = class %>% select(site,DMP.wgbs,DMP.hz,DMP.cg)
ap = ap %>% mutate(across(starts_with("DMP"), ~ ifelse(. == "MISSING", "Unknown", .)))
names(ap) = c('site','ComGar.WGBS','HybZon.RRBS','ComGar.RRBS')
ap %>% count(ComGar.WGBS,HybZon.RRBS,ComGar.RRBS)

# Grab the levesl which will be our xp_var IDs
all_levels <- ap %>%
  select(-site) %>%  # Exclude 'site' column if not relevant
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Level") %>%
  distinct(Variable, Level)  # Get all unique variable-level combinations

# Binarize all xp * variable combinations
binary_matrix <- ap %>%
  select(-site) %>%  # Exclude 'site' if it's not a factor for analysis
  mutate(row_id = row_number()) %>%  # Add a unique ID for each row
  pivot_longer(cols = -row_id, names_to = "Variable", values_to = "Level") %>% 
  left_join(all_levels, by = c("Variable", "Level")) %>%
  mutate(Presence = 1) %>%
  pivot_wider(names_from = c(Variable, Level), values_from = Presence, values_fill = 0) %>%
  select(-row_id)  # Remove row ID after creating the binary matrix

# Extract the coutns corresponding to each xp * variable combination
aggregated_counts <- binary_matrix %>%
  group_by(across(everything())) %>%
  summarise(count = n(), .groups = 'drop')

# Extract the column names for the counts, this is eg. xp1_var1&xp1_var2
new_df <- aggregated_counts %>%
  rowwise() %>%
  mutate(new_column = paste(names(df)[which(c_across(everything()) == 1)], collapse = "&")) %>%
  select(new_column, count) %>%
  ungroup() %>%
  filter(new_column != "")  # Remove rows with empty new_column

# Transpose into the upset format
counts <- data.frame(t(new_df$count))
names(counts) <- new_df$new_column

up <- upset(fromExpression(counts), order.by = "freq",nsets = 13)

pdf('20250107_Upset_AllVariables_Plot.pdf',height=7,width=8)
up
dev.off()

```

## Distance to SNPs

```R
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_202504/')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(data.table)
library(viridis)
library(meRo)
library(rstatix)
library(ggpubr)

pval <- 0.01

#### CG
datsCG <- read.table('CG_DSS.10x.txt',header=TRUE)
datsCG <- datsCG %>% dplyr::rename(fdrs_Species.CG = fdrs_SpeciesS.CG, fdrs_Sex.CG = fdrs_SexM.CG)

# With interacting effects, exclude variables significant for other variables
cg <- datsCG %>% mutate(DMP.CG = ifelse(fdrs_Species.CG < pval & fdrs_Sex.CG > pval & fdrs_StageCHK.CG > pval & fdrs_StageYRL.CG > pval, 'Taxon',  #species
                                        ifelse(fdrs_Sex.CG < pval & fdrs_Species.CG > pval & fdrs_StageCHK.CG > pval & fdrs_StageYRL.CG > pval, 'Sex', #sex
                                               ifelse((fdrs_StageCHK.CG < pval | fdrs_StageYRL.CG < pval) & (fdrs_Sex.CG > pval & fdrs_Species.CG > pval),'Age',
                                                      'Unknown'))))

cg %>% dplyr::count(DMP.CG)


#merge with HZ data..
hz <- read.table('HZ_DSS.10x.txt',header=T)
hz = hz %>% mutate(DMP.HZ = ifelse(fdrs_Chr18_Hybrid_Index.HZ < pval & fdrs_Year2013.HZ > pval & fdrs_Year2014.HZ > pval, 'Taxon',
                                   ifelse((fdrs_Year2014.HZ < pval | fdrs_Year2013.HZ < pval) & fdrs_Chr18_Hybrid_Index.HZ,'Year',
                                          'Unknown')))
hz %>% dplyr::count(DMP.HZ)
dmps <- rbind(cg %>% mutate(Experiment = 'ComGar',end=start,start=start-1) %>% select(chr,start,end,DMP=DMP.CG,Experiment),
      hz %>% mutate(Experiment = 'HybZon',end=start,start=start-1) %>% select(chr,start,end,DMP=DMP.HZ,Experiment)) %>% 
  filter(DMP != 'Unknown') %>% arrange(chr,start) %>% filter(!grepl('scaf',chr))
write.table(dmps,file='CG-HZ_DMPsTaxon.sorted.bed',quote=F,sep='\t',row.names=F,col.names=F)
# bedtools closest -a CG-HZ_DMPsTaxon.sorted.bed -b GCF_000738735.5_ASM73873v5_genomic.AllSNPs.sorted.bed -d > CG-HZ_DMPsTaxon.sorted-SNPDistance.bed

dist <- read_tsv('CG-HZ_DMPsTaxon.sorted-SNPDistance.bed',col_names = F)
names(dist) <- c('chr','start','end','DMP','Experiment','d1','SNPa','SNPb','Alt','Ref','Distance')

# Exclude these as they intersect directly and the 5mC call is unreliable! 
dist <- dist %>% filter(Distance != 0)

# Summarize
dists <- dist %>% 
  group_by(Experiment,DMP) %>%
  sum_stats(Distance)

# Experiment DMP    mean     sd     se median   iqr conf_low conf_high
# <chr>      <chr> <dbl>  <dbl>  <dbl>  <dbl> <dbl>    <dbl>     <dbl>
#   1 ComGar     Age   1361. 45876.  540.     168  397     303.      2419.
# 2 ComGar     Sex    303.   640.   11.8    129  306.    280.       327.
# 3 ComGar     Taxon 1796. 53272. 1478.     152  326.  -1103.      4694.
# 4 HybZon     Taxon  546.  1969.  293.      99  198     -45.7     1137.
# 5 HybZon     Year   544.   741.   76.0    221  671     393.       695.

# pairwise wilcoxon test
pw <- dist %>%
  group_by(Experiment) %>%
  pairwise_wilcox_test(Distance ~ DMP, p.adjust.method = "bonferroni")
ypos <- dists %>% mutate(max_y = conf_high+1e3) %>% select(Experiment,DMP,max_y)
pw <- pw %>%
  left_join(ypos, by = c("Experiment", "group1" = "DMP")) %>%
  rename(y1 = max_y) %>%
  left_join(ypos, by = c("Experiment", "group2" = "DMP")) %>%
  rename(y2 = max_y) %>%
  mutate(y.position = pmax(y1, y2) + 0.05 * max(pmax(y1, y2), na.rm = TRUE),
         label = ifelse(p.adj < 0.001, "***",
                        ifelse(p.adj < 0.01, "**",
                               ifelse(p.adj < 0.05, "*", "ns"))))

dp <- dists %>% 
  ggplot(aes(x = DMP, y = mean, ymin = pmax(conf_low, 0), ymax = conf_high)) +
  geom_point() +
  geom_errorbar() +
  facet_grid(. ~ Experiment, scales = 'free', space = 'free') +
  ylab('Distance from DMP to SNP') +
  theme_bw(base_size=8) +
  stat_pvalue_manual(pw, 
                     label = "label", 
                     y.position = "y.position", 
                     xmin = "group1", 
                     xmax = "group2", 
                     tip.length = 0.01,inherit.aes=FALSE)

ggsave('~/symlinks/hzone/Methylation_Analyses/figures/20250411_SNP_Distance_HZ-CG.pdf',dp,height=2.5,width=4)

```





# Decision Trees 

Take the classified data, and prepare it for RF:

Great tutorial: https://www.kirenz.com/post/2021-02-17-r-classification-tidymodels/

## Prep

```R
.libPaths('~/mambaforge/envs/tidymodels/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_202504/RFS_202504/')
library(dplyr)
library(ggplot2)
library(ggpubr)
library(viridis)
library(GGally)
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(matrixStats)
library(corrr)
options(scipen=999)

all <- read.table('../Full-HZ-CG.DSS.BP-10x.Classified-Annotated-SNPless_20250410.txt',header=T)
all <- all %>% arrange(chr,start)
all <- all %>% mutate_at(grep("HapD|FuLi|Tajima|Dxy|FST|GC",colnames(.)),funs(as.numeric))
all <- subset(all,Region != 'Missing' & Region != '.')
c = read.table('../CG_DSS.10x.stat',header=TRUE) %>% select(site, mstat.cg = 'stat_SpeciesS')
h = read.table('../HZ_DSS.10x.stat',header=TRUE) %>% select(site, mstat.hz = 'stat_Chr18_Hybrid_Index')

#merge data frames
f = left_join(all,c) %>% 
  left_join(.,h) 
s <- f %>% mutate(end = start) %>% select(chr,start,end,mstat.cg,mstat.hz,Region) 
write.table(s,file='../5mC_10x_STAT-Region_202504.txt',quote=F,sep='\t',row.names=F)

#cycle through the experiments 
xps = c('hz','cg')
rfdat = NULL

#add chr length information
genome <- read.table('../genome5.7.bed',header=F) %>% filter(!grepl('scaf|MT',V1)) %>% arrange(desc(V2))
names(genome) <- c('chr','Length')

#count missing by experiment
da = NULL
for (xp in xps){
  d1 = f %>% dplyr::rename(meth = as.symbol(paste0('mstat.',xp))) %>% select(meth,TajimaD,GC,chr,GenStart) %>% 
    drop_na(meth) %>% group_by(chr,GenStart) %>% 
    summarize(GC = mean(GC,na.rm=TRUE),
              val = mean(meth),
              taj = mean(TajimaD))
  gcall = d1 %>% ungroup %>% drop_na(taj) %>% summarize(gca = mean(GC,na.rm=TRUE))
  d1 =  d1 %>% ungroup %>% mutate(Total = n()) %>% filter(is.na(taj)) %>% summarize(nas = n(),
                                                                                    total = Total,
                                                                                    gca = mean(GC,na.rm=TRUE)) %>% unique %>%
    mutate(pct = nas/total * 100,xp = xp,gcall = gcall$gca)
  da = rbind(da,d1)
}
da 
# # A tibble: 2 × 6
# nas total   gca   pct xp    gcall
# <int> <int> <dbl> <dbl> <chr> <dbl>
#   1   228 37883  67.9 0.602 hz     60.0
# 2   356 34441  68.4 1.03  cg     61.4


rfdat <- NULL
for (xp in xps){
  
  title <- case_when(
    xp == 'cg' ~ 'ComGar.RRBS',
    xp == 'hz' ~ 'HybZon.RRBS',
    TRUE ~ NA_character_
  )
  
  cat('Creating output for: ',xp,', ',title,'\n')
  #add chr length and position
  kb = f %>% select(chr,Region,GenStart,GenEnd,HapD,FuLiD,TajimaD,Dxy,FST,GC,as.symbol(paste0('mstat.',xp))) %>% 
    dplyr::rename(val = as.symbol(paste0('mstat.',xp))) %>% 
    mutate(val = abs(val)) %>%   #ensure it's simply the absolute value since we don't care about direction 
    left_join(.,genome) %>%  #add length 
    group_by(chr,Length,Region,GenStart,GenEnd,HapD,FuLiD,TajimaD,Dxy,FST) %>% #calculate mean values within these windows 
    summarize(GC = mean(GC,na.rm=TRUE),
              Statm = mean(val,na.rm=TRUE),
              Statx = max(val,na.rm=TRUE)) %>% 
    arrange(chr,GenStart)
  kb %>% nrow
  kb %>% drop_na(Statm) %>% nrow
  kb %>% drop_na(Statm,TajimaD) %>% nrow
  #on the assumption of LD, impute the closest genetic data for regions which are missing
  kb$Dxy = zoo::na.fill(kb$Dxy, "extend")
  kb$FST = zoo::na.fill(kb$FST, "extend")
  kb$HapD = zoo::na.fill(kb$HapD, "extend")
  kb$TajimaD = zoo::na.fill(kb$TajimaD, "extend")
  kb$FuLiD = zoo::na.fill(kb$FuLiD, "extend")
  kb1 = kb %>% na.omit
  
  #set classifications
  kb1 = kb1 %>% mutate(chr = gsub('chr','',chr),
                       chr = gsub('4A','4.1',chr),
                       chr = gsub('1A','1.1',chr)) %>% 
    filter(chr != 'Z')
  order = kb1 %>% ungroup %>% select(chr,Length) %>% unique %>% arrange(desc(Length))
  kb1$chr = factor(kb1$chr,levels=order$chr)
  cgm <- kb1 %>% ungroup %>% pivot_longer(!c(chr,Length,Region,GenStart,GenEnd,HapD,FuLiD,TajimaD,Dxy,FST,GC),values_to = 'raw',names_to = 'variable')
  cgm = cgm %>% mutate(log = log(ifelse(raw == 0,0.0001,raw)),
                       sqrt = sqrt(raw),
                       cube = raw^(1/3))
  cgm1 = cgm %>% pivot_longer(!c(chr,Length,Region,GenStart,GenEnd,HapD,FuLiD,TajimaD,Dxy,FST,variable,GC),names_to='transformation') %>% 
    mutate(variable = gsub('Statm','Mean',variable),
           variable = gsub('Statx','Max',variable))
  cgm1$transformation = factor(cgm1$transformation,levels=c('raw','log','sqrt','cube')) 
  rawp = cgm1 %>% ggplot(aes(y=value,fill=transformation))+
    geom_density()+scale_fill_manual(values=viridis(4,option='inferno'))+
    facet_wrap(variable~transformation,scales='free',ncol = 4)+
    theme_classic()
  
  #apply transformations to SD and variance
  xfm = cgm %>% mutate(valueuse = sqrt(raw)) %>% 
    mutate(variable = gsub('Statm','Mean',variable),
           variable = gsub('Statx','Max',variable))
  #replot
  xfmp  = xfm %>% select(variable,valueuse,raw) %>% 
    pivot_longer(!c(variable)) %>% 
    mutate(name = case_when(
      name == 'valueuse' ~ 'log(Raw)',
      name == 'raw' ~ 'Raw')) %>% 
    ggplot(aes(x=value,fill=name))+
    geom_density()+scale_fill_manual(values=viridis(4,option='inferno'))+
    facet_wrap(name~variable,scales='free',ncol=4)+
    ggtitle(paste0('Continuous: ',title))+
    theme_classic()
  assign(paste0(xp,'_values'),xfmp)
  
  #set classifications for each variable
  ##### BINARY
  dat = xfm %>% ungroup %>% select(!c(log,sqrt,cube))
  bdat = dat
  thresh = bdat %>% group_by(variable) %>% summarize(Threshold = quantile(valueuse,0.8)) %>% unique
  thresh
  bdat = left_join(bdat,thresh)
  bdat = bdat %>% mutate(Classified = ifelse(valueuse >=Threshold, 'High','Low'))
  binp = bdat %>% ggplot(aes(x=valueuse,fill=variable)) + 
    facet_wrap(~variable,scales='free')+
    geom_histogram()+
    scale_fill_viridis(discrete=TRUE,option='mako')+
    geom_vline(data=thresh, aes(xintercept=Threshold),lty=2,lwd=0.5,col='yellow')+
    ylab('')+theme_classic()+
    ggtitle(paste0('Binary: ',title))
  assign(paste0(xp,'_thresh'),binp)
  
  #merge them 
  bd = bdat %>% select(-contains('Thresh')) %>% dplyr::rename(Continuous = valueuse)
  ad = bd %>% mutate(Relativepos = GenStart / Length)
  ad$Subset = toupper(xp)
  rfdat = rbind(rfdat,ad)
}

pdf('~/symlinks/hzone/Methylation_Analyses/figures/20250411_Raw_Distributions.pdf',height=5,width=5)
ggarrange(cg_values,hz_values,common.legend = TRUE,nrow=3)
dev.off()

pdf('~/symlinks/hzone/Methylation_Analyses/figures/20250411_Raw_Thresholds.pdf',height=6,width=4)
ggarrange(cg_thresh,hz_thresh,common.legend = TRUE,nrow=3)
dev.off()

#save frame 
write.table(rfdat,file='RF-Input-250411.txt',quote=F,sep='\t',row.names=F)

#output correlations among variables
cord = NULL
xps = c('HZ','CG')
for (xp in xps){
  d1 = rfdat %>% filter(Subset == xp & variable == 'Mean')
  np = d1 %>% 
    select(Length,Relativepos,Region,HapD,FuLiD,TajimaD,Dxy,FST,GC) %>% 
    correlate(method='spearman') %>%    # Create correlation data frame (cor_df)
    corrr::network_plot(min_cor = 0,legend='range')
  assign(paste0(xp,'_cor'),np)
  
  x = d1 %>%
    select(Length,Relativepos,Region,HapD,FuLiD,TajimaD,Dxy,FST,GC) %>%
    corrr::correlate() %>%    # Create correlation data frame (cor_df)
    rearrange()  # rearrange by correlations
  x$experiment = xp
  cord = rbind(cord,x)
}

pdf('~/symlinks/hzone/Methylation_Analyses/figures/20250411_Correlations_RF_Input.pdf',height=6,width=12)
ggarrange(CG_cor,HZ_cor,labels=c('ComGar.RRBS','HybZon.RRBS'),nrow=1)
dev.off()

write.table(cord,file='~/symlinks/hzone/Methylation_Analyses/figures/20250411_Correlations_RF_Input.txt',quote=F,sep='\t',row.names=F)

```

## Run models

```R
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
.libPaths('~/mambaforge/envs/tidymodels/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_202504/RFS_202504/')
#### Randomforest Classification
library(tidymodels)
library(GGally)
library(gt)
library(visdat)
library(themis)
library(skimr)
library(viridis)
library(finetune)
library(xgboost)
require(doParallel)
library(kknn)
library(forcats)
library(vip)

#Parallelisation
cores <- 5
cl <- makePSOCKcluster(cores)
registerDoParallel(cl)

dats <- read.table('RF-Input-250411.txt',header=T)
type = args[1]  #e.g. Binary,Regression
vari = args[2] #e.g. Mean,Max
xp = args[3] #subset experiment, CG HZ
iter = args[4] #iteration

dat = dats %>% filter(variable == vari & Subset == xp)
if (type == 'Regression') {
  cat('Response type is Continuous, taking continuous values\n')
  rf = dat %>% select(c(Continuous,Length,Relativepos,Region,HapD,FuLiD,TajimaD,Dxy,FST,GC))
  rf = rf %>% dplyr::rename(Response = Continuous)
} else {
  cat('Response type is Binary, taking classified values\n')
  rf = dat %>% select(c(as.symbol(type),Length,Relativepos,Region,HapD,FuLiD,TajimaD,Dxy,FST,GC))
  rf = rf %>% dplyr::rename(Response = Classified)
}

#for exploration, run on subset
#rf = rf %>% sample_n(5000)
# convert to numeric
rf = rf %>% mutate_at(grep("Length|Relative|HapD|Fu|Taj|Dxy|FST|GC",colnames(.)),funs(as.numeric))

# convert all remaining character variables to factors
rf = rf %>% mutate(across(where(is.character), as.factor))
skim(rf)

#split into training and testing datasets
data_split <- initial_split(rf, strata = "Response", prop = 0.75)

train_explore <- training(data_split)
test_explore  <- testing(data_split)
data_explore = train_explore

# Generate resamples and repeat
cv_folds = vfold_cv(train_explore, v = 5, strata = Response)

# define the recipe
mc_recipe =  # which consists of the formula (outcome ~ predictors)
  recipe(Response ~ ., data = rf) %>%
  # and some pre-processing steps
  step_novel(all_nominal(), -all_outcomes()) %>% # converts all nominal variables to factors and takes care of other issues related to categorical variables
  step_normalize(Length, Relativepos, HapD, FuLiD, TajimaD, Dxy, FST, GC, -all_outcomes()) %>% # normalizes (center and scales) the numeric variables to have a standard deviation of one and a mean of zero. (i.e., z-standardization)
  step_dummy(all_nominal(), -all_outcomes()) %>% # converts our factor column ocean_proximity into numeric binary (0 and 1) variables
  step_zv(all_numeric(), -all_outcomes()) %>% # removes any numeric variables that have zero variance
  step_corr(all_predictors(), threshold = 0.8, method = "spearman") # will remove predictor variables that have large correlations with other predictor variables

#check the pre-processing
examine_preprocessed <- mc_recipe %>% # apply the recipe to the training data
  prep(rf) %>%  # extract the pre-processed training dataset
  juice()
examine_preprocessed

#set the models
if (type == 'Regression') {
  cat('Response type is Regression, using regression engines\n')
  #boosted trees
  xgb_mod = boost_tree(trees = tune(),min_n = tune(),mtry = tune(),learn_rate = tune()) %>% set_engine("xgboost",importance="permutation") %>% set_mode("regression")
  xgb_wf <- workflow(mc_recipe, xgb_mod)
  xgb_rs <- tune_race_anova(xgb_wf,resamples = cv_folds,grid = 15, metrics = metric_set(rsq),control = control_race(verbose_elim = TRUE))
  #plot_race(xgb_rs)
  show_best(xgb_rs)
  xgb_last = xgb_wf %>% finalize_workflow(select_best(xgb_rs, "rsq")) %>% last_fit(data_split)
  xgb_pred = collect_predictions(xgb_last) %>% select(.pred,Response)
  xgb_pred$Engine = 'XGBoost'
  xgb_vi = extract_workflow(xgb_last) %>% extract_fit_parsnip() %>% vi()
  xgb_vi$Engine = 'XGBoost'
  xgb_fit = xgb_last %>% collect_metrics()
  xgb_fit$Engine = 'XGBoost'
  
  #randomforest
  rf_mod = rand_forest(trees = tune(),min_n = tune(),mtry = tune()) %>% set_engine("ranger",importance="permutation") %>% set_mode("regression")
  rf_wf = workflow(mc_recipe, rf_mod)
  rf_rs = tune_race_anova(rf_wf,resamples = cv_folds,grid = 15, metrics = metric_set(rsq),control = control_race(verbose_elim = TRUE))
  show_best(rf_rs)
  rf_last = rf_wf %>% finalize_workflow(select_best(rf_rs, "rsq")) %>% last_fit(data_split)
  rf_pred = collect_predictions(rf_last) %>% select(.pred,Response)
  rf_pred$Engine = 'Random forest'
  rf_vi = extract_workflow(rf_last) %>% extract_fit_parsnip() %>% vi()
  rf_vi$Engine = 'Random forest'
  rf_fit = rf_last %>% collect_metrics()
  rf_fit$Engine = 'Random forest'
  
} else {
  cat('Response type is Binary, using classification engines\n')
  rf_fit = NULL
  rf_vi = NULL
  rf_pred = NULL
  #boosted trees
  xgb_mod = boost_tree(trees = tune(),min_n = tune(),mtry = tune(),learn_rate = tune()) %>% set_engine("xgboost",importance="permutation") %>% set_mode("classification")
  xgb_wf <- workflow(mc_recipe, xgb_mod)
  xgb_rs <- tune_race_anova(xgb_wf,resamples = cv_folds,grid = 15, metrics = metric_set(roc_auc),control = control_race(verbose_elim = TRUE))
  #plot_race(xgb_rs)
  show_best(xgb_rs)
  xgb_last = xgb_wf %>% finalize_workflow(select_best(xgb_rs, "roc_auc")) %>% last_fit(data_split)
  xgb_pred = collect_predictions(xgb_last) %>% select(.pred_class,Response)
  xgb_pred$Engine = 'XGBoost'
  xgb_vi = extract_workflow(xgb_last) %>% extract_fit_parsnip() %>% vi()
  xgb_vi$Engine = 'XGBoost'
  xgb_fit = xgb_last %>% collect_metrics()
  xgb_fit$Engine = 'XGBoost'
  
}

#save fit
fit = rbind(xgb_fit,rf_fit)
names(fit) = c('Metric','Estimator','Estimate','Config','Engine')
fitdf = data.frame(fit,Type = type,Response = vari,Subset = xp,Subset = xp, Iteration=iter)
fitdf
write.table(fitdf,file=paste0('results/',type,'_',vari,'_',xp,'_',iter,'_FIT.txt'),quote=F,sep='\t',row.names=F)

#save vi
vi = rbind(xgb_vi,rf_vi)
vidf = data.frame(vi,Type = type,Response = vari,Subset = xp,Iteration=iter)
vidf
write.table(vidf,file=paste0('results/',type,'_',vari,'_',xp,'_',iter,'_VI.txt'),quote=F,sep='\t',row.names=F)

#save predictions
preds = rbind(xgb_pred,rf_pred)
names(preds) = c('Predicted','Truth','Engine')
preddf = data.frame(preds,Type = type,Response = vari,Subset = xp,Iteration=iter)
#preddf
write.table(preddf,file=paste0('results/',type,'_',vari,'_',xp,'_',iter,'_PRED.txt'),quote=F,sep='\t',row.names=F)

```

Prepare the sbatch:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --time=48:00:00

T=$1
V=$2
S=$3
I=$4

echo "Running Randomforest.R with: Type=$T, Variable=$V, Subset=$S, Iteration=$I"

Rscript Randomforest.R $T $V $S $I

```



Submit the many loops with:

```bash
for a in $(cat Types.list); do  for b in $(cat Variables.list); do for c in $(cat Subsets.list); do for d in $(cat Iterations.list); do  sbatch -J ${a}_${b}_${c}_${d} randomforest.sh $a $b $c $d;  done;  done; done; done
```

Merge results:

```bash
mergem *VI*txt > Master.vi
mergem *FIT*txt > Master.fit
mergem Class*PRED* > Classification_Predictions.txt
mergem Regre*PRED* > Regression_Predictions.txt
```

## Plot Results

```R
.libPaths('~/mambaforge/envs/r/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_202504/RFS_202504/results')
setwd('G:/My Drive/Research/Crow_Hybrid_Epigenetics/rf/2025_04/')
#### Plot decision tree results
library(tidyverse)
library(viridis)
library(forcats)
library(ggpubr)
library(rstatix)
library(RColorBrewer)
library(gmodels)

#evaluate fit
fit = read.table('Master.fit',header=T,sep='\t') %>% select(-c(Subset.1))
fs = fit %>% group_by(Subset,Response,Type,Metric) %>% 
  summarize(mean = ci(Estimate)[1],
            hi = ci(Estimate)[2],
            lo = ci(Estimate)[3],
            sd = ci(Estimate)[4])

# Subset Response Type       Metric     mean     hi     lo        sd
# <chr>  <chr>    <chr>      <chr>     <dbl>  <dbl>  <dbl>     <dbl>
#   1 CG     Max      Classified accuracy 0.800  0.798  0.803  0.000507 
# 2 CG     Max      Classified roc_auc  0.682  0.674  0.690  0.00178  
# 3 CG     Max      Regression rmse     0.370  0.369  0.372  0.000620 
# 4 CG     Max      Regression rsq      0.102  0.101  0.104  0.000720 
# 5 CG     Mean     Classified accuracy 0.800  0.800  0.801  0.000127 
# 6 CG     Mean     Classified roc_auc  0.704  0.700  0.709  0.00105  
# 7 CG     Mean     Regression rmse     0.243  0.241  0.245  0.000831 
# 8 CG     Mean     Regression rsq      0.194  0.189  0.198  0.00171  
# 9 HZ     Max      Classified accuracy 0.800  0.800  0.800  0.0000969
# 10 HZ     Max      Classified roc_auc  0.644  0.631  0.657  0.00299  
# 11 HZ     Max      Regression rmse     0.292  0.291  0.293  0.000404 
# 12 HZ     Max      Regression rsq      0.0964 0.0920 0.101  0.00171  
# 13 HZ     Mean     Classified accuracy 0.800  0.800  0.800  0        
# 14 HZ     Mean     Classified roc_auc  0.663  0.645  0.682  0.00434  
# 15 HZ     Mean     Regression rmse     0.167  0.164  0.171  0.00129  
# 16 HZ     Mean     Regression rsq      0.0794 0.0791 0.0798 0.000132 

fs = fs %>% 
  mutate(Subset = gsub('CG','ComGar',Subset),
         Subset = gsub('HZ','HybZon',Subset))
fs$Subset <- factor(fs$Subset,levels=c('ComGar','HybZon'))
fs %>% filter(Response == 'Mean' & Type == 'Regression' & Metric == 'rsq')
fpr = fs %>% 
  filter(Type == 'Regression') %>% 
  ggplot(aes(x=Response,y=mean,col=Subset,shape=Subset))+
  geom_point(size=3,position=position_dodge(width=0.5))+
  geom_errorbar(aes(ymin=lo,ymax=hi),width=0.5,position=position_dodge(width=0.5))+
  ylab('Evaluation Metric')+xlab('Response (DNA Methylation Divergence Within 5KB Window)')+
  scale_color_manual(values=brewer.pal(3,'Set2'))+
  ggtitle('Regression')+
  theme_bw()+
  theme(legend.position='none')+
  facet_grid(Metric~Type,scales='free')
fpr
fpc = fs %>% 
  filter(Type != 'Regression') %>% 
  ggplot(aes(x=Response,y=mean,col=Subset,shape=Subset))+
  geom_point(size=3,position=position_dodge(width=0.5))+
  geom_errorbar(aes(ymin=lo,ymax=hi),width=0.5,position=position_dodge(width=0.5))+
  ylab('Evaluation Metric')+xlab('Response (DNA Methylation Divergence Within 5KB Window)')+
  scale_color_manual(values=brewer.pal(3,'Set2'))+
  ggtitle('Classification')+
  theme_bw()+
  #theme(legend.position='none')+
  facet_grid(Metric~Type,scales='free')
fpc
fpa = ggarrange(fpr,fpc,widths = c(0.75,1.2))
fpa

#save
pdf('~/symlinks/hzone/Methylation_Analyses/figures/20250414_RF_Fit.pdf',height=5,width=8)
fpa
dev.off()

#evalute VI
vi = read.table('Master.vi',header=T,sep='\t')
#this won't change anything for XGBoost, but to plot together we will scale the RF covariate importance to sum to 1 (as with XGBoost)
vi = vi %>% group_by(Subset,Response,Engine,Type,Iteration) %>% mutate(Importance = ifelse(Importance < 0,0,Importance),
                                                                       TotalImp = sum(Importance),
                                                                       PctImp = Importance / TotalImp,
                                                                       Variable = gsub('GC','GC-Content',Variable),
                                                                       Variable = gsub('Length','Chromosome Length',Variable),
                                                                       Variable = gsub('Relativepos','Chromosomal Position',Variable),
                                                                       Variable = gsub('HapD','Haplotype Diversity',Variable),
                                                                       Variable = gsub('TajimaD',"Tajima's D",Variable),
                                                                       Variable = gsub("FuLiD","Fu & Li's D*",Variable))
vis = vi %>% group_by(Subset,Response,Type,Variable) %>% 
  summarize(mean = ci(PctImp)[1],
            hi = ci(PctImp)[2],
            lo = ci(PctImp)[3],
            sd = ci(PctImp)[4])
vis = vis %>% 
  mutate(Subset = gsub('CG','ComGar',Subset),
         Subset = gsub('HZ','HybZon',Subset))
vis$Subset <- factor(vis$Subset,levels=c('ComGar','HybZon'))
#find most important overall, we will sort by this
imps = vi %>% group_by(Response,Type,Variable) %>% summarize(mean=mean(PctImp)) %>% arrange(desc(mean)) %>% mutate(ORDER = row_number())
vis = left_join(vis,imps %>% select(Response,Type,Variable,ORDER))
vip = vis %>% 
  #filter(Response == 'Mean' & Type == 'Regression') %>% 
  mutate(NewOrd = fct_reorder(Variable,ORDER,.desc = TRUE)) %>% 
  ggplot(aes(y=NewOrd,x=mean,fill=Subset))+
  geom_bar(stat='identity',position=position_dodge(width=0.9))+
  geom_errorbar(aes(xmin=ifelse(mean-sd<0,0,mean-sd),xmax=mean+sd),width=0.5,position=position_dodge(width=0.9))+
  scale_fill_manual(values=c('#b8de29','#33638d'))+
  theme_bw()+
  #theme(legend.position='none')+
  ylab('')+xlab('Permutation Importance')+
  facet_grid(Type~Response,scales='free')+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))
vip

#save
pdf('20250414_RF_VIP_Regression-All.pdf',height=6,width=6)
vip 
dev.off()


viS = vi %>% group_by(Subset,Response,Type,Variable,Engine) %>% 
  #filter(Response == 'Mean' & Type == 'Regression') %>% 
  summarize(mean = ci(PctImp)[1],
            lo = ci(PctImp)[2],
            hi = ci(PctImp)[3],
            se = ci(PctImp)[4])
write.table(viS,file='~/symlinks/hzone/Methylation_Analyses/figures/20250414_Permutation_Importance.txt',quote=F,sep='\t',row.names=F)

# Across all experiments, as reported in main text
vi %>% group_by(Variable) %>% 
  filter(Response == 'Mean' & Type == 'Regression') %>% 
  summarize(mean = ci(PctImp)[1],
            lo = ci(PctImp)[2],
            hi = ci(PctImp)[3],
            se = ci(PctImp)[4])

# Variable                 mean       lo      hi       se
# <chr>                   <dbl>    <dbl>   <dbl>    <dbl>
#   1 Chromosomal Position 0.0556   0.0492   0.0620  0.00289 
# 2 Chromosome Length    0.0728   0.0653   0.0803  0.00341 
# 3 Dxy                  0.0959   0.0638   0.128   0.0146  
# 4 FST                  0.0172   0.0110   0.0234  0.00283 
# 5 Fu & Li's D*         0.0480   0.0390   0.0570  0.00409 
#  6 GC-Content           0.251    0.214    0.287   0.0166  
#  7 Haplotype Diversity  0.0895   0.0669   0.112   0.0103  
#  8 Region_Intergenic    0.0276   0.0200   0.0352  0.00345 
#  9 Region_Intron        0.0179   0.0115   0.0243  0.00290 
# 10 Region_Promoter      0.267    0.219    0.316   0.0219  
# 11 Region_Repeat        0.000589 0.000146 0.00103 0.000201
# 12 Tajima's D           0.0566   0.0429   0.0704  0.00625 

#save
pdf('~/symlinks/hzone/Methylation_Analyses/figures/20250414_RF_VIP_Regression.pdf',height=6,width=6)
vip 
dev.off()


#Also plot for ALL subsets, experiments, types 
ALL = vi %>% group_by(Subset,Response,Type,Variable) %>% 
  summarize(mean = ci(PctImp)[1],
            hi = ci(PctImp)[2],
            lo = ci(PctImp)[3],
            sd = ci(PctImp)[4])
ALL = ALL %>% 
  mutate(Subset = gsub('CG','ComGar',Subset),
         Subset = gsub('HZ','HybZon',Subset))
ALL$Subset <- factor(ALL$Subset,levels=c('ComGar','HybZon'))
#find most important overall, we will sort by this
ALL = left_join(ALL,imps %>% ungroup %>% select(Variable,ORDER))
ALL$Type = factor(ALL$Type,levels=c('Regression','Classified'))
ALLp = ALL %>% 
  mutate(NewOrd = fct_reorder(Variable,ORDER,.desc = TRUE)) %>% 
  ggplot(aes(y=NewOrd,x=mean,fill=Subset))+
  geom_bar(stat='identity',position=position_dodge(width=0.9))+
  geom_errorbar(aes(xmin=ifelse(mean-sd<0,0,mean-sd),xmax=mean+sd),width=0.5,position=position_dodge(width=0.9))+
  scale_fill_manual(values=brewer.pal(3,'Set2'))+
  theme_bw()+
  #theme(legend.position='none')+
  ylab('')+xlab('Permutation Importance')+
  facet_grid(Response~Type,scales='free')+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))
ALLp

#save
pdf('~/symlinks/hzone/Methylation_Analyses/figures/20250414_RF_VIP_All.pdf',height=6,width=10)
ALLp
dev.off()

#evalute predictions
conf = read.table('Regression_Predictions.txt',header=T,sep='\t')
cp = conf %>% 
  ggplot(aes(x=Predicted,y=Truth,col=Subset))+
  #geom_point(alpha = 0.5)+
  geom_smooth()+
  #geom_density_2d()+
  scale_color_manual(values=c(viridis(3)))+
  theme_bw()+xlab('Predicted')+ylab('Truth')+
  theme(legend.position = "none")+
  facet_wrap(Subset~Response,scales='free',nrow=3)
cp

#evalute predictions, classification
confclass = read_tsv('Classification_Predictions.txt')
confclass %>% group_by(Predicted,Truth,Engine,Response,Subset) %>% 
  summarize(count = n()) %>% 
  ggplot(aes(x=Predicted,y=Truth,fill=count))+
  geom_tile()+
  scale_fill_continuous(low='yellow',high='red')+
  theme_bw()+xlab('Predicted')+ylab('Truth')+
  facet_grid(Response+Engine~Subset,scales='free')


#save
png('~/symlinks/hzone/Methylation_Analyses/figures/20250414_RF_Predictions.png',height=5,width=3.5,res=600,units='in')
cp
dev.off()

```



## Rho: 5mC PopGen Correlations

```bash
.libPaths('~/mambaforge/envs/r/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_202504/RFS_202504')
#### PopGen ~ 5mC spearman correlations
library(tidyverse)
library(viridis)
library(forcats)
library(ggpubr)
library(rstatix)
library(gmodels)
library(RColorBrewer)

#plot some of the interesting comparisons
rf = read_tsv('RF-Input-250411.txt')

### Correlations: 5mC / genetic variation
cat('Working on correlations for genetic variation, 5mc \n')
rf1 = rf %>% filter(variable == 'Mean') %>% select(chr,GenStart,GenEnd,raw,HapD,FuLiD,TajimaD,Dxy,FST,Subset)
rf1 = rf1 %>% mutate(Focal_Region = ifelse(chr == '18' & GenStart > 8070001 & GenStart < 10070000,'Focal','Undiverged'))
rf2 = rf1 %>%
  pivot_longer(!c(chr,GenStart,GenEnd,raw,Subset,Focal_Region)) %>%
  na.omit()
xps = c('CG','HZ','WGBS')

fullboots = NULL
options(dplyr.summarise.inform = FALSE)
set.seed(123)
B <- 999 # number of bootstrap replicates

for (xp in xps) {
  for (vr in unique(rf2$name)) {
    cat('Working on Experiment & Variable: ',xp,' ',vr,'\n')
    bd = rf2 %>% filter(Subset == xp & name == vr)
    n = bd %>% filter(Focal_Region == 'Focal') %>% nrow
    bdat=NULL
    ## Run the bootstrap procedure:
    for (b in 1:B) {
      cat('Iteration: ',b,'\n')
      d1 = bd %>% group_by(Focal_Region) %>% slice_sample(n=n,replace=TRUE) %>% summarize(cor = cor(raw,value,method='spearman'))
      bdat = rbind(bdat,d1 )
    }
    bdat$Experiment = xp
    bdat$Variable = vr
    names(bdat) = c('Focal_Region','Boots','Experiment','Variable')
    fullboots = rbind(fullboots,bdat)
  }
}

write.table(fullboots,file='Genetic5mC_Bootstraped_Results-2020250411.txt',quote=F,sep='\t',row.names=F)
fullboots = read.table('Genetic5mC_Bootstraped_Results-2020250411.txt',header=TRUE)

#add counts for windows 
features = rf2 %>% group_by(Subset,Focal_Region,name) %>% count()
names(features)= c('Experiment','Focal_Region','Variable')

cdm = fullboots %>% 
  mutate(Experiment = gsub('CG','ComGar',Experiment),
         Experiment = gsub('HZ','HybZon',Experiment))
cdm$Experiment <- factor(cdm$Experiment,levels=c('HybZon','ComGar'))
#calculate CIs
cis = cdm %>% 
  group_by(Experiment,Variable,Focal_Region) %>% 
  summarize(LowerBound = quantile(Boots,0.025), #can change based on your significance threshold 
            UpperBound = quantile(Boots,0.975))
cis = cis %>% group_by(Experiment,Variable) %>% 
  mutate(Signif = ifelse(LowerBound > 0 | UpperBound < 0,'*','n.s.'))

genp = cdm %>% 
  ggplot(aes(x=Variable,y=Boots,fill=Focal_Region))+
  geom_boxplot(alpha=0.75)+
  geom_text(data=cis,aes(x=Variable,col=Focal_Region,group=Focal_Region,y=Inf,label=Signif,size=Signif),col='black',vjust=1.5,position=position_dodge(width=1))+
  scale_size_manual(values=c(6,4))+
  scale_fill_manual(values=c('darkorchid2','grey60'))+
  facet_grid(Experiment~.,scales='free')+
  geom_hline(yintercept=0,lty=2)+
  ylab("Spearman's Rho: Taxon Methylation Divergence") + xlab('')+
  scale_y_continuous(expand = expansion(mult = .25)) + #expand y axis slightly 
  theme_bw()
genp

pdf('~/symlinks/hzone/Methylation_Analyses/figures/20250414_GenCorBootBootstraps.pdf',height=4,width=5)
genp
dev.off()

write.table(cis,file='../../Correlations_95CIs.txt',quote=F,sep='\t',row.names=F)
```



## Including ComGar.WGBS

### Prep 

```R
.libPaths('~/mambaforge/envs/tidymodels/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_Sept/RF_JAN25')
library(dplyr)
library(ggplot2)
library(ggpubr)
library(viridis)
library(GGally)
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(matrixStats)
options(scipen=999)

all <- read.table('../Full.DSS.BP-10x.Classified-Annotated_20250107.txt',header=T)
all <- all %>% arrange(chr,start)
all <- all %>% mutate_at(grep("HapD|FuLi|Tajima|Dxy|FST|GC",colnames(.)),funs(as.numeric))
all <- subset(all,Region != 'Missing' & Region != '.')
c = read.table('../CG.DSS-10x.stat',header=TRUE) %>% select(site, mstat.cg = 'CG_stat')
w = read.table('../WGBS.DSS-10x.stat',header=TRUE) %>% select(site, mstat.wgbs = 'stat_SpeciesS')
h = read.table('../HZ.DSS.BP-Parentals-YearChr18HI_10x.stat',header=TRUE) %>% select(site, mstat.hz = 'stat_Chr18_Hybrid_Index')

#merge data frames
f = left_join(all,c) %>% 
  left_join(.,w) %>% 
  left_join(.,h) 
s <- f %>% mutate(end = start) %>% select(chr,start,end,mstat.wgbs,mstat.cg,mstat.hz,Region) 
write.table(s,file='../5mC_10x_STAT-Region_20250107.txt',quote=F,sep='\t',row.names=F)

#examine some top hits from the methylation test statistic 
# f2 %>% select(chr,start,contains(c('mstat.cg'))) %>% slice_max(mstat.cg,n = 10)
# f2 %>% filter(chr == 'chr18' & start == 770644) %>% select(chr,contains(c('S_','D_'),ignore.case = F)) %>% 
#   select(chr,contains(c('.cg'))) %>% 
#   pivot_longer(!c(chr)) %>% 
#   separate(name, into=c('species','locality','ID','tissue','age','sex','xp')) %>% 
#   ggplot(aes(x=species,y=value,col=age,shape=sex))+
#   geom_point(position=position_jitter())+
#   theme_bw()

#cycle through the experiments 
xps = c('hz','wgbs','cg')
rfdat = NULL

#add chr length information
genome <- read.table('../genome5.7.bed',header=F) %>% filter(!grepl('scaf|MT',V1)) %>% arrange(desc(V2))
names(genome) <- c('chr','Length')

#count missing by experiment
da = NULL
for (xp in xps){
  d1 = f %>% dplyr::rename(meth = as.symbol(paste0('mstat.',xp))) %>% select(meth,TajimaD,GC,chr,GenStart) %>% 
    drop_na(meth) %>% group_by(chr,GenStart) %>% 
    summarize(GC = mean(GC,na.rm=TRUE),
              val = mean(meth),
              taj = mean(TajimaD))
  gcall = d1 %>% ungroup %>% drop_na(taj) %>% summarize(gca = mean(GC,na.rm=TRUE))
  d1 =  d1 %>% ungroup %>% mutate(Total = n()) %>% filter(is.na(taj)) %>% summarize(nas = n(),
                                                                                    total = Total,
                                                                                    gca = mean(GC,na.rm=TRUE)) %>% unique %>%
    mutate(pct = nas/total * 100,xp = xp,gcall = gcall$gca)
  da = rbind(da,d1)
}
da 

rfdat <- NULL
for (xp in xps){
  
  title <- case_when(
    xp == 'cg' ~ 'ComGar.RRBS',
    xp == 'wgbs' ~ 'ComGar.WGBS',
    xp == 'hz' ~ 'HybZon.RRBS',
    TRUE ~ NA_character_
  )
  
  cat('Creating output for: ',xp,', ',title,'\n')
  #add chr length and position
  kb = f %>% select(chr,Region,GenStart,GenEnd,HapD,FuLiD,TajimaD,Dxy,FST,GC,as.symbol(paste0('mstat.',xp))) %>% 
    dplyr::rename(val = as.symbol(paste0('mstat.',xp))) %>% 
    mutate(val = abs(val)) %>%   #ensure it's simply the absolute value since we don't care about direction 
    left_join(.,genome) %>%  #add length 
    group_by(chr,Length,Region,GenStart,GenEnd,HapD,FuLiD,TajimaD,Dxy,FST) %>% #calculate mean values within these windows 
    summarize(GC = mean(GC,na.rm=TRUE),
              Statm = mean(val,na.rm=TRUE),
              Statx = max(val,na.rm=TRUE)) %>% 
    arrange(chr,GenStart)
  kb %>% nrow
  kb %>% drop_na(Statm) %>% nrow
  kb %>% drop_na(Statm,TajimaD) %>% nrow
  #on the assumption of LD, impute the closest genetic data for regions which are missing
  kb$Dxy = zoo::na.fill(kb$Dxy, "extend")
  kb$FST = zoo::na.fill(kb$FST, "extend")
  kb$HapD = zoo::na.fill(kb$HapD, "extend")
  kb$TajimaD = zoo::na.fill(kb$TajimaD, "extend")
  kb$FuLiD = zoo::na.fill(kb$FuLiD, "extend")
  kb1 = kb %>% na.omit
  
  #set classifications
  kb1 = kb1 %>% mutate(chr = gsub('chr','',chr),
                       chr = gsub('4A','4.1',chr),
                       chr = gsub('1A','1.1',chr)) %>% 
    filter(chr != 'Z')
  order = kb1 %>% ungroup %>% select(chr,Length) %>% unique %>% arrange(desc(Length))
  kb1$chr = factor(kb1$chr,levels=order$chr)
  cgm <- kb1 %>% ungroup %>% pivot_longer(!c(chr,Length,Region,GenStart,GenEnd,HapD,FuLiD,TajimaD,Dxy,FST,GC),values_to = 'raw',names_to = 'variable')
  cgm = cgm %>% mutate(log = log(ifelse(raw == 0,0.0001,raw)),
                       sqrt = sqrt(raw),
                       cube = raw^(1/3))
  cgm1 = cgm %>% pivot_longer(!c(chr,Length,Region,GenStart,GenEnd,HapD,FuLiD,TajimaD,Dxy,FST,variable,GC),names_to='transformation') %>% 
    mutate(variable = gsub('Statm','Mean',variable),
           variable = gsub('Statx','Max',variable))
  cgm1$transformation = factor(cgm1$transformation,levels=c('raw','log','sqrt','cube')) 
  rawp = cgm1 %>% ggplot(aes(y=value,fill=transformation))+
    geom_density()+scale_fill_manual(values=viridis(4,option='inferno'))+
    facet_wrap(variable~transformation,scales='free',ncol = 4)+
    theme_classic()
  
  #apply transformations to SD and variance
  xfm = cgm %>% mutate(valueuse = sqrt(raw)) %>% 
    mutate(variable = gsub('Statm','Mean',variable),
           variable = gsub('Statx','Max',variable))
  #replot
  xfmp  = xfm %>% select(variable,valueuse,raw) %>% 
    pivot_longer(!c(variable)) %>% 
    mutate(name = case_when(
      name == 'valueuse' ~ 'log(Raw)',
      name == 'raw' ~ 'Raw')) %>% 
    ggplot(aes(x=value,fill=name))+
    geom_density()+scale_fill_manual(values=viridis(4,option='inferno'))+
    facet_wrap(name~variable,scales='free',ncol=4)+
    ggtitle(paste0('Continuous: ',title))+
    theme_classic()
  assign(paste0(xp,'_values'),xfmp)
  
  #set classifications for each variable
  ##### BINARY
  dat = xfm %>% ungroup %>% select(!c(log,sqrt,cube))
  bdat = dat
  thresh = bdat %>% group_by(variable) %>% summarize(Threshold = quantile(valueuse,0.8)) %>% unique
  thresh
  bdat = left_join(bdat,thresh)
  bdat = bdat %>% mutate(Classified = ifelse(valueuse >=Threshold, 'High','Low'))
  binp = bdat %>% ggplot(aes(x=valueuse,fill=variable)) + 
    facet_wrap(~variable,scales='free')+
    geom_histogram()+
    scale_fill_viridis(discrete=TRUE,option='mako')+
    geom_vline(data=thresh, aes(xintercept=Threshold),lty=2,lwd=0.5,col='yellow')+
    ylab('')+theme_classic()+
    ggtitle(paste0('Binary: ',title))
  assign(paste0(xp,'_thresh'),binp)
  
  #merge them 
  bd = bdat %>% select(-contains('Thresh')) %>% dplyr::rename(Continuous = valueuse)
  ad = bd %>% mutate(Relativepos = GenStart / Length)
  ad$Subset = toupper(xp)
  rfdat = rbind(rfdat,ad)
}

pdf('20250107_Raw_Distributions.pdf',height=5,width=5)
ggarrange(wgbs_values,cg_values,hz_values,common.legend = TRUE,nrow=3)
dev.off()

pdf('20250107_Raw_Thresholds.pdf',height=6,width=4)
ggarrange(wgbs_thresh,cg_thresh,hz_thresh,common.legend = TRUE,nrow=3)
dev.off()

#save frame 
write.table(rfdat,file='RF-Input-250107.txt',quote=F,sep='\t',row.names=F)

#output correlations among variables
cord = NULL
xps = c('HZ','CG','WGBS')
for (xp in xps){
  d1 = rfdat %>% filter(Subset == xp & variable == 'Mean')
  np = d1 %>% 
    select(Length,Relativepos,Region,HapD,FuLiD,TajimaD,Dxy,FST,GC) %>% 
    correlate(method='spearman') %>%    # Create correlation data frame (cor_df)
    corrr::network_plot(min_cor = 0,legend='range')
  assign(paste0(xp,'_cor'),np)
  
  x = d1 %>%
    select(Length,Relativepos,Region,HapD,FuLiD,TajimaD,Dxy,FST,GC) %>%
    corrr::correlate() %>%    # Create correlation data frame (cor_df)
    rearrange()  # rearrange by correlations
  x$experiment = xp
  cord = rbind(cord,x)
}

pdf('20250107_Correlations_RF_Input.pdf',height=6,width=12)
ggarrange(WGBS_cor,CG_cor,HZ_cor,labels=c('ComGar.WGBS','ComGar.RRBS','HybZon.RRBS'),nrow=1)
dev.off()

write.table(cord,file='20250107_Correlations_RF_Input.txt',quote=F,sep='\t',row.names=F)

```

### Run RF & XGB

Once the data is prepared, grow the forests: 

```R
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
.libPaths('~/mambaforge/envs/tidymodels/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_Sept/RF_JAN25')
#### Randomforest Classification
library(tidymodels)
library(GGally)
library(gt)
library(visdat)
library(themis)
library(skimr)
library(viridis)
library(finetune)
library(xgboost)
require(doParallel)
library(kknn)
library(forcats)
library(vip)

#Parallelisation
cores <- 5
cl <- makePSOCKcluster(cores)
registerDoParallel(cl)

dats <- read.table('RF-Input-250107.txt',header=T)
type = args[1]  #e.g. Binary,Regression
vari = args[2] #e.g. Mean,Max
xp = args[3] #subset experiment, CG WGBS HZ
iter = args[4] #iteration

dat = dats %>% filter(variable == vari & Subset == xp)
if (type == 'Regression') {
  cat('Response type is Continuous, taking continuous values\n')
  rf = dat %>% select(c(Continuous,Length,Relativepos,Region,HapD,FuLiD,TajimaD,Dxy,FST,GC))
  rf = rf %>% dplyr::rename(Response = Continuous)
} else {
  cat('Response type is Binary, taking classified values\n')
  rf = dat %>% select(c(as.symbol(type),Length,Relativepos,Region,HapD,FuLiD,TajimaD,Dxy,FST,GC))
  rf = rf %>% dplyr::rename(Response = Classified)
}

#for exploration, run on subset
#rf = rf %>% sample_n(5000)
# convert to numeric
rf = rf %>% mutate_at(grep("Length|Relative|HapD|Fu|Taj|Dxy|FST|GC",colnames(.)),funs(as.numeric))

# convert all remaining character variables to factors
rf = rf %>% mutate(across(where(is.character), as.factor))
skim(rf)

#split into training and testing datasets
data_split <- initial_split(rf, strata = "Response", prop = 0.75)

train_explore <- training(data_split)
test_explore  <- testing(data_split)
data_explore = train_explore

# Generate resamples and repeat
cv_folds = vfold_cv(train_explore, v = 5, strata = Response)

# define the recipe
mc_recipe =  # which consists of the formula (outcome ~ predictors)
  recipe(Response ~ ., data = rf) %>%
  # and some pre-processing steps
  step_novel(all_nominal(), -all_outcomes()) %>% # converts all nominal variables to factors and takes care of other issues related to categorical variables
  step_normalize(Length, Relativepos, HapD, FuLiD, TajimaD, Dxy, FST, GC, -all_outcomes()) %>% # normalizes (center and scales) the numeric variables to have a standard deviation of one and a mean of zero. (i.e., z-standardization)
  step_dummy(all_nominal(), -all_outcomes()) %>% # converts our factor column ocean_proximity into numeric binary (0 and 1) variables
  step_zv(all_numeric(), -all_outcomes()) %>% # removes any numeric variables that have zero variance
  step_corr(all_predictors(), threshold = 0.8, method = "spearman") # will remove predictor variables that have large correlations with other predictor variables

#check the pre-processing
examine_preprocessed <- mc_recipe %>% # apply the recipe to the training data
  prep(rf) %>%  # extract the pre-processed training dataset
  juice()
examine_preprocessed

#set the models
if (type == 'Regression') {
  cat('Response type is Regression, using regression engines\n')
  #boosted trees
  xgb_mod = boost_tree(trees = tune(),min_n = tune(),mtry = tune(),learn_rate = tune()) %>% set_engine("xgboost",importance="permutation") %>% set_mode("regression")
  xgb_wf <- workflow(mc_recipe, xgb_mod)
  xgb_rs <- tune_race_anova(xgb_wf,resamples = cv_folds,grid = 15, metrics = metric_set(rsq),control = control_race(verbose_elim = TRUE))
  #plot_race(xgb_rs)
  show_best(xgb_rs)
  xgb_last = xgb_wf %>% finalize_workflow(select_best(xgb_rs, "rsq")) %>% last_fit(data_split)
  xgb_pred = collect_predictions(xgb_last) %>% select(.pred,Response)
  xgb_pred$Engine = 'XGBoost'
  xgb_vi = extract_workflow(xgb_last) %>% extract_fit_parsnip() %>% vi()
  xgb_vi$Engine = 'XGBoost'
  xgb_fit = xgb_last %>% collect_metrics()
  xgb_fit$Engine = 'XGBoost'
  
  #randomforest
  rf_mod = rand_forest(trees = tune(),min_n = tune(),mtry = tune()) %>% set_engine("ranger",importance="permutation") %>% set_mode("regression")
  rf_wf = workflow(mc_recipe, rf_mod)
  rf_rs = tune_race_anova(rf_wf,resamples = cv_folds,grid = 15, metrics = metric_set(rsq),control = control_race(verbose_elim = TRUE))
  show_best(rf_rs)
  rf_last = rf_wf %>% finalize_workflow(select_best(rf_rs, "rsq")) %>% last_fit(data_split)
  rf_pred = collect_predictions(rf_last) %>% select(.pred,Response)
  rf_pred$Engine = 'Random forest'
  rf_vi = extract_workflow(rf_last) %>% extract_fit_parsnip() %>% vi()
  rf_vi$Engine = 'Random forest'
  rf_fit = rf_last %>% collect_metrics()
  rf_fit$Engine = 'Random forest'
  
} else {
  cat('Response type is Binary, using classification engines\n')
  rf_fit = NULL
  #boosted trees
  xgb_mod = boost_tree(trees = tune(),min_n = tune(),mtry = tune(),learn_rate = tune()) %>% set_engine("xgboost",importance="permutation") %>% set_mode("classification")
  xgb_wf <- workflow(mc_recipe, xgb_mod)
  xgb_rs <- tune_race_anova(xgb_wf,resamples = cv_folds,grid = 15, metrics = metric_set(roc_auc),control = control_race(verbose_elim = TRUE))
  #plot_race(xgb_rs)
  show_best(xgb_rs)
  xgb_last = xgb_wf %>% finalize_workflow(select_best(xgb_rs, "roc_auc")) %>% last_fit(data_split)
  xgb_pred = collect_predictions(xgb_last) %>% select(.pred_class,Response)
  xgb_pred$Engine = 'XGBoost'
  xgb_vi = extract_workflow(xgb_last) %>% extract_fit_parsnip() %>% vi()
  xgb_vi$Engine = 'XGBoost'
  xgb_fit = xgb_last %>% collect_metrics()
  xgb_fit$Engine = 'XGBoost'
  
}

#save fit
fit = rbind(xgb_fit,rf_fit)
names(fit) = c('Metric','Estimator','Estimate','Config','Engine')
fitdf = data.frame(fit,Type = type,Response = vari,Subset = xp,Subset = xp, Iteration=iter)
fitdf
write.table(fitdf,file=paste0('results/',type,'_',vari,'_',xp,'_',iter,'_FIT.txt'),quote=F,sep='\t',row.names=F)

#save vi
vi = rbind(xgb_vi,rf_vi)
vidf = data.frame(vi,Type = type,Response = vari,Subset = xp,Iteration=iter)
vidf
write.table(vidf,file=paste0('results/',type,'_',vari,'_',xp,'_',iter,'_VI.txt'),quote=F,sep='\t',row.names=F)

#save predictions
preds = rbind(xgb_pred,rf_pred)
names(preds) = c('Predicted','Truth','Engine')
preddf = data.frame(preds,Type = type,Response = vari,Subset = xp,Iteration=iter)
#preddf
write.table(preddf,file=paste0('results/',type,'_',vari,'_',xp,'_',iter,'_PRED.txt'),quote=F,sep='\t',row.names=F)

```

Launch:

```bash
for a in $(cat Types.list); do  for b in $(cat Variables.list); do for c in $(cat Subsets.list); do for d in $(cat Iterations.list); do  sbatch -J ${a}_${b}_${c}_${d} randomforest.sh $a $b $c $d;  done;  done; done; done
```

Merge results:

```bash
mergem *VI*txt > Master.vi
mergem *FIT*txt > Master.fit
mergem Class*PRED* > Classification_Predictions.txt
mergem Regre*PRED* > Regression_Predictions.txt
```

### Plot

And plot results:

```R
.libPaths('~/mambaforge/envs/r/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_Sept/RF_JAN25/results')
#### Plot decision tree results
library(tidyverse)
library(viridis)
library(forcats)
library(ggpubr)
library(rstatix)
library(gmodels)
library(RColorBrewer)

#evaluate fit
fit = read.table('Master.fit',header=T,sep='\t') %>% select(-c(Subset.1))
fs = fit %>% group_by(Subset,Response,Type,Metric) %>% 
  summarize(mean = ci(Estimate)[1],
            hi = ci(Estimate)[2],
            lo = ci(Estimate)[3],
            sd = ci(Estimate)[4])
fs = fs %>% 
  mutate(Subset = gsub('CG','ComGar.RRBS',Subset),
         Subset = gsub('WGBS','ComGar.WGBS',Subset),
         Subset = gsub('HZ','HybZon.RRBS',Subset))
fs$Subset <- factor(fs$Subset,levels=c('ComGar.WGBS','ComGar.RRBS','HybZon.RRBS'))
fs %>% filter(Response == 'Mean' & Type == 'Regression' & Metric == 'rsq')
fpr = fs %>% 
  filter(Type == 'Regression') %>% 
  ggplot(aes(x=Response,y=mean,col=Subset,shape=Subset))+
  geom_point(size=3,position=position_dodge(width=0.5))+
  geom_errorbar(aes(ymin=lo,ymax=hi),width=0.5,position=position_dodge(width=0.5))+
  ylab('Evaluation Metric')+xlab('Response (DNA Methylation Divergence Within 5KB Window)')+
  scale_color_manual(values=brewer.pal(3,'Set2'))+
  ggtitle('Regression')+
  theme_bw()+
  theme(legend.position='none')+
  facet_grid(Metric~Type,scales='free')
fpr
fpc = fs %>% 
  filter(Type != 'Regression') %>% 
  ggplot(aes(x=Response,y=mean,col=Subset,shape=Subset))+
  geom_point(size=3,position=position_dodge(width=0.5))+
  geom_errorbar(aes(ymin=lo,ymax=hi),width=0.5,position=position_dodge(width=0.5))+
  ylab('Evaluation Metric')+xlab('Response (DNA Methylation Divergence Within 5KB Window)')+
  scale_color_manual(values=brewer.pal(3,'Set2'))+
  ggtitle('Classification')+
  theme_bw()+
  #theme(legend.position='none')+
  facet_grid(Metric~Type,scales='free')
fpc
fpa = ggarrange(fpr,fpc,widths = c(0.75,1.2))
fpa

#save
pdf('20250107_RF_Fit.pdf',height=5,width=8)
fpa
dev.off()

#evalute VI
vi = read.table('Master.vi',header=T,sep='\t')
#this won't change anything for XGBoost, but to plot together we will scale the RF covariate importance to sum to 1 (as with XGBoost)
vi = vi %>% group_by(Subset,Response,Engine,Type,Iteration) %>% mutate(Importance = ifelse(Importance < 0,0,Importance),
                                                                       TotalImp = sum(Importance),
                                                                       PctImp = Importance / TotalImp,
                                                                       Variable = gsub('GC','GC-Content',Variable),
                                                                       Variable = gsub('Length','Chromosome Length',Variable),
                                                                       Variable = gsub('Relativepos','Chromosomal Position',Variable),
                                                                       Variable = gsub('HapD','Haplotype Diversity',Variable),
                                                                       Variable = gsub('TajimaD',"Tajima's D",Variable),
                                                                       Variable = gsub("FuLiD","Fu & Li's D*",Variable))
vis = vi %>% group_by(Subset,Response,Type,Variable) %>% 
  filter(Response == 'Mean') %>% 
  summarize(mean = ci(PctImp)[1],
            hi = ci(PctImp)[2],
            lo = ci(PctImp)[3],
            sd = ci(PctImp)[4])
vis = vis %>% 
  mutate(Subset = gsub('CG','ComGar.RRBS',Subset),
         Subset = gsub('WGBS','ComGar.WGBS',Subset),
         Subset = gsub('HZ','HybZon.RRBS',Subset))
vis$Subset <- factor(vis$Subset,levels=c('ComGar.WGBS','ComGar.RRBS','HybZon.RRBS'))
#find most important overall, we will sort by this
imps = vi %>% filter(Response == 'Mean' & Type == 'Regression') %>% group_by(Response,Type,Variable) %>% summarize(mean=mean(PctImp)) %>% arrange(desc(mean)) %>% mutate(ORDER = row_number())
vis = left_join(vis,imps %>% select(Response,Type,Variable,ORDER))
vip = vis %>% 
  filter(Response == 'Mean' & Type == 'Regression') %>% 
  mutate(NewOrd = fct_reorder(Variable,ORDER,.desc = TRUE)) %>% 
  ggplot(aes(y=NewOrd,x=mean,fill=Subset))+
  geom_bar(stat='identity',position=position_dodge(width=0.9))+
  geom_errorbar(aes(xmin=ifelse(mean-sd<0,0,mean-sd),xmax=mean+sd),width=0.5,position=position_dodge(width=0.9))+
  scale_fill_manual(values=brewer.pal(3,'Set2'))+
  theme_bw()+
  #theme(legend.position='none')+
  ylab('')+xlab('Permutation Importance')+
  #facet_grid(.~Response,scales='free')+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))
vip

viS = vi %>% group_by(Subset,Response,Type,Variable,Engine) %>% 
  #filter(Response == 'Mean' & Type == 'Regression') %>% 
  summarize(mean = ci(PctImp)[1],
            lo = ci(PctImp)[2],
            hi = ci(PctImp)[3],
            se = ci(PctImp)[4])
write.table(viS,file='20250107_Permutation_Importance.txt',quote=F,sep='\t',row.names=F)

# Across all experiments, as reported in main text
vi %>% group_by(Variable) %>% 
  filter(Response == 'Mean' & Type == 'Regression') %>% 
  summarize(mean = ci(PctImp)[1],
            lo = ci(PctImp)[2],
            hi = ci(PctImp)[3],
            se = ci(PctImp)[4])

# Variable               mean     lo     hi      se
# <chr>                 <dbl>  <dbl>  <dbl>   <dbl>
#   1 Chromosomal Position 0.0542 0.0478 0.0605 0.00320
# 2 Chromosome Length    0.0731 0.0670 0.0792 0.00307
# 3 Dxy                  0.104  0.0919 0.116  0.00612
# 4 FST                  0.0238 0.0202 0.0274 0.00182
# 5 Fu & Li's D*         0.0558 0.0514 0.0601 0.00218
#  6 GC-Content           0.284  0.260  0.308  0.0121 
#  7 Haplotype Diversity  0.120  0.109  0.132  0.00598
#  8 Region_Intergenic    0.0315 0.0272 0.0357 0.00214
#  9 Region_Intron        0.0251 0.0210 0.0292 0.00207
# 10 Region_Promoter      0.145  0.111  0.178  0.0168 
# 11 Region_Repeat        0.0168 0.0122 0.0215 0.00231
# 12 Tajima's D           0.0674 0.0616 0.0732 0.00291

#save
pdf('20250107_RF_VIP_Regression.pdf',height=6,width=6)
vip 
dev.off()


#Also plot for ALL subsets, experiments, types 
ALL = vi %>% group_by(Subset,Response,Type,Variable) %>% 
  summarize(mean = ci(PctImp)[1],
            hi = ci(PctImp)[2],
            lo = ci(PctImp)[3],
            sd = ci(PctImp)[4])
ALL = ALL %>% 
  mutate(Subset = gsub('CG','ComGar.RRBS',Subset),
         Subset = gsub('WGBS','ComGar.WGBS',Subset),
         Subset = gsub('HZ','HybZon.RRBS',Subset))
ALL$Subset <- factor(ALL$Subset,levels=c('ComGar.WGBS','ComGar.RRBS','HybZon.RRBS'))
#find most important overall, we will sort by this
ALL = left_join(ALL,imps %>% ungroup %>% select(Variable,ORDER))
ALL$Type = factor(ALL$Type,levels=c('Regression','Classification'))
ALLp = ALL %>% 
  mutate(NewOrd = fct_reorder(Variable,ORDER,.desc = TRUE)) %>% 
  ggplot(aes(y=NewOrd,x=mean,fill=Subset))+
  geom_bar(stat='identity',position=position_dodge(width=0.9))+
  geom_errorbar(aes(xmin=ifelse(mean-sd<0,0,mean-sd),xmax=mean+sd),width=0.5,position=position_dodge(width=0.9))+
  scale_fill_manual(values=brewer.pal(3,'Set2'))+
  theme_bw()+
  #theme(legend.position='none')+
  ylab('')+xlab('Permutation Importance')+
  facet_grid(Response~Type,scales='free')+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))
ALLp

#save
pdf('20250107_RF_VIP_All.pdf',height=6,width=10)
ALLp
dev.off()

#evalute predictions
conf = read.table('Regression_Predictions.txt',header=T,sep='\t')
cp = conf %>% 
  ggplot(aes(x=Predicted,y=Truth,col=Subset))+
  #geom_point(alpha = 0.5)+
  geom_smooth()+
  #geom_density_2d()+
  scale_color_manual(values=c(viridis(3)))+
  theme_bw()+xlab('Predicted')+ylab('Truth')+
  theme(legend.position = "none")+
  facet_wrap(Subset~Response,scales='free',nrow=3)
cp

#evalute predictions, classification
confclass = read_tsv('Classification_Predictions.txt')
confclass %>% group_by(Predicted,Truth,Engine,Response,Subset) %>% 
  summarize(count = n()) %>% 
  ggplot(aes(x=Predicted,y=Truth,fill=count))+
  geom_tile()+
  scale_fill_continuous(low='yellow',high='red')+
  theme_bw()+xlab('Predicted')+ylab('Truth')+
  facet_grid(Response+Engine~Subset,scales='free')


#save
png('20250107_RF_Predictions.png',height=5,width=3.5,res=600,units='in')
cp
dev.off()

```

### Rho: 5mC & Population Genetic

For genetic v 5mC correlations using the RF input (within 5kb windows):

```R
.libPaths('~/mambaforge/envs/r/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_Sept/RF_JAN25')
#### PopGen ~ 5mC spearman correlations
library(tidyverse)
library(viridis)
library(forcats)
library(ggpubr)
library(rstatix)
library(gmodels)
library(RColorBrewer)

#plot some of the interesting comparisons
rf = read_tsv('RF-Input-250107.txt')

### Correlations: 5mC / genetic variation
cat('Working on correlations for genetic variation, 5mc \n')
rf1 = rf %>% filter(variable == 'Mean') %>% select(chr,GenStart,GenEnd,raw,HapD,FuLiD,TajimaD,Dxy,FST,Subset)
rf1 = rf1 %>% mutate(Focal_Region = ifelse(chr == '18' & GenStart > 8070001 & GenStart < 10070000,'Focal','Undiverged'))
rf2 = rf1 %>%
  pivot_longer(!c(chr,GenStart,GenEnd,raw,Subset,Focal_Region)) %>%
  na.omit()
xps = c('CG','HZ','WGBS')

fullboots = NULL
options(dplyr.summarise.inform = FALSE)
set.seed(123)
B <- 999 # number of bootstrap replicates

for (xp in xps) {
  for (vr in unique(rf2$name)) {
    cat('Working on Experiment & Variable: ',xp,' ',vr,'\n')
    bd = rf2 %>% filter(Subset == xp & name == vr)
    n = bd %>% filter(Focal_Region == 'Focal') %>% nrow
    bdat=NULL
    ## Run the bootstrap procedure:
    for (b in 1:B) {
      cat('Iteration: ',b,'\n')
      d1 = bd %>% group_by(Focal_Region) %>% slice_sample(n=n,replace=TRUE) %>% summarize(cor = cor(raw,value,method='spearman'))
      bdat = rbind(bdat,d1 )
    }
    bdat$Experiment = xp
    bdat$Variable = vr
    names(bdat) = c('Focal_Region','Boots','Experiment','Variable')
    fullboots = rbind(fullboots,bdat)
  }
}

write.table(fullboots,file=paste0('Genetic5mC_Bootstraped_Results-20250107.txt'),quote=F,sep='\t',row.names=F)
fullboots = read.table('Genetic5mC_Bootstraped_Results-20250107.txt',header=TRUE)

#add counts for windows 
features = rf2 %>% group_by(Subset,Focal_Region,name) %>% count()
names(features)= c('Experiment','Focal_Region','Variable')

cdm = fullboots %>% 
  mutate(Experiment = gsub('CG','ComGar.RRBS',Experiment),
         Experiment = gsub('WGBS','ComGar.WGBS',Experiment),
         Experiment = gsub('HZ','HybZon.RRBS',Experiment))
cdm$Experiment <- factor(cdm$Experiment,levels=c('HybZon.RRBS','ComGar.RRBS','ComGar.WGBS'))
#calculate CIs
cis = cdm %>% 
  group_by(Experiment,Variable,Focal_Region) %>% 
  summarize(LowerBound = quantile(Boots,0.025), #can change based on your significance threshold 
            UpperBound = quantile(Boots,0.975))
cis = cis %>% group_by(Experiment,Variable) %>% 
  mutate(Signif = ifelse(LowerBound > 0 | UpperBound < 0,'*','n.s.'))

genp = cdm %>% 
  ggplot(aes(x=Variable,y=Boots,fill=Focal_Region))+
  geom_boxplot(alpha=0.75)+
  geom_text(data=cis,aes(x=Variable,col=Focal_Region,group=Focal_Region,y=Inf,label=Signif,size=Signif),col='black',vjust=1.5,position=position_dodge(width=1))+
  scale_size_manual(values=c(6,4))+
  scale_fill_manual(values=c('darkorchid2','grey60'))+
  facet_grid(Experiment~.,scales='free')+
  geom_hline(yintercept=0,lty=2)+
  ylab("Spearman's Rho: Taxon Methylation Divergence") + xlab('')+
  scale_y_continuous(expand = expansion(mult = .25)) + #expand y axis slightly 
  theme_bw()
genp

pdf('../../RF_GenCorBootBootstraps.pdf',height=4,width=5)
genp
dev.off()

write.table(cis,file='../../Correlations_95CIs.txt',quote=F,sep='\t',row.names=F)
```

# Focal v BG Bootstraps

## Bootstrap

Bootstrap sampling compare DMP test statistics within and outside focal region

```bash
# Bootstrap distributions between 5mC test statistics in focal region vs genomic background
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_202504')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)

set.seed(123)
# Read file
d2 = read.table('5mC_10x_STAT-Region_202504.txt',header=TRUE)
d3 <- d2 %>% 
  pivot_longer(matches('mstat'),values_to = 'Stat', names_to = 'experiment') %>% 
  mutate(experiment = toupper(gsub('mstat.','',experiment)),
         chr = gsub('chr','',chr)) %>% 
  mutate(Focal_Region = ifelse(chr == '18' & start > 8070001 & start < 10070000,'Focal','Undiverged'), #our focal region lies on chr18 between these coordinates
         Stat = abs(Stat)) %>% na.omit %>% filter(Region != 'Repeat') #removing any NAs and we won't analyze repeats because there are so few 

# #store your bootstrapped data in this object 
# fullboots = NULL
# options(dplyr.summarise.inform = FALSE)
# 
# 
# bdat=NULL
# for (xp in unique(d3$experiment)) {
#   
# subd3 <- d3 %>% filter(experiment == xp)  
d4 = d3 %>% mutate(Group = paste0(experiment,'_',Region,'_',Focal_Region)) %>% select(Group,Stat)
groupNa = d3 %>% filter(Focal_Region == 'Focal') %>% count(Focal_Region,experiment,Region)
groupNb = d3 %>% filter(Focal_Region == 'Focal') %>% count(Focal_Region,experiment,Region) %>% mutate(Focal_Region = gsub('Focal','Undiverged',Focal_Region))
groupN = rbind(groupNa,groupNb) %>% mutate(Group = paste0(experiment,'_',Region,'_',Focal_Region)) %>% ungroup %>% select(Group,n) %>% arrange(Group)
d4 %>% 
  group_split(Group) %>% 
  map2_dfr(c(groupN$n), ~ slice_sample(.x, n = .y,replace=TRUE)) %>% 
  group_by(Group) %>% 
  summarize(Mean = mean(Stat,na.rm=TRUE),
            Median = median(Stat,na.rm=TRUE),
            Count = n())
bdat=NULL
B <- 999 # number of bootstrap replicates, realistically at least 999

## Run the bootstrap procedure:
for (b in 1:B) {
  cat('Iteration: ',b,'\n')
  groupN = rbind(groupNa,groupNb) %>% mutate(Group = paste0(experiment,'_',Region,'_',Focal_Region)) %>% ungroup %>% select(Group,n) %>% arrange(Group)
  d5 = d4 %>% 
    group_split(Group) %>% 
    map2_dfr(c(groupN$n), ~ slice_sample(.x, n = .y,replace=TRUE)) %>% 
    group_by(Group) %>% 
    summarize(Mean = mean(Stat,na.rm=TRUE),
              Count = n())
  bdat = rbind(bdat,d5) #calculate the difference in means between the 2 comparisons 
}

#save 
write.table(bdat,file=paste0('5mC_Bootstraped_ResultsBoxplots-20250414.txt'),quote=F,sep='\t',row.names=F)
bdat = read.table('5mC_Bootstraped_ResultsBoxplots-20250414.txt',header=TRUE)

#fix ugly names
bdatplot = bdat %>% separate(Group, into=c('Experiment','Region','Focal_Region')) %>% 
  mutate(Experiment = gsub('CG','ComGar',Experiment),
         Experiment = gsub('HZ','HybZon',Experiment)) 
#proper order
bdatplot$Experiment <- factor(bdatplot$Experiment,levels=c('HybZon','ComGar'))
bdatplot$Region = factor(bdatplot$Region,levels=c('Promoter','CDS','Intron','Intergenic'))
#calculate CIs
cis = bdatplot %>% 
  group_by(Experiment,Region,Focal_Region,Count) %>% 
  summarize(LowerBound = quantile(Mean,0.025), #can change based on your significance threshold 
            UpperBound = quantile(Mean,0.975))
cis = cis %>% group_by(Experiment,Region) %>% 
  mutate(Signif = ifelse(LowerBound[Focal_Region == 'Focal'] > UpperBound[Focal_Region == 'Undiverged'],'*',
                         ifelse(UpperBound[Focal_Region == 'Focal'] < LowerBound[Focal_Region == 'Undiverged'],'*','n.s.')))

divplot = bdatplot %>% 
  ggplot(aes(x=Region,y=Mean,fill=Focal_Region))+
  geom_boxplot(alpha=0.8,outlier.alpha = 0.25)+
  geom_text(data=cis,aes(x=Region,group=Region,y=Inf,label=Signif,size=Signif),vjust=1.5,position=position_dodge(width=1))+
  geom_text(data=cis %>% filter(Focal_Region == 'Focal'),aes(x=Region,group=Region,y=-Inf,label=Count),vjust=-1,size=3,position=position_dodge(width=1))+
  scale_size_manual(values=c(6,4))+
  scale_fill_manual(values=c('darkorchid2','grey60'))+
  facet_grid(Experiment~.,scales='free')+
  ylab('Bootstrapped DNA Methylation Divergence') + xlab('')+
  scale_y_continuous(expand = expansion(mult = .25)) + #expand y axis slightly 
  theme_bw(base_size=10)
divplot

pdf('~/symlinks/hzone/Methylation_Analyses/figures/20250414_Bootstrapped_MethylationDivergenceBoxplot.pdf',height=4,width=4)
divplot
dev.off()

write.table(cis,file='~/symlinks/hzone/Methylation_Analyses/figures/MethylationDivergence_95CIs.txt',quote=F,sep='\t',row.names=F)

```

## Including ComGar.WGBS

```R
# Bootstrap distributions between 5mC test statistics in focal region vs genomic background
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_202504')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)

set.seed(123)
# Read file
d2 = read.table('5mC_10x_STAT-Region_202504.txt',header=TRUE)
d3 <- d2 %>% 
  pivot_longer(matches('mstat'),values_to = 'Stat', names_to = 'experiment') %>% 
  mutate(experiment = toupper(gsub('mstat.','',experiment)),
         chr = gsub('chr','',chr)) %>% 
  mutate(Focal_Region = ifelse(chr == '18' & start > 8070001 & start < 10070000,'Focal','Undiverged'), #our focal region lies on chr18 between these coordinates
         Stat = abs(Stat)) %>% na.omit %>% filter(Region != 'Repeat') #removing any NAs and we won't analyze repeats because there are so few 

# #store your bootstrapped data in this object 
# fullboots = NULL
# options(dplyr.summarise.inform = FALSE)
# 
# 
# bdat=NULL
# for (xp in unique(d3$experiment)) {
#   
# subd3 <- d3 %>% filter(experiment == xp)  
d4 = d3 %>% mutate(Group = paste0(experiment,'_',Region,'_',Focal_Region)) %>% select(Group,Stat)
groupNa = d3 %>% filter(Focal_Region == 'Focal') %>% count(Focal_Region,experiment,Region)
groupNb = d3 %>% filter(Focal_Region == 'Focal') %>% count(Focal_Region,experiment,Region) %>% mutate(Focal_Region = gsub('Focal','Undiverged',Focal_Region))
groupN = rbind(groupNa,groupNb) %>% mutate(Group = paste0(experiment,'_',Region,'_',Focal_Region)) %>% ungroup %>% select(Group,n) %>% arrange(Group)
d4 %>% 
  group_split(Group) %>% 
  map2_dfr(c(groupN$n), ~ slice_sample(.x, n = .y,replace=TRUE)) %>% 
  group_by(Group) %>% 
  summarize(Mean = mean(Stat,na.rm=TRUE),
            Median = median(Stat,na.rm=TRUE),
            Count = n())
bdat=NULL
B <- 999 # number of bootstrap replicates, realistically at least 999

## Run the bootstrap procedure:
for (b in 1:B) {
  cat('Iteration: ',b,'\n')
  groupN = rbind(groupNa,groupNb) %>% mutate(Group = paste0(experiment,'_',Region,'_',Focal_Region)) %>% ungroup %>% select(Group,n) %>% arrange(Group)
  d5 = d4 %>% 
    group_split(Group) %>% 
    map2_dfr(c(groupN$n), ~ slice_sample(.x, n = .y,replace=TRUE)) %>% 
    group_by(Group) %>% 
    summarize(Mean = mean(Stat,na.rm=TRUE),
              Count = n())
  bdat = rbind(bdat,d5) #calculate the difference in means between the 2 comparisons 
}

#save 
write.table(bdat,file=paste0('5mC_Bootstraped_ResultsBoxplots-20250107.txt'),quote=F,sep='\t',row.names=F)
bdat = read.table('5mC_Bootstraped_ResultsBoxplots-20250107.txt',header=TRUE)

#fix ugly ass names
bdatplot = bdat %>% separate(Group, into=c('Experiment','Region','Focal_Region')) %>% 
  mutate(Experiment = gsub('CG','ComGar.RRBS',Experiment),
         Experiment = gsub('WGBS','ComGar.WGBS',Experiment),
         Experiment = gsub('HZ','HybZon.RRBS',Experiment)) 
#proper order
bdatplot$Experiment <- factor(bdatplot$Experiment,levels=c('HybZon.RRBS','ComGar.RRBS','ComGar.WGBS'))
bdatplot$Region = factor(bdatplot$Region,levels=c('Promoter','CDS','Intron','Intergenic'))
#calculate CIs
cis = bdatplot %>% 
  group_by(Experiment,Region,Focal_Region,Count) %>% 
  summarize(LowerBound = quantile(Mean,0.025), #can change based on your significance threshold 
            UpperBound = quantile(Mean,0.975))
cis = cis %>% group_by(Experiment,Region) %>% 
  mutate(Signif = ifelse(LowerBound[Focal_Region == 'Focal'] > UpperBound[Focal_Region == 'Undiverged'],'*',
                         ifelse(UpperBound[Focal_Region == 'Focal'] < LowerBound[Focal_Region == 'Undiverged'],'*','n.s.')))

divplot = bdatplot %>% 
  ggplot(aes(x=Region,y=Mean,fill=Focal_Region))+
  geom_boxplot(alpha=0.8,outlier.alpha = 0.25)+
  geom_text(data=cis,aes(x=Region,group=Region,y=Inf,label=Signif,size=Signif),vjust=1.5,position=position_dodge(width=1))+
  geom_text(data=cis %>% filter(Focal_Region == 'Focal'),aes(x=Region,group=Region,y=-Inf,label=Count),vjust=-1,size=3,position=position_dodge(width=1))+
  scale_size_manual(values=c(6,4))+
  scale_fill_manual(values=c('darkorchid2','grey60'))+
  facet_grid(Experiment~.,scales='free')+
  ylab('Bootstrapped DNA Methylation Divergence') + xlab('')+
  scale_y_continuous(expand = expansion(mult = .25)) + #expand y axis slightly 
  theme_bw()
divplot

pdf('20250107_Bootstrapped_MethylationDivergenceBoxplot.pdf',height=5,width=4)
divplot
dev.off()

write.table(cis,file='MethylationDivergence_95CIs.txt',quote=F,sep='\t',row.names=F)

```

# Ordinations

## Methylation Ordinations

5mC: dbRDA 

```bash
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

```

### Including ComGar.WGBS

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



# Manhattan Plots

## Taxon 

Manhattan plot, showing the FDR-corrected pvalues for 5mC and pop gen data:

```R
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
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_202504')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(karyoploteR)
library(viridis)
library(LICORS)
library(ggplot2)
library(scales)
library(RColorBrewer)

class = read.table('Full-HZ-CG.DSS.BP-10x.Classified-Annotated-SNPless_20250410.txt',header=T,sep = "\t",comment.char = '')#volcano
v1 = class %>% select(site,Class,contains(c('divergence'),ignore.case = FALSE))  %>% pivot_longer(!c(site,Class),values_to = 'Divergence') %>% na.omit %>% mutate(name = gsub('divergence_','',name))
v1 = v1 %>% separate(name,into=c('Variable','Experiment'))
v2 = class %>% select(site,Class,contains(c('fdrs'),ignore.case = FALSE))  %>% pivot_longer(!c(site,Class),values_to = 'FDRS') %>% na.omit %>% mutate(name = gsub('fdrs_','',name))
v2 = v2 %>% separate(name,into=c('Variable','Experiment'))
v3 = left_join(v1,v2)
v4 <- v3 %>% filter(!grepl('Hybrid|hz',Experiment)) %>% na.omit
v4 = v4 %>% mutate(Experiment = gsub('^CG','ComGar',Experiment),
                   Variable = gsub('Species','Taxon',Variable),
                   Variable = gsub('StageCHK','StageChick',Variable),
                   Variable = gsub('StageYRL','StageYearling',Variable))

# Checkpoint
write.table(v4,file='20250415_Volcano_Input.txt',quote=F,sep='\t',row.names=F)
v4 = read.table('20250415_Volcano_Input.txt',header=T,sep = "\t",comment.char = '')

# Plot
png('~/symlinks/hzone/Methylation_Analyses/figures/20250415_VolcanoPlots-Known.png',units='in',height=5,width=7,res=600)
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
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Methylation_Analyses/Classification_202504')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(ggplot2)
library(scales)
library(RColorBrewer)

# Import 5mC data
class <- read.table('Full-HZ-CG.DSS.BP-10x.Classified-Annotated-SNPless_20250410.txt',header=TRUE)

##### Plot Hybrid Index ~ 5mC #####
#Plot Sites of Interest 
md = read.table('HZ_Metadata.txt',header=T)
hz = class %>% filter(Class == 'Taxon' & SupportCount >= 2 & DMP.HZ == 'Taxon') %>% select(site,contains('.HZ'))

#first hybrid index
hi = hz %>% 
  select(-c(DMP.HZ,Divergence.HZ,contains('divergence'),contains('fdrs'))) %>%
  pivot_longer(!site,names_to = 'Short_ID',values_to = 'Methylation') %>% 
  mutate(Short_ID = gsub('.HZ','',Short_ID))
him = merge(hi,md %>% select(Short_ID = variable,Chr18_Hybrid_Index,Species),by='Short_ID')
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

ggsave('~/symlinks/hzone/Methylation_Analyses/figures/20250414_Hybrid5mC_Distributions.pdf',hip,dpi=300,height=2,width=4)


sid=NULL
for (si in unique(him$site)) {
  s1 = him %>% filter(site == si) %>% select(id = Short_ID, mC = Methylation) %>% na.omit()
  mcd = s1 %>% dplyr::select(mC) %>% data.frame
  rownames(mcd) = s1$id
  
  #kmean's clustering
  k2 = kmeans(mcd, centers = 3, nstart = 5)
  mcd$cluster = k2$cluster
  mcd$Short_ID = rownames(mcd)
  mcd$site = si
  
  sid = rbind(sid,mcd)
}

#plot
hild = merge(sid,md %>% dplyr::select(Short_ID = variable,Chr18_Hybrid_Index,Species),by='Short_ID')

# Assign genotypes based on unmeth / het / meth
clusters <- hild %>%
  group_by(site,cluster) %>% 
  summarize(mean = mean(mC)) %>% 
  ungroup %>% 
  group_by(site) %>%
  arrange(mean, .by_group = TRUE) %>% # Sort within each group by mean
  mutate(Genotype = case_when(
    row_number() == 1 ~ "U/U", # Lowest mean
    row_number() == 2 ~ "U/M", # Middle mean
    row_number() == 3 ~ "M/M"  # Highest mean
  )) %>%
  ungroup() %>% select(site,cluster,Genotype)
him2 <- left_join(hild,clusters)

# Plot and confirm genotypes are consistent
hip = him2 %>% dplyr::rename(Taxon = Species) %>% 
  ggplot(aes(x=Chr18_Hybrid_Index,y=mC,fill=Genotype,shape=Taxon))+
  geom_jitter(width=0.00,size=2.5,show.legend = T)+
  scale_fill_manual(values=rev(viridis(4,option='mako')))+
  scale_shape_manual(values=c(24,21,22))+ylim(c(-.01,1.01))+
  xlab('Genetic Hybrid Index')+ylab('Methylation (%)')+
  facet_wrap(site~.,scales='free')+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  guides(fill=guide_legend(override.aes=list(shape=21)))
hip

pdf('~/symlinks/hzone/Methylation_Analyses/figures/20250414_Clustered_Genotypes.pdf',height=2,width=5)
hip
dev.off()

##### LD Between outlier sites ##### 
counter = 0
for (si in unique(him2$site)) {
  counter = counter + 1
  gts = subset(him2,site==si) %>% dplyr::select(Genotype)
  gt = genotype(gts$Genotype)
  assign(paste0('g',counter),gt)
}

LD(g1,g2)

# Pairwise LD
# -----------
#   D        D'      Corr
# Estimates: 0.1362936 0.9994866 0.6321308
# 
#               X^2      P-value  N
# LD Test: 17.58193 2.751904e-05 22

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



# Sensitivity Analyses

## Library Effects: ComGar.WGBS

Examine WGBS inter-library variation:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=24:00:00

SCRATCH=/tmp/$SLURM_JOB_ID

genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Genomes/Corvus_ASM73873v5
genome=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.fa
cgi=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.CGI.bed
repeats=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.REPEATS.bed
cutsites=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.CUT-SITES.bed
ctsnp=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.TS_SNPs.bed
chrsbed=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.CHRS.bed

#working directories
rawdataWGBS=/dss/dsslegfs01/pr53da/pr53da-dss-0018/rawdata/Bisulfite-Seq/WGBS
rawdataRRBS=/dss/dsslegfs01/pr53da/pr53da-dss-0018/rawdata/Bisulfite-Seq/RRBS
basedir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/crow_scratch/2021-07_HybridZone
outdir=${basedir}/scratch/wgbs/lib_check
workdir=${basedir}/scratch/wgbs/lib_check

#positional arguments
RUN=$1

mkdir ${workdir}
mkdir ${outdir}

cat ${rawdataWGBS}/${RUN}.fastq.gz > $SCRATCH/${RUN}.R1.fastq.gz
cat ${rawdataWGBS}${RUN}.fastq.gz > $SCRATCH/${RUN}.R2.fastq.gz

cd $workdir

#Trim adapters
trim_galore --fastqc -j 10 --clip_r1 9 --clip_r2 9 --three_prime_clip_R2 1 --paired --quality 20 --output_dir ${workdir} ${SCRATCH}/${RUN}.R1.fastq.gz ${SCRATCH}/${RUN}.R2.fastq.gz

#Map reads
bismark --parallel 7 --output_dir ${workdir} --genome ${genomedir} -1 ${workdir}/${RUN}.R1_val_1.fq.gz -2 ${workdir}/${RUN}.R2_val_2.fq.gz

#Deduplicate
#deduplicate_bismark --paired ${workdir}/${RUN}.R1_val_1_bismark_bt2_pe.bam --output_dir ${workdir}

#Extract methylation
bismark_methylation_extractor --parallel 7 --gzip --bedgraph --ignore_3prime_r2 1 --ignore_r2 3 --buffer_size 20G --cytosine_report --genome_folder ${genomedir} --output ${workdir} ${workdir}/${RUN}.R1_val_1_bismark_bt2_pe.deduplicated.bam
```

Subset reads:

```bash
#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=2:00:00

genomedir=/dss/dsslegfs01/pr53da/pr53da-dss-0018/projects/2021__HybridZone_Epigenetics/Genomes/Corvus_ASM73873v5
genome=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.fa
cgi=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.CGI.bed
repeats=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.REPEATS.bed
cutsites=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.CUT-SITES.bed
ctsnp=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.TS_SNPs.bed
chrsbed=${genomedir}/GCF_000738735.5_ASM73873v5_genomic.CHRS.bed

#working directories
basedir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/crow_scratch/2021-07_HybridZone
outdir=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/crow_scratch/2021-07_HybridZone/scratch/wgbs/lib_check
workdir=${basedir}/scratch/wgbs

#positional arguments
RUN=$1

#conver to temporary file
zcat ${workdir}/${RUN}*deduplicated.bismark.cov.gz | awk '{OFS="\t"}{print $1, $2, $3, $4, $5, $6}' | sort | uniq > ${workdir}/${RUN}.tmp

#Remove any positions overlapping a transition
bedtools subtract -A -a ${workdir}/${RUN}.tmp -b ${ctsnp} > ${workdir}/${RUN}.filt_ctsnp.tmp

#Only keep sites on the major chromosomes
bedtools intersect -wb -a ${chrsbed} -b ${workdir}/${RUN}.filt_ctsnp.tmp | awk '{OFS="\t"}{print $4, $5, $6, $7, $8, $9}' | sort | uniq | bgzip -c > ${outdir}/${RUN}.CpG_5mC.cov.gz
```

Average:

```R
library(dplyr)
library(tidyr)
library(reshape2)
library(matrixStats)

#grab all files
files <- list.files(path=".", pattern="*.cov.gz", full.names=TRUE, recursive=FALSE)
allnames <- NULL

#loop through them..
for (samp in files) {

  cat('Reading in file: ',samp,'\n')
  #read data
  ztab <- gzfile(samp,'rt');
  tab <- read.table(ztab,header=F);
  #only keep positions with more than 4x coverage
  tab <- subset(tab, (V5 + V6) > 4)
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

#filter according to missingness, if there's 3 more more missing individuals, drop the site
covdat <- masterWGBS[rowSums(is.na(masterWGBS[grepl('^Cs_', names(masterWGBS))])) == 0, ]

#and filter for maximum coverage:
covs <- covdat[grepl('^site|^Cs_|^Ts_',names(covdat))]
covs$total <- rowSums(covs[grepl('Ts_|Cs_',names(covs))],na.rm=TRUE)
thresh <- quantile(covs$total,0.999)
keep <- subset(covs,total < thresh)
keep <- keep[grepl('^site$',names(keep))]

#calculate variable divergence
p5mC <- masterWGBS[grepl('^D_|^H_|^S_',names(masterWGBS))]
masterWGBS$Species <- (rowMeans(p5mC[grepl('^D_', names(p5mC))],na.rm=TRUE) - rowMeans(p5mC[grepl('^S_', names(p5mC))],na.rm=TRUE))
masterWGBS$Tissue <- (rowMeans(p5mC[grepl('_BL_', names(p5mC))],na.rm=TRUE) - rowMeans(p5mC[grepl('_M_', names(p5mC))],na.rm=TRUE))
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

#362,655 windows here, equates to 181,327,BP.p, or around 16% of the genome
masterWGBS <- masterWGBSf[!grepl('^minsite|^Sites_',names(masterWGBSf))]

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
write.table(masterWGBS,'WGBS_PCA_Input.txt',quote=F,row.names=FALSE)

#now ready for model
write.table(CT ,file="WGBS_GLMM.txt",row.names=FALSE,quote=F,col.names=TRUE)
```

And plot:

```R
library(dplyr)
library(tidyr)
library(reshape2)
library(matrixStats)
library(ggplot2)

pc1 <- read.table('WGBS_PCA_Input.txt',header=T)

#pca-rda...
#remove missing data
dats <- na.omit(pc1[grepl('^S_|^D_',names(pc1))])
dats$var <- apply(dats, 1, var)
dats2 <- subset(dats,var > 0)
dats3 <- t(dats2[!grepl('var',names(dats2))])
pca <- prcomp(dats3,scale=TRUE)
var_explained <- pca$sdev^2/sum(pca$sdev^2)
pcd <- data.frame(pca$x)
pcd$ID <- row.names(pcd)

write.table(var_explained,file='PCA_Eigenvals.txt',quote=F,row.names=F)
write.table(pcd,file='PCA_Eigenvecs.txt',quote=F,row.names=F)

###if necessary, download the var/pcd frames
setwd('F:/Research/scratch/crow_hybrid_paper/2021-11/PCA')

library(tidyverse)
library(ggplot2)
library(ggsci)
library(gridExtra)
library(patchwork)

vec <- read.table('PCA_Eigenvecs.txt',header=T)
val <- read.table('PCA_Eigenvals.txt',header=T)

vec$BAK <- vec$ID
vec <- vec %>% separate(BAK, into = c('Species','Locality','IDNUM','Tissue','Stage','Sex','Batch'))
vec <- vec %>% replace_na(list(Batch='SciLife'))
vec$Batch <- gsub('NO','Novogene',vec$Batch)

p12 <- ggplot(vec,aes(x=PC1,y=PC2,col=Tissue,shape=Batch,size=Sex,label=IDNUM))+
  geom_point(show.legend=T)+  
  geom_text(col='black',size=4)+
  scale_color_jco()+
  xlab(paste0("PC1 (", signif(val$x[1], 4)*100, "%)")) +
  ylab(paste0("PC2 (", signif(val$x[2], 4)*100, "%)")) +
  scale_size_manual(values=c(4,8))+
  theme_bw(base_size=14)
p34 <- ggplot(vec,aes(x=PC3,y=PC4,col=Tissue,shape=Batch,size=Sex,label=IDNUM))+
  geom_point(show.legend=T)+  
  geom_text(col='black',size=4)+
  scale_color_jco()+
  xlab(paste0("PC1 (", signif(val$x[3], 4)*100, "%)")) +
  ylab(paste0("PC2 (", signif(val$x[4], 4)*100, "%)")) +
  scale_size_manual(values=c(4,8))+
  theme_bw(base_size=14)

combined <- p12 + p34 & theme(legend.position = "right")
combined <- combined + plot_layout(guides = "collect")
combined

png('WGBS-Batch-Effects.png',units='in',height=4,width=10,res=600,bg='transparent')
combined
dev.off()
```



## HybZon & ComGar.RRBS Batch

Try only chicks males from the common garden, to the PCA stage:

```r
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)

ztab <- gzfile('D_Ko_C22_BL_CHK_M.CpG_5mC.cov.gz','rt'); tab <- read.table(ztab,header=F);tab <- subset(tab, (V5 + V6) > 5);D_Ko_C22<- tab[,c(1,2,5,6,4)]; names(D_Ko_C22) <- c('chr','start','Cs_D_Ko_C22','Ts_D_Ko_C22','D_Ko_C22');
ztab <- gzfile('D_Ko_C40_BL_CHK_M.CpG_5mC.cov.gz','rt'); tab <- read.table(ztab,header=F);tab <- subset(tab, (V5 + V6) > 5);D_Ko_C40<- tab[,c(1,2,5,6,4)]; names(D_Ko_C40) <- c('chr','start','Cs_D_Ko_C40','Ts_D_Ko_C40','D_Ko_C40');
ztab <- gzfile('D_Ko_C42_BL_CHK_M.CpG_5mC.cov.gz','rt'); tab <- read.table(ztab,header=F);tab <- subset(tab, (V5 + V6) > 5);D_Ko_C42<- tab[,c(1,2,5,6,4)]; names(D_Ko_C42) <- c('chr','start','Cs_D_Ko_C42','Ts_D_Ko_C42','D_Ko_C42');
ztab <- gzfile('D_Ba_H02_BL_CHK_M.CpG_5mC.cov.gz','rt'); tab <- read.table(ztab,header=F);tab <- subset(tab, (V5 + V6) > 5);H_Ba_H02<- tab[,c(1,2,5,6,4)]; names(H_Ba_H02) <- c('chr','start','Cs_H_Ba_H02','Ts_H_Ba_H02','H_Ba_H02');
ztab <- gzfile('D_Ba_H09_BL_CHK_M.CpG_5mC.cov.gz','rt'); tab <- read.table(ztab,header=F);tab <- subset(tab, (V5 + V6) > 5);H_Ba_H09<- tab[,c(1,2,5,6,4)]; names(H_Ba_H09) <- c('chr','start','Cs_H_Ba_H09','Ts_H_Ba_H09','H_Ba_H09');
ztab <- gzfile('D_Ba_H19_BL_CHK_M.CpG_5mC.cov.gz','rt'); tab <- read.table(ztab,header=F);tab <- subset(tab, (V5 + V6) > 5);H_Ba_H19<- tab[,c(1,2,5,6,4)]; names(H_Ba_H19) <- c('chr','start','Cs_H_Ba_H19','Ts_H_Ba_H19','H_Ba_H19');
ztab <- gzfile('D_Ba_H21_BL_CHK_M.CpG_5mC.cov.gz','rt'); tab <- read.table(ztab,header=F);tab <- subset(tab, (V5 + V6) > 5);H_Ba_H21<- tab[,c(1,2,5,6,4)]; names(H_Ba_H21) <- c('chr','start','Cs_H_Ba_H21','Ts_H_Ba_H21','H_Ba_H21');
ztab <- gzfile('D_Hi_C03_BL_CHK_M.CpG_5mC.cov.gz','rt'); tab <- read.table(ztab,header=F);tab <- subset(tab, (V5 + V6) > 5);H_Hi_C03<- tab[,c(1,2,5,6,4)]; names(H_Hi_C03) <- c('chr','start','Cs_H_Hi_C03','Ts_H_Hi_C03','H_Hi_C03');
ztab <- gzfile('D_Ku_H02_BL_CHK_M.CpG_5mC.cov.gz','rt'); tab <- read.table(ztab,header=F);tab <- subset(tab, (V5 + V6) > 5);H_Ku_H02<- tab[,c(1,2,5,6,4)]; names(H_Ku_H02) <- c('chr','start','Cs_H_Ku_H02','Ts_H_Ku_H02','H_Ku_H02');
ztab <- gzfile('D_Lo_C19_BL_CHK_M.CpG_5mC.cov.gz','rt'); tab <- read.table(ztab,header=F);tab <- subset(tab, (V5 + V6) > 5);H_Lo_C19<- tab[,c(1,2,5,6,4)]; names(H_Lo_C19) <- c('chr','start','Cs_H_Lo_C19','Ts_H_Lo_C19','H_Lo_C19');
ztab <- gzfile('D_Lo_C20_BL_CHK_M.CpG_5mC.cov.gz','rt'); tab <- read.table(ztab,header=F);tab <- subset(tab, (V5 + V6) > 5);H_Lo_C20<- tab[,c(1,2,5,6,4)]; names(H_Lo_C20) <- c('chr','start','Cs_H_Lo_C20','Ts_H_Lo_C20','H_Lo_C20');
ztab <- gzfile('D_Ne_Y03_BL_CHK_M.CpG_5mC.cov.gz','rt'); tab <- read.table(ztab,header=F);tab <- subset(tab, (V5 + V6) > 5);H_Ne_Y03<- tab[,c(1,2,5,6,4)]; names(H_Ne_Y03) <- c('chr','start','Cs_H_Ne_Y03','Ts_H_Ne_Y03','H_Ne_Y03');
ztab <- gzfile('D_Ne_Y14_BL_CHK_M.CpG_5mC.cov.gz','rt'); tab <- read.table(ztab,header=F);tab <- subset(tab, (V5 + V6) > 5);H_Ne_Y14<- tab[,c(1,2,5,6,4)]; names(H_Ne_Y14) <- c('chr','start','Cs_H_Ne_Y14','Ts_H_Ne_Y14','H_Ne_Y14');
ztab <- gzfile('D_Ne_Y36_BL_CHK_M.CpG_5mC.cov.gz','rt'); tab <- read.table(ztab,header=F);tab <- subset(tab, (V5 + V6) > 5);H_Ne_Y36<- tab[,c(1,2,5,6,4)]; names(H_Ne_Y36) <- c('chr','start','Cs_H_Ne_Y36','Ts_H_Ne_Y36','H_Ne_Y36');
ztab <- gzfile('D_Ne_Y42_BL_CHK_M.CpG_5mC.cov.gz','rt'); tab <- read.table(ztab,header=F);tab <- subset(tab, (V5 + V6) > 5);H_Ne_Y42<- tab[,c(1,2,5,6,4)]; names(H_Ne_Y42) <- c('chr','start','Cs_H_Ne_Y42','Ts_H_Ne_Y42','H_Ne_Y42');
ztab <- gzfile('D_Rb_Y07_BL_CHK_M.CpG_5mC.cov.gz','rt'); tab <- read.table(ztab,header=F);tab <- subset(tab, (V5 + V6) > 5);H_Rb_Y07<- tab[,c(1,2,5,6,4)]; names(H_Rb_Y07) <- c('chr','start','Cs_H_Rb_Y07','Ts_H_Rb_Y07','H_Rb_Y07');
ztab <- gzfile('D_Rb_Y08_BL_CHK_M.CpG_5mC.cov.gz','rt'); tab <- read.table(ztab,header=F);tab <- subset(tab, (V5 + V6) > 5);H_Rb_Y08<- tab[,c(1,2,5,6,4)]; names(H_Rb_Y08) <- c('chr','start','Cs_H_Rb_Y08','Ts_H_Rb_Y08','H_Rb_Y08');
ztab <- gzfile('D_Rb_Y15_BL_CHK_M.CpG_5mC.cov.gz','rt'); tab <- read.table(ztab,header=F);tab <- subset(tab, (V5 + V6) > 5);H_Rb_Y15<- tab[,c(1,2,5,6,4)]; names(H_Rb_Y15) <- c('chr','start','Cs_H_Rb_Y15','Ts_H_Rb_Y15','H_Rb_Y15');
ztab <- gzfile('S_Up_H77_BL_CHK_M.CpG_5mC.cov.gz','rt'); tab <- read.table(ztab,header=F);tab <- subset(tab, (V5 + V6) > 5);S_Up_H77<- tab[,c(1,2,5,6,4)]; names(S_Up_H77) <- c('chr','start','Cs_S_Up_H77','Ts_S_Up_H77','S_Up_H77');
ztab <- gzfile('S_Up_H80_BL_CHK_M.CpG_5mC.cov.gz','rt'); tab <- read.table(ztab,header=F);tab <- subset(tab, (V5 + V6) > 5);S_Up_H80<- tab[,c(1,2,5,6,4)]; names(S_Up_H80) <- c('chr','start','Cs_S_Up_H80','Ts_S_Up_H80','S_Up_H80');
ztab <- gzfile('S_Up_H81_BL_CHK_M.CpG_5mC.cov.gz','rt'); tab <- read.table(ztab,header=F);tab <- subset(tab, (V5 + V6) > 5);S_Up_H81<- tab[,c(1,2,5,6,4)]; names(S_Up_H81) <- c('chr','start','Cs_S_Up_H81','Ts_S_Up_H81','S_Up_H81');
#and the 4 CG samples
ztab <- gzfile('D_Ko_C31_BL_CHK_M.CpG_5mC.cov.gz','rt'); tab <- read.table(ztab,header=F);tab <- subset(tab, (V5 + V6) > 5);D_Ko_C31<- tab[,c(1,2,5,6,4)]; names(D_Ko_C31) <- c('chr','start','Cs_D_Ko_C31','Ts_D_Ko_C31','D_Ko_C31');
ztab <- gzfile('D_Ko_C36_BL_CHK_M.CpG_5mC.cov.gz','rt'); tab <- read.table(ztab,header=F);tab <- subset(tab, (V5 + V6) > 5);D_Ko_C36<- tab[,c(1,2,5,6,4)]; names(D_Ko_C36) <- c('chr','start','Cs_D_Ko_C36','Ts_D_Ko_C36','D_Ko_C36');
ztab <- gzfile('S_Up_H75_BL_CHK_M.CpG_5mC.cov.gz','rt'); tab <- read.table(ztab,header=F);tab <- subset(tab, (V5 + V6) > 5);S_Up_H75<- tab[,c(1,2,5,6,4)]; names(S_Up_H75) <- c('chr','start','Cs_S_Up_H75','Ts_S_Up_H75','S_Up_H75');
ztab <- gzfile('S_Up_H59_BL_CHK_M.CpG_5mC.cov.gz','rt'); tab <- read.table(ztab,header=F);tab <- subset(tab, (V5 + V6) > 5);S_Up_H59<- tab[,c(1,2,5,6,4)]; names(S_Up_H59) <- c('chr','start','Cs_S_Up_H59','Ts_S_Up_H59','S_Up_H59');
# merge the points, keeping all
masterHZ <- Reduce(function(x,y) merge(x = x, y = y, by=c('chr','start'),all=TRUE),
                   list(H_Ba_H02,H_Ba_H09,H_Ba_H19,H_Ba_H21,H_Hi_C03,D_Ko_C22,D_Ko_C40,D_Ko_C42,H_Ku_H02,H_Lo_C19,H_Lo_C20,H_Ne_Y03,H_Ne_Y14,H_Ne_Y36,H_Ne_Y42,H_Rb_Y07,H_Rb_Y08,H_Rb_Y15,S_Up_H77,S_Up_H80,S_Up_H81,D_Ko_C31,D_Ko_C36,S_Up_H75,S_Up_H59))
head(masterHZ)

#count missing by individual
na_count <-sapply(masterHZ[!grepl('^Ts|Cs',names(masterHZ))], function(y) sum(length(which(is.na(y)))))

#filter according to missingness
covdat <- masterHZ[rowSums(is.na(masterHZ[grepl('^Cs_S|^Cs_D', names(masterHZ))])) == 0, ]
covdat2 <- covdat[rowSums(is.na(covdat[grepl('^Cs_H', names(covdat))])) <= 3, ]
nrow(covdat)
nrow(covdat2)
covdat <- covdat2

#and filter for maximum coverage:
covs <- covdat[grepl('^chr|^start|^Cs_|^Ts_',names(covdat))]
covs$total <- rowSums(covs[grepl('Ts_|Cs_',names(covs))],na.rm=TRUE)
thresh <- quantile(covs$total,0.99)
keep <- subset(covs,total < thresh)
keep$chr_pos <- paste0(keep$chr,"_",keep$start)
keepCOV <- keep[,c('chr_pos','chr')]

#count positions by chromosome
covdat %>% dplyr::count(chr)
masterHZ <- covdat
masterHZ <- masterHZ[!grepl('^Cs_|^Ts_',names(masterHZ))]

HZdats <- masterHZ %>% arrange(chr,start)

#add site and chrZ marker
HZdats$chr_pos <- paste0(HZdats$chr,"_",HZdats$start)
HZdats <- HZdats %>% mutate(AvZ = ifelse(chr == 'chrZ','chrZ','Autosome'))

#remove scaffolds
HZdatsC <- HZdats[!grepl('scaf',HZdats$chr),]

#reorganize to bed
HZdatsCB <- HZdatsC[,c(1,2,2,3:ncol(HZdatsC))]
names(HZdatsCB)[names(HZdatsCB) == 'start.1'] <- 'end'

#and finally, only keep the sites that are within the 99% coverage limits
HZfin <- merge(keepCOV,HZdatsCB,by=c('chr_pos','chr'),all.x=TRUE)
nrow(HZfin)
#reorganize to bed
HZdatsCBdivBED <- HZfin[,c(2,3,4,1,5:ncol(HZfin))]

#also output PCA input
masterHZ$chr_pos <- paste0(masterHZ$chr,"_",masterHZ$start)
pcainput <- merge(keepCOV,masterHZ,by=c('chr','chr_pos'),all.x=TRUE)
write.table(pcainput,'Variance_PCA_Input.txt',quote=F,row.names=FALSE)

#perform PCA
pc1 <- pcainput[!grepl('^Cs_|^Ts|^chr|^start',names(pcainput))]
pc2 <- pc1[rowSums(is.na(pc1[grepl('_', names(pc1))])) == 0, ]
nrow(pc1)
nrow(pc2)
#remove sites with no variation
pc2$var <- apply(pc2, 1, var, na.rm=TRUE)
pc3 <- subset(pc2,var > 0)
pc4 <- t(pc3[!grepl('var',names(pc3))])
pca <- prcomp(pc4[,c(1:ncol(pc4))],scale=TRUE)
var_explained <- pca$sdev^2/sum(pca$sdev^2)
pcd <- data.frame(pca$x)
pcd$ID <- row.names(pcd)
pcdm <- as.data.frame(pcd) %>% separate(ID, into = c('SPECIES','ID1','ID2'))
pcdm$ID <- row.names(pcd)
pcdm <- pcdm %>% mutate(Experiment = ifelse(ID == 'D_Ko_C31' | ID == 'D_Ko_C36' | ID == 'S_Up_H75' | ID == 'S_Up_H59','CG','HZ'))
l <- ggplot(pcdm,aes(x=PC1,y=PC2,col=SPECIES,size=Experiment,shape=SPECIES))+
  geom_point(show.legend=TRUE)+
  scale_size_manual(values=c(10,4))+
  scale_color_manual(values=c('grey10','darkorchid2','grey70'))+
  xlab(paste0('PC1 : ',round(var_explained[1],3)*100,'% Explained'))+
  ylab(paste0('PC2 : ',round(var_explained[2],3)*100,'% Explained'))+
  theme_classic()
l

pdf('Batch-Effects.pdf',width=5,height=4)
l
dev.off()
```



## Technical Replicates

RRBS vs WGBS, interindividual vs tech rep correlation: 

```R
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/crow_scratch/2021-07_HybridZone/Nov/tech_reps/2023_09')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(matrixStats)
library(viridis)
library(ggpubr)

# list of all .cov.gz files in the working directory
files <- list.files(path=".", pattern="*.cov.gz", full.names=TRUE, recursive=FALSE)

# Initialize empty  dat frames to store names and coverage data
master = NULL
cov = NULL

# Loop through all .cov.gz files
for (samp in files) {
  cat('Reading in file: ',samp,'\n') # Display current file being read
  
  # Read in .cov.gz file
  ztab <- gzfile(samp,'rt');
  tab <- read.table(ztab,header=F);
  names(tab) = c('chr','start','end','Percent_5mC','countM','countU')
  
  # Extract name for new variable, trim all the garbage, modify as necessary depending on your file extensions....
  name <- gsub('./','',samp); name <- gsub('_val_1.*','',name); name <- gsub('_trimmed.*','',name); name <- gsub('.cov.gz','',name)
  
  # Calculate quick coverage and %5mc summary statistics per sample and store them
  tab = tab %>% mutate(Coverage=countM+countU,ID=name)

  # Keep only positions with 10x coverage or higher ### IMPORTANT!
  tab <- subset(tab, Coverage > 9)
  tab$site <- paste0(tab$chr,'_',tab$start)
  tab$Percent_5mC <- tab$Percent_5mC/100
  
  # Grab only site / countM / countU / %5mC
  tab2 <- tab %>% select(site,Percent_5mC)
  names(tab2) <- c('site',name);
  
  cat('Saved to variable: ',name,'\n') # Indicate current file has been processed, then bind with all the data
  if(is.null(master)){
    master <- tab2
  }else{
    master <- full_join(master, tab2, by = "site")
  }
}

class = read_tsv('classifications.bed',col_names = F)
names(class) = c('chr','start','class','region')
class = class %>% mutate(site = paste0(chr,'_',start))
mc = left_join(master,class)
#samples for analysis
samps = c('C31','C45','H59','H60')
dat=NULL
for (samp in samps) {
  cat('workin on sample: ',samp,'\n')
  #random non-tech sample
  avail = samps[!grepl(samp,samps)]
  othersamp = sample(avail, 1)
  
  m = mc %>% select(contains(c('site','chr','class','region',samp,othersamp))) %>% as_tibble %>% drop_na
  names(m) = c('site','chr','class','region','s1','s2','d3','d4')
  c1 = m %>% summarize(cor = cor(s1,s2,method='spearman'),
                       bcor1 = cor(s1,d4,method='spearman'),
                       bcor2 = cor(s2,d3,method='spearman'),
                       nocor = mean(bcor1,bcor2),
                       sample = samp,
                       Sites = paste0(round(nrow(.)/1000,0),'K'))
  dat = rbind(dat,c1)
}
dat = dat %>% select(cor,nocor,sample,Sites)
names(dat) = c('RRBS_WGBS','Inter_Individual','Sample','Sites')
rp = dat %>% pivot_longer(!c(Sample,Sites),names_to = 'Comparison') %>% 
  ggplot(aes(x=Sample,y=value,col=Comparison,shape=Comparison,label=Sites))+
  geom_text(aes(y=0.915),col='black')+
  geom_point(size=3)+xlab('Sample')+ylab("Spearman's Rho")+
  theme_bw()
rp

png('TechRep_2023SEPT28.png',units='in',res=600,height=5,width=7)
rp
dev.off()


```

