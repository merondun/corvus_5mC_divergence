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

