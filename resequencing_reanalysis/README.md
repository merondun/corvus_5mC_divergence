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

