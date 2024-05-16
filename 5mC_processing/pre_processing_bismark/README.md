## Obtain bismark cov.gz files

### ComGar.RRBS

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

### HybZon.RRBS

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

### ComGar.WGBS

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

## QC

### Bisulfite Conversion

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

### M-Bias & Counts

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
```


