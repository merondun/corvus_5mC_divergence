## Resequencing Ordinations

Resequencing PCA. I run this once for the background, and then again for the focal region: 

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

#read metadat
variables <- read.table('../k28.metdat',header=F)
names(variables) <- c('ID','Locality','Sex','Taxon')

#Plot 
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
```

![PCA background](plotting_files/Reseq_PCA-Background.pdf)

* PCA on autosomal background

![PCA focal](plotting_files/Reseq_PCA-Focal.pdf)

* PCA focal region
