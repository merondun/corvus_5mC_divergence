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


