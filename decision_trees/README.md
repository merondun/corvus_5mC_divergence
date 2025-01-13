# Decision Trees 

Take the classified data, and prepare it for RF:

Great tutorial: https://www.kirenz.com/post/2021-02-17-r-classification-tidymodels/

## Prep

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

## Run RF & XGB

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

## Plot

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


