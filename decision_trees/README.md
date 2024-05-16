# Decision Trees Assessing 5mC ~ Chromosomal Variation

This uses as input largely the `RF-Input-23JAN31.txt.gz` found in `/plotting_files/`. 

Figures at bottom. 

Take the classified data, and prepare it for RF/XGB:

Great tutorial: https://www.kirenz.com/post/2021-02-17-r-classification-tidymodels/

Jan 31 2023:

```R
.libPaths('~/miniconda3/envs/tidymodels/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/crow_scratch/2021-07_HybridZone/Nov/Classification_Sept/RF_JAN23/')
library(dplyr)
library(ggplot2)
library(ggsci)
library(reshape2)
library(ggpubr)
library(viridis)
library(GGally)
library(matrixStats)
.libPaths('~/miniconda3/envs/r/lib/R/library')
library(tidyverse)

options(scipen=999)

all <- read.table('../Full.DSS.BP-10x.Classified-Annotated_23JAN04.txt',header=T)
all <- all %>% arrange(chr,start)
all <- all %>% mutate_at(grep("HapD|FuLi|Tajima|Dxy|FST|GC",colnames(.)),funs(as.numeric))
all <- subset(all,Region != 'Missing' & Region != '.')
c = read.table('../CG_DSS-10x.stat',header=TRUE)
w = read.table('../WGBS_DSS-5x.stat',header=TRUE)
h = read.table('../HZ_DSS-10x.stat',header=TRUE)

#merge data frames
f = left_join(all,c) %>% 
  left_join(.,w) %>% 
  left_join(.,h) %>% 
  dplyr::rename(mstat.cg = CG_stat,
                mstat.wgbs = WGBS_stat,
                mstat.hz = HZ_stat)  

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

for (xp in xps){
  cat('Creating output for: ',xp,'\n')
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
    ggplot(aes(x=value,fill=name))+
    geom_density()+scale_fill_manual(values=viridis(4,option='inferno'))+
    facet_wrap(name~variable,scales='free',ncol=4)+
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
    ylab('')+theme_classic()+ggtitle('Binary')
  ##### TRINARY
  tdat = dat
  thresh = tdat %>% group_by(variable) %>% summarize(ThresholdHi = quantile(valueuse,0.8),
                                                     ThresholdLo = quantile(valueuse,0.2))
  tdat = left_join(tdat,thresh)
  tdat = tdat %>% mutate(Classified = ifelse(valueuse >=ThresholdHi, 'High',
                                             ifelse(valueuse >=ThresholdLo,'Intermediate',
                                                    ifelse(valueuse < ThresholdLo,'Low','Miss'))))
  tdat %>% dplyr::count(Classified)
  trip = tdat %>% ggplot(aes(x=valueuse,fill=variable)) + 
    facet_wrap(~variable,scales='free')+
    geom_histogram()+
    scale_fill_viridis(discrete=TRUE,option='mako')+
    geom_vline(data=thresh, aes(xintercept=ThresholdHi),lty=2,lwd=0.5,col='yellow')+
    geom_vline(data=thresh, aes(xintercept=ThresholdLo),lty=2,lwd=0.5,col='yellow')+
    ylab('')+theme_classic()+ggtitle('Trinary')
  gp = ggarrange(binp,trip,common.legend = TRUE,nrow = 2)
  assign(paste0(xp,'_thresh'),gp)
  
  #merge them 
  td = tdat %>% select(-contains(c('Thresh','valueuse'))) %>% dplyr::rename(Trinary=Classified)
  bd = bdat %>% select(-contains('Thresh')) %>% dplyr::rename(Binary=Classified)
  ad = full_join(td,bd)
  ad = ad %>% dplyr::rename(Continuous = valueuse)
  ad = ad %>% mutate(Relativepos = GenStart / Length)
  ad$Subset = toupper(xp)
  rfdat = rbind(rfdat,ad)
}

pdf('Raw_Distributions.pdf',height=5,width=5)
ggarrange(wgbs_values,cg_values,hz_values,common.legend = TRUE,nrow=3)
dev.off()

pdf('Raw_Thresholds.pdf',height=10,width=4)
ggarrange(wgbs_thresh,cg_thresh,hz_thresh,common.legend = TRUE,nrow=3)
dev.off()

#save frame 
write.table(rfdat,file='RF-Input-23JAN31.txt',quote=F,sep='\t',row.names=F)

#output correlations among variables
cord = NULL
xps = c('HZ','CG','WGBS')
for (xp in xps){
  d1 = rfdat %>% filter(Subset == xp & variable == 'Mean')
  np = d1 %>% 
    select(Length,Relativepos,Region,HapD,FuLiD,TajimaD,Dxy,FST,GC) %>% 
    correlate(method='spearman') %>%    # Create correlation data frame (cor_df)
    network_plot(min_cor = 0,legend='range')
  assign(paste0(xp,'_cor'),np)
  
  x = d1 %>% 
    select(Length,Relativepos,Region,HapD,FuLiD,TajimaD,Dxy,FST,GC) %>% 
    correlate() %>%    # Create correlation data frame (cor_df)
    rearrange() %>%  # rearrange by correlations
    shave() # Shave off the upper triangle for a clean result
  x$experiment = xp
  cord = rbind(cord,x)
}

pdf('Correlations_RF_Input.pdf',height=6,width=12)
ggarrange(WGBS_cor,CG_cor,HZ_cor,labels=c('ComGar.WGBS','ComGar.RRBS','HybZon.RRBS'),nrow=1)
dev.off()

write.table(cord,file='Correlations_RF_Input.txt',quote=F,sep='\t',row.names=F)
```

Once the data is prepared, grow the forests: 

```R
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
.libPaths('~/miniconda3/envs/tidymodels/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/crow_scratch/2021-07_HybridZone/Nov/Classification_Sept/RF_JAN23')
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

dats <- read.table('RF-Input-23JAN13.txt',header=T)
type = args[1]  #e.g. Binary,Trinary,Regression
vari = args[2] #e.g. Mean,Max
xp = args[3] #subset experiment, CG WGBS HZ
iter = args[4] #iteration

dat = dats %>% filter(variable == vari & Subset == xp)
if (type == 'Regression') {
  cat('Response type is Continuous, taking continuous values\n')
  rf = dat %>% select(c(Continuous,Length,Relativepos,Region,HapD,FuLiD,TajimaD,Dxy,FST,GC))
  rf = rf %>% dplyr::rename(Response = Continuous)
} else {
  cat('Response type is Binary/Trinary, taking classified values\n')
  rf = dat %>% select(c(as.symbol(type),Length,Relativepos,Region,HapD,FuLiD,TajimaD,Dxy,FST,GC))
  rf = rf %>% dplyr::rename(Response = as.symbol(type))
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
  cat('Response type is Binary/Trinary, using classification engines\n')
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
```

And plot results:

```R
.libPaths('~/mambaforge/envs/r/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/crow_scratch/2021-07_HybridZone/Nov/Classification_Sept/RF_JAN23/results')
#setwd('F:/Research/scratch/crow_hybrid_paper/rf/2023_01/')
#### Randomforest Classification
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
fpa = ggarrange(fpr,fpc,widths = c(0.6,1.2))

#save
pdf('RF_Fit.pdf',height=5,width=8)
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
            hi = ci(PctImp)[2],
            lo = ci(PctImp)[3],
            sd = ci(PctImp)[4])
write.table(viS,file='Permutation_Importance.txt',quote=F,sep='\t',row.names=F)

#save
pdf('RF_VIP_Regression.pdf',height=6,width=6)
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
ALL$Type = factor(ALL$Type,levels=c('Regression','Binary','Trinary'))
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
pdf('RF_VIP_All.pdf',height=6,width=10)
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
  facet_wrap(Subset~Response,scales='free',nrow=2)
cp

#save
png('RF_Predictions.png',height=7,width=11,res=600,units='in')
pdf('RF_Predictions.pdf',height=5,width=2.5)
cp
dev.off()

#Spearman's rank correlation 
conf %>% 
  group_by(Response,Type,Subset) %>% 
  summarize(Spearman = cor(Predicted,Truth,method='spearman'))

# # A tibble: 6 × 4
# # Groups:   Response, Type [2]
# Response Type       Subset Spearman
# <chr>    <chr>      <chr>     <dbl>
#   1 Max      Regression CG        0.309
# 2 Max      Regression HZ        0.288
# 3 Max      Regression WGBS      0.529
# 4 Mean     Regression CG        0.477
# 5 Mean     Regression HZ        0.341
# 6 Mean     Regression WGBS      0.193
```

## Rho: 5mC & Population Genetic

For genetic v 5mC correlations using the RF input (within 5kb windows):

```R
.libPaths('~/miniconda3/envs/r/lib/R/library')
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/crow_scratch/2021-07_HybridZone/Nov/Classification_Sept/RF_JAN23/init4/')
#setwd('F:/Research/scratch/crow_hybrid_paper/rf/2023_01/')
#### Randomforest Classification
library(tidyverse)
library(viridis)
library(forcats)
library(ggpubr)
library(rstatix)
library(gmodels)
library(RColorBrewer)

#plot some of the interesting comparisons
rf = read.table('../RF-Input-23JAN31.txt',header=TRUE)

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

write.table(fullboots,file=paste0('Genetic5mC_Bootstraped_Results-2023FEB14.txt'),quote=F,sep='\t',row.names=F)
fullboots = read.table('Genetic5mC_Bootstraped_Results-2023FEB14.txt',header=TRUE)

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
  geom_violin(alpha=0.75,draw_quantiles = c(0.025,0.975))+
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


![Correlations of input variables](plotting_files/Correlations_RF_Input.pdf)

* Correlations of input variables 


![Raw Distributions](plotting_files/Raw_Distributions.pdf)

* Raw distributions of input 5mC 


![Raw Distributions Thresholds](plotting_files/Raw_Thresholds.pdf)

* Thresholds used for binary / trinary classification of distributions of input 5mC 


![Regression Predictions](plotting_files/RF_Predictions.png)

* Predictions from randomforest regression on testing set. 
