rm(list=ls())
library(brms)
library(stringr)
library(dplyr)



# colors ---------------------------------------------------------

cola<-c(
  rgb(1,0,0,alpha=1),
  rgb(0.2,0.2,1,alpha=1),
  rgb(0, 0.7, 0,alpha=1),
  rgb(0.7,0.7,0.1,alpha=1))

alp=0.18;colb<-c(
  rgb(1,0,0,alpha=alp),
  rgb(0.2,0.2,1,alpha=alp),
  rgb(0, 0.9, 0,alpha=alp),
  rgb(0.7,0.7,0.1,alpha=alp))

colc<-c(
  rgb(1,0,0,alpha=0.8),
  rgb(0.2,0.2,1,alpha=0.8),
  rgb(0, 0.9, 0,alpha=0.8),
  rgb(0.7,0.7,0.1,alpha=0.8))

# Defining priors -------------------------
prior_group <- c(
set_prior("normal(0,30)", class = "b", coef = "groupFinland.Males"),
set_prior("normal(0,30)", class = "b", coef = "groupNorway.Males"),
set_prior("normal(0,30)", class = "b", coef = "groupSweden.Males"),
set_prior("normal(0,30)", class = "b", coef = "groupDenmark.Males"),
set_prior("normal(0,30)", class = "b", coef = "groupNorway.Females"),
set_prior("normal(0,30)", class = "b", coef = "groupFinland.Females"),
set_prior("normal(0,30)", class = "b", coef = "groupSweden.Females"))


# Data, models and posterior sample extraction ---------------------------------------------------------

## melanoma data -------------
urlfile <- "https://raw.githubusercontent.com/filip-tichanek/nord_melanoma/main/source_code/melanoma_1y.csv"
melanoma_1y <- read.csv(url(urlfile), sep = ';')

urlfile <- "https://raw.githubusercontent.com/filip-tichanek/nord_melanoma/main/source_code/melanoma_5y.csv"
melanoma_5y <- read.csv(url(urlfile), sep = ';')


## Data wrangling
x<-2;melanoma_1y_est<-data.frame(melanoma_1y[,1])
repeat{melanoma_1y_est[,x]<-str_sub(melanoma_1y[,x],1,4);x<-x+1
if(x>9){break}}
repeat{melanoma_1y_est[,x]<-str_sub(melanoma_1y[,x-8],6,9);x<-x+1;if(x>17){break}}
x<-2;melanoma_5y_est<-data.frame(melanoma_5y[,1])
repeat{melanoma_5y_est[,x]<-str_sub(melanoma_5y[,x],1,4);x<-x+1
if(x>9){break}}
repeat{melanoma_5y_est[,x]<-str_sub(melanoma_5y[,x-8],6,9);x<-x+1;if(x>17){break}}

melanoma<-(data.frame(unlist(melanoma_1y_est[,2:9])));colnames(melanoma)<-"surv_1y"
melanoma$cil_1y<-as.numeric(unlist(melanoma_1y_est[,10:17]))
melanoma$surv_5y<-as.numeric(unlist(melanoma_5y_est[,2:9]))
melanoma$cil_5y<-as.numeric(unlist(melanoma_5y_est[,10:17]))
melanoma$year<-rep(seq(1973,2018,by=5),8)
melanoma$sex<-c(rep("Males",40),rep("Females",40))
melanoma$country<-c(rep(c(rep("Denmark",10),rep("Finland",10),rep("Norway",10),rep("Sweden",10)),2))
melanoma$shou<-c(rep(c(rep("den_mal_",10),rep("fin_mal_",10),rep("nor_mal_",10),rep("swe_mal_",10)
                     ,rep("den_fem_",10),rep("fin_fem_",10),rep("nor_fem_",10),rep("swe_fem_",10)),1))
melanoma$group<-interaction(melanoma$country,melanoma$sex)
melanoma$years10cen<-(melanoma$year-1995.5)/10
melanoma$surv_1y<-as.numeric(melanoma$surv_1y)
melanoma$surv_5y<-as.numeric(melanoma$surv_5y)
melanoma$cil_1y<-as.numeric(melanoma$cil_1y)
melanoma$cil_5y<-as.numeric(melanoma$cil_5y)
melanoma$se_1y<-(melanoma$surv_1y-melanoma$cil_1y)/1.96
melanoma$se_5y<-(melanoma$surv_5y-melanoma$cil_5y)/1.96
melanoma$surv_cond<-(melanoma$surv_5y/melanoma$surv_1y)*100

## melanoma model ------------------------------------------------------------
set.seed(17)
melanoma_1y_model<-brm(surv_1y|se(se_1y)~group+s(years10cen,by=group,k=5)
                    ,family="Gaussian",save_pars = save_pars(all = TRUE),prior=prior_group,seed=17
                    ,data=melanoma,iter=7000, warmup=2000,chains=2,cores=1,control = list(adapt_delta = 0.98))
summary(melanoma_1y_model)
plot(melanoma_1y_model)
fixef(melanoma_1y_model)
conditional_effects(melanoma_1y_model,effects='years10cen:group')

melanoma_5y_model<-brm(surv_5y|se(se_5y)~group+s(years10cen,by=group,k=5)
                    ,family="Gaussian",save_pars = save_pars(all = TRUE),prior=prior_group,seed=17
                    ,data=melanoma,iter=7000, warmup=2000,chains=2,cores=1,control = list(adapt_delta = 0.98))
conditional_effects(melanoma_5y_model,effects='years10cen:group')
fixef(melanoma_5y_model)
conditional_effects(melanoma_5y_model,effects='years10cen:group')

## Model diagnostics ---------------------

pp_check(melanoma_1y_model, type='dens_overlay',ndraws = 50)
pp_check(melanoma_1y_model, type='scatter_avg',ndraws = 50) 
pp_check(melanoma_1y_model, type="stat_2d", stat = c("max", "min"),ndraws=200)
pp_check(melanoma_1y_model, type="stat_2d", stat = c("mean", "sd"),ndraws=200)


pp_check(melanoma_5y_model, type='dens_overlay',ndraws = 50)
pp_check(melanoma_5y_model, type='scatter_avg',ndraws = 50) 
pp_check(melanoma_5y_model, type="stat_2d", stat = c("max", "min"),ndraws=200)
pp_check(melanoma_5y_model, type="stat_2d", stat = c("mean", "sd"),ndraws=200)
## >> Everything seems fine: data generated from estimated parameters have similar 
## >> distribution as are real data


# Functions ---------------
# 
# All the functions below serve for plotting of nonlinear trends in survival and related uncertainty. 
# All the functions below use extracted posterior samples as an input *data*. 
# The samples represent estimation of survival over the 50 years (1971 to 2020) for specific country, 
# cancer, and type of survival (1-year, 5-years, conditional [5y/1y]. 
# The time period of 50 years must be divided to sequence of 500 numbers (1/10 of a year).


## brekapo function ---------

## Function serves to identify *breaking points*, i.e. times when the annual change of survival 
## changed with at least 95% plausibility. This was assessed by calculation of the 
## 2nd derivation of the given survival measure and its 95% CI; the ‘breaking point’ 
## was defined as a peak value within at least a 3-year interval where 95% CI for 
## the 2nd derivation did not cross zero. If the 2nd derivation is plausibly non-zero 
## for at least 3 years, the function takes the peak in the 2nd derivation 
## (within the identified time interval) as the *breaking point* 
## (must be between the years 1976 and 2016)
## 
## There is one additional argument 'arb' which should be zero except for the situation when 
## the multiple breaking points overlap. The argument only move the breaking points 
## by given value to avoid overlapping. 

## The function often gives warning: "Error in xy.coords(x, y): 'x' and 'y' lengths differ". 
## Warning occurs when there is no 'breaking point' in a given group and survival measure, 
## ie, it is not real error. 

breakpo<-function(data,arb){ 
  data<-data.frame((data[,-1] - data[,-ncol(data)])*10)
  data<-data.frame(data[,-1] - data[,-ncol(data)])
  data=sapply(data, function(p) quantile(p, probs = c(0.025,0.975,0.5)))
  cbinl<-c()
  x=1
  repeat{
    cbinl[x]<-
      if(data[1,x]>0|data[2,x]<0){cbinl[x]=1}else{cbinl[x]=0}
    x<-x+1
    if(x>length(data[1,])){break}}
  cbin=1;x<-1
  repeat{
    cbin[x+1]<-abs(cbinl[x]-cbinl[x+1])
    x=x+1
    if(x>(length(cbinl)-1)){break}}
  data<-data.frame(rbind(data,cbin));data[5,]<-yreal[c(2:499)]-0.049
  data[6,]<-1:498;data[4,51]<-1;data[4,449]<-1
  
  row.names(data)<-c("cil","ciu","est","stat_change","year","timepoint")
  data2=t(data[,data[4,]==1])
  data2<-data.frame(data2)
  tr<-subset(data2,data2$timepoint>50&data2$timepoint<450)
  y=1
  bp<-c()
  bx<-1
  repeat{
    if( (tr[y,1]<0) & (tr[y,2]<0) & (  (tr$year[y+1]-tr$year[y])>3   )  ) {
      tr2<-data[,tr[y,6]:(tr[y+1,6]-2)]
      bp[bx]<- tr2[,order(tr2[3,],decreasing=F)[1]][5]
      bx=bx+1
    }
    if( (tr[y,1]>0) & (tr[y,2]>0) & (  (tr$year[y+1]-tr$year[y])>3   )) {
      tr2<-data[,tr[y,6]:(tr[y+1,6]-2)]
      bp[bx]<- tr2[,order(tr2[3,],decreasing=T)[1]][5]
      bx=bx+1
    }
    y=y+1;if(y>(dim(tr)[1]-1)){break}}
  y<-1;repeat{
    lines(c(bp[y]+arb,bp[y]+arb),c(range[1],range[1]+0.025*scal),col=cola[xx],lwd=3.5,lend=1)
    lines(c(bp[y]+arb,bp[y]+arb),c(range[2],range[2]-0.025*scal),col=cola[xx],lwd=3.5,lend=1)
    lines(c(bp[y]+arb,bp[y]+arb),c(range[1],range[1]+0.999*scal),col=colc[xx],lwd=1,lend=1,lty=2)
    y=y+1;if(y>length(bp)){break}}
  print(bp)
  print(tr)}


## 'polyg_surv' function -------------

## for drawing 95% credible interval for survival indicators

polyg_surv<-function(data){ 
  data<-data.frame(data)
  data<-sapply(data, function(p) quantile(p, probs = c(0.025,0.975,0.5)))
  cis<-c(data[1,],rev(data[2,]))
  x<-c(yreal[1:500],yreal[500:1])
  cis[cis<range[1]]<-range[1]
  cis[cis>range[2]]<-range[2]
  polygon(x,cis,border=NA,col=colb[xx],xpd=F)}

## 'surv_fit' function

## Fit curve of survival trend over the 50 years. Solid lines imply 
## that that the curve is increasing or decreasing with 
## at least 95% plausibility (95% credible interval for 1st derivation 
## of estimated survival does not cross zero) for at least 5 years. 
## Dashed lines show otherwise (no detectable change in the survival)

surv_fit<-function(dat2){
  dat2=data.frame(dat2)
  cis= data.frame((dat2[,-1] - dat2[,-ncol(dat2)])*10)
  cis=sapply(cis, function(p) quantile(p, probs = c(0.025,0.975)))
  est=sapply(dat2, function(p) quantile(p, probs = c(0.5)))
  data=cis
  cbinl<-c()
  x=1
  repeat{
    cbinl[x]<-
      if(data[1,x]>0|data[2,x]<0){cbinl[x]=1}else{cbinl[x]=0}
    x<-x+1
    if(x>length(data[1,])){break}}
  cbin=1;x<-1
  repeat{
    cbin[x+1]<-abs(cbinl[x]-cbinl[x+1])
    x=x+1
    if(x>(length(cbinl)-1)){break}}
  data<-data.frame(rbind(data,cbin));data[4,]<-yreal[c(2:500)]-0.049;data[5,]<-1:499
  row.names(data)<-c("cil","ciu","stat_change","year","timepoint");data[3,499]<-1
  data2=t(data[,data[3,]==1])
  data2<-data.frame(data2)
  tr<-data2;y=1
  yreal=seq(1971,2020,length=500)
  repeat{
    if((((tr[y,1]<0)&(tr[y,2])<0))|(((tr[y,1]>0)&(tr[y,2])>0))&(tr[y+1,4]-tr[y,4])>5)
    {lines(yreal[tr[y,5]:tr[y+1,5]],
           est[tr[y,5]:tr[y+1,5]],lwd=1.9,col=cola[xx],lend=1)}
    y<-y+1;if(y>dim(tr)[1]-1){break}}
  lines(yreal[1:500],est,lwd=1,col=cola[xx],lty=2,lend=1)}

## 'polyg_slope' function -----------

## Draws 95% credible interval for slope of the 
## survival trend (1st derivation of the estimated survival trend)
polyg_slope<-function(data){ 
  x<-c(yreal[1:499],yreal[499:1])
  data<-data.frame((data[,-1] - data[,-ncol(data)]))*10
  data=sapply(data, function(p) quantile(p, probs = c(0.025,0.975,0.5)))
  cis<-c(data[1,],rev(data[2,]))
  cis[cis<range[1]]<-range[1]
  cis[cis>range[2]]<-range[2]
  polygon(x,cis,border=NA,col=colb[xx])}

## 'slope_fit' function

## Fit curve of slope of the survival trend over the 50 years. 
## Solid lines imply that that the curve is increasing or decreasing with 
## at least 95% plausibility (95% credible interval for 2st derivation of 
## estimated survival does not cross zero) for at least 3 years. Dashed 
## lines show otherwise (no detectable change in the slope of the survival trend)

slope_fit<-function(data){
  data = data.frame((data[,-1] - data[,-ncol(data)])*10)
  dar2<-sapply(data, function(p) quantile(p, probs = c(0.5)))
  data<-data.frame(data[,-1] - data[,-ncol(data)])
  data=sapply(data, function(p) quantile(p, probs = c(0.025,0.975,0.5)))
  cbinl<-c()
  x=1
  repeat{
    cbinl[x]<-
      if(data[1,x]>0|data[2,x]<0){cbinl[x]=1}else{cbinl[x]=0}
    x<-x+1
    if(x>length(data[1,])){break}}
  cbin=1;x<-1
  repeat{
    cbin[x+1]<-abs(cbinl[x]-cbinl[x+1])
    x=x+1
    if(x>(length(cbinl)-1)){break}}
  data<-data.frame(rbind(data,cbin));data[5,]<-yreal[c(2:499)]-0.049;data[6,]<-1:498
  row.names(data)<-c("cil","ciu","est","stat_change","year","timepoint")
  data[4,450]<-1
  data[4,50]<-1
  data=data[,50:450]
  dar2[dar2<range[1]]<-range[1]
  lines(yreal[1:499],dar2,lwd=1,col=cola[xx],lty=2,lend=1)
  data2=(t(data[,data[4,]==1]))
  data2<-data.frame(data2)
  tr<-data2;y=1
  yreal=seq(1971,2020,length=498)
  repeat{
    if(((((tr[y,1]<0)&(tr[y,2])<0))|(((tr[y,1]>0)&(tr[y,2])>0)))&(tr[y+1,5]-tr[y,5])>3)
    {lines(yreal[tr[y,6]:tr[y+1,6]],
           dar2[tr[y,6]:tr[y+1,6]]
           ,lwd=1.9,col=cola[xx],lend=1)}
    y<-y+1;if(y>dim(tr)[1]-1){break}}
  return(data2)}




# Posterior extraction ----------------------------------------------------
## melanoma posterior extraction ---------------------------------------------------
yreal <- seq(1971,2020,length=500)
first <- expand.grid(years10cen = ((yreal-1995.5)/10),
                     group = c('Denmark.Females','Finland.Females','Norway.Females','Sweden.Females',
                               'Denmark.Males','Finland.Males','Norway.Males','Sweden.Males'), y = 0)

ms_1y<-posterior_smooths(melanoma_1y_model,smooth="s(years10cen,by=group,k=5)",newdata = first)
post_fix<-as.data.frame(melanoma_1y_model, variable = c("b_Intercept","b_groupFinland.Females",
                                                     "b_groupNorway.Females","b_groupSweden.Females","b_groupDenmark.Males",
                                                     "b_groupFinland.Males","b_groupNorway.Males","b_groupSweden.Males"))
fixef(melanoma_1y_model)
post_melanoma_den_fem_1y<-ms_1y[,1:500]     +post_fix[,1]
post_melanoma_fin_fem_1y<-ms_1y[,501:1000]  +post_fix[,1]+post_fix[,2]
post_melanoma_nor_fem_1y<-ms_1y[,1001:1500] +post_fix[,1]+post_fix[,3]
post_melanoma_swe_fem_1y<-ms_1y[,1501:2000] +post_fix[,1]+post_fix[,4]
post_melanoma_den_mal_1y<-ms_1y[,2001:2500] +post_fix[,1]+post_fix[,5]
post_melanoma_fin_mal_1y<-ms_1y[,2501:3000] +post_fix[,1]+post_fix[,6]
post_melanoma_nor_mal_1y<-ms_1y[,3001:3500] +post_fix[,1]+post_fix[,7]
post_melanoma_swe_mal_1y<-ms_1y[,3501:4000] +post_fix[,1]+post_fix[,8]

post_fix<-as.data.frame(melanoma_5y_model, variable = c("b_Intercept","b_groupFinland.Females",
                                                     "b_groupNorway.Females","b_groupSweden.Females","b_groupDenmark.Males",
                                                     "b_groupFinland.Males","b_groupNorway.Males","b_groupSweden.Males"))
ms_5y<-posterior_smooths(melanoma_5y_model,smooth="s(years10cen,by=group,k=5)",newdata = first)
post_melanoma_den_fem_5y<-ms_5y[,1:500]     +post_fix[,1]
post_melanoma_fin_fem_5y<-ms_5y[,501:1000]  +post_fix[,1]+post_fix[,2]
post_melanoma_nor_fem_5y<-ms_5y[,1001:1500] +post_fix[,1]+post_fix[,3]
post_melanoma_swe_fem_5y<-ms_5y[,1501:2000] +post_fix[,1]+post_fix[,4]
post_melanoma_den_mal_5y<-ms_5y[,2001:2500] +post_fix[,1]+post_fix[,5]
post_melanoma_fin_mal_5y<-ms_5y[,2501:3000] +post_fix[,1]+post_fix[,6]
post_melanoma_nor_mal_5y<-ms_5y[,3001:3500] +post_fix[,1]+post_fix[,7]
post_melanoma_swe_mal_5y<-ms_5y[,3501:4000] +post_fix[,1]+post_fix[,8]

post_melanoma_den_fem_cond<-(post_melanoma_den_fem_5y/post_melanoma_den_fem_1y)*100
post_melanoma_den_mal_cond<-(post_melanoma_den_mal_5y/post_melanoma_den_mal_1y)*100
post_melanoma_fin_fem_cond<-(post_melanoma_fin_fem_5y/post_melanoma_fin_fem_1y)*100
post_melanoma_fin_mal_cond<-(post_melanoma_fin_mal_5y/post_melanoma_fin_mal_1y)*100
post_melanoma_nor_fem_cond<-(post_melanoma_nor_fem_5y/post_melanoma_nor_fem_1y)*100
post_melanoma_nor_mal_cond<-(post_melanoma_nor_mal_5y/post_melanoma_nor_mal_1y)*100
post_melanoma_swe_fem_cond<-(post_melanoma_swe_fem_5y/post_melanoma_swe_fem_1y)*100
post_melanoma_swe_mal_cond<-(post_melanoma_swe_mal_5y/post_melanoma_swe_mal_1y)*100


# -----------

  
# Plotting 1 - survival trends ----------------


## general setting -------------
m= matrix(c( 25, 1, 2, 3, 4
              ,5, 9,13,17,21
              ,6,10,14,18,22
              ,7,11,15,19,23
              ,8,12,16,20,24), nrow = 5, ncol=5, byrow = TRUE)
layout(mat = m,heights = c(0.03,0.6/2,0.37/2,0.6/2,0.37/2),
       widths = c(0.04,rep(0.96/4,4)))
par(mgp=c(1.6,0.62,0))
par(mar=c(0,0,0,0))


plot(NULL, axes=FALSE,xlab="",ylab="",xlim=c(-1,1),ylim=c(-0.85,0.85))
text(0,-0.2,"Denmark",cex=1.6,font=3,xpd=TRUE)

plot(NULL, axes=FALSE,xlab="",ylab="",xlim=c(-1,1),ylim=c(-0.85,0.85))
text(0,-0.2,"Finland",cex=1.6,font=3,xpd=TRUE)

plot(NULL, axes=FALSE,xlab="",ylab="",xlim=c(-1,1),ylim=c(-0.85,0.85))
text(0,-0.2,"Norway",cex=1.6,font=3,xpd=TRUE)

plot(NULL, axes=FALSE,xlab="",ylab="",xlim=c(-1,1),ylim=c(-0.85,0.85))
text(0,-0.2,"Sweden",cex=1.6,font=3,xpd=TRUE)



range_b<-c(40,110);scal_b<-range_b[2]-range_b[1]
range_c<-c(0,2.2);scal_c<-range_c[2]-range_c[1]

range<-range_b;scal=scal_b
plot(NULL, axes=FALSE,xlab="",ylab="",xlim=c(-1,1),
     ylim=c(range_b[1],range_b[2]))
text(0,range_b[1]+scal_b*0.5,"Relative survival (%) in males", cex=1.35,srt=90)


range<-range_c;scal=scal_c
par(mar=c(2.5,0,0,0))
plot(NULL, axes=FALSE,xlab="",ylab="",xlim=c(-1,1),
     ylim=c(range_c[1],range_c[2]))
text(0,range_c[1]+scal_c*0.5,expression(paste(delta, " (survival)")),
     cex=1.35,srt=90)

range<-range_b;scal=scal_b
par(mar=c(0,0,0,0))
plot(NULL, axes=FALSE,xlab="",ylab="",xlim=c(-1,1),
     ylim=c(range_b[1],range_b[2]))
text(0,range_b[1]+scal_b*0.5,"Relative survival (%) in females", cex=1.4,srt=90)


range<-range_c;scal=scal_c
par(mar=c(2.5,0,0.2,0))
plot(NULL, axes=FALSE,xlab="",ylab="",xlim=c(-1,1),
     ylim=c(range_c[1],range_c[2]))
text(0,range_c[1]+scal_c*0.5,expression(paste(delta, " (survival)")),
     cex=1.4,srt=90)


### DEN  % ------------------------------------------------------
#### males --------------
range<-range_b;scal=scal_b
par(mgp=c(1.6,0.4,0))
par(mar=c(0,1.4,0,0.3))

plot(NULL,xlim=c(1971,2020),ylim=c(range[1],range[2]),xlab="",ylab=""
     ,las=1, axes=FALSE)
rect(1971,range[2],2020,range[1],col="grey90",border=NA)
x<-10
repeat{
  lines(c(1971,2020),c(x,x),col="white",lwd=0.7)
  x=x+10;if(x>100){break}}
lines(c(1971,2020),c(50,50),col="white",lwd=1.7)


x<-1980
repeat{
  lines(c(x,x),c(range[1],range[2]),col="white",lwd=0.7)
  x=x+10;if(x>2020){break}}

# breaking points identification
xx=1
breakpo(post_melanoma_den_mal_1y,0)
xx=xx+1
breakpo(post_melanoma_den_mal_cond,-0.7)
xx=xx+1
breakpo(post_melanoma_den_mal_5y,0.7)

# 95% credible interval for survival
xx=1
polyg_surv(post_melanoma_den_mal_1y);xx=xx+1
polyg_surv(post_melanoma_den_mal_cond);xx=xx+1
polyg_surv(post_melanoma_den_mal_5y)

# data points of estimated survival
xx=1
points(melanoma[melanoma$shou=="den_mal_",]$surv_1y~melanoma[melanoma$shou=="den_mal_",]$year
       ,pch=17,col=colc[xx],cex=1);xx=xx+1
points(melanoma[melanoma$shou=="den_mal_",]$surv_cond~melanoma[melanoma$shou=="den_mal_",]$year
       ,pch=17,col=colc[xx],cex=1);xx=xx+1
points(melanoma[melanoma$shou=="den_mal_",]$surv_5y~melanoma[melanoma$shou=="den_mal_",]$year
       ,pch=17,col=colc[xx],cex=1)

tckk=-0.016
# fitted lines for survival
xx=1
surv_fit(post_melanoma_den_mal_1y);xx=xx+1
surv_fit(post_melanoma_den_mal_cond);xx=xx+1
surv_fit(post_melanoma_den_mal_5y);xx=xx+1

axis(2,las=2,cex.axis=1.3,at=seq(range[1],range[2],by=10),
     labels=c(rep("",length(seq(range[1],range[2],by=10)))),pos=1971,tck=tckk)

axis(2,las=2,cex.axis=1.3,at=c(seq(range[1],range[2],by=20)),pos=1971,tck=tckk)

axis(side=1,las=1,cex.axis=1.3,at=c(seq(1980,2020,by=10)),pos=range[1],
     labels = c(rep("",5)),tck=tckk)
lines(c(1971,2020),c(range[1],range[1]))

# legend
xpo=-8
ypo=0.3
xx=1;yy=range[1]+scal*ypo
rect(2004.3+xpo,yy+0.035*scal,2024+xpo,yy-0.18*scal,col="white",border="grey50",lwd=0.8)
repeat{
  points(2007+xpo,yy,pch=17,col=cola[xx],cex=1.2)
  lines(c(2005+xpo,2009+xpo),c(yy,yy),lwd=1.6,col=cola[xx],lend=1)
  xx<-xx+1;yy=yy-(scal*0.07);if(xx>3){break}}
xx=1;yy=range[1]+scal*ypo
text(2016.7+xpo,yy,"1-year",col=cola[xx],cex=1.2);xx=xx+1;yy=yy-(scal*0.07)
text(2016.7+xpo,yy,"5/1-year",col=cola[xx],cex=1.2);xx=xx+1;yy=yy-(scal*0.07)
text(2016.7+xpo,yy,"5-year",col=cola[xx],cex=1.2);xx=xx+1;yy=yy-(scal*0.07)
text(1974,range[2]-0.05*scal,"a",cex=2.2)

##### slope plot --------------
range<-range_c;scal=scal_c
par(mar=c(2.5,1.4,0,0))

plot(NULL,xlim=c(1971,2020),ylim=c(range[1],range[2]),xlab="",ylab=""
     ,las=1, axes=FALSE)
rect(1971,range[2],2020,range[1],col="grey90",border=NA)
x<-range[1]
repeat{
  lines(c(1971,2020),c(x,x),col="white",lwd=0.7)
  x=x+0.5;if(x>100){break}}

x<-1980
repeat{
  lines(c(x,x),c(range[1],range[2]),col="white",lwd=0.7)
  x=x+10;if(x>2020){break}}

xx=1
polyg_slope(post_melanoma_den_mal_1y);xx=xx+1
polyg_slope(post_melanoma_den_mal_cond);xx=xx+1
polyg_slope(post_melanoma_den_mal_5y)

xx=1
slope_fit(post_melanoma_den_mal_1y);xx=xx+1
slope_fit(post_melanoma_den_mal_cond);xx=xx+1
slope_fit(post_melanoma_den_mal_5y)

axis(2,las=2,cex.axis=1.3,at=c(seq(range[1],range[2],by=0.5)),
     labels=c(rep("",length(seq(range[1],range[2],by=0.5)))),pos=1971,tck=tckk*2)
axis(2,las=2,cex.axis=1.3,at=c(seq(-1,range[2],by=1)),pos=1971,tck=tckk*2)
lines(c(1971,2020),c(0,0),col="grey50")

axis(side=1,las=1,cex.axis=1.3,at=c(seq(1980,2020,by=10)),pos=range[1],labels=c(rep("",5)),tck=tckk*2)
text(c(seq(1977.3,2020,by=20)),c(rep(range[1]-0.1*scal,5)),c(seq(1980,2020,by=20)),xpd=TRUE,cex=1.4);
lines(c(1971,2020),c(range[1],range[1]))
title(xlab="Year", line=1.4, cex.lab=1.4,xpd=TRUE)



#### females --------------
range<-range_b;scal=scal_b
par(mgp=c(1.6,0.4,0))
par(mar=c(0,1.4,0,0.3))

plot(NULL,xlim=c(1971,2020),ylim=c(range[1],range[2]),xlab="",ylab=""
     ,las=1, axes=FALSE)
rect(1971,range[2],2020,range[1],col="grey90",border=NA)
x<-10
repeat{
  lines(c(1971,2020),c(x,x),col="white",lwd=0.7)
  x=x+10;if(x>100){break}}
lines(c(1971,2020),c(50,50),col="white",lwd=1.7)


x<-1980
repeat{
  lines(c(x,x),c(range[1],range[2]),col="white",lwd=0.7)
  x=x+10;if(x>2020){break}}

# breaking points identification
xx=1
breakpo(post_melanoma_den_fem_1y,0)
xx=xx+1
breakpo(post_melanoma_den_fem_cond,1)
xx=xx+1
breakpo(post_melanoma_den_fem_5y,0)

# 95% credible interval for survival
xx=1
polyg_surv(post_melanoma_den_fem_1y);xx=xx+1
polyg_surv(post_melanoma_den_fem_cond);xx=xx+1
polyg_surv(post_melanoma_den_fem_5y)

# data points of estimated survival
xx=1
points(melanoma[melanoma$shou=="den_fem_",]$surv_1y~melanoma[melanoma$shou=="den_fem_",]$year
       ,pch=17,col=colc[xx],cex=1);xx=xx+1
points(melanoma[melanoma$shou=="den_fem_",]$surv_cond~melanoma[melanoma$shou=="den_fem_",]$year
       ,pch=17,col=colc[xx],cex=1);xx=xx+1
points(melanoma[melanoma$shou=="den_fem_",]$surv_5y~melanoma[melanoma$shou=="den_fem_",]$year
       ,pch=17,col=colc[xx],cex=1)

tckk=-0.016
# fitted lines for survival
xx=1
surv_fit(post_melanoma_den_fem_1y);xx=xx+1
surv_fit(post_melanoma_den_fem_cond);xx=xx+1
surv_fit(post_melanoma_den_fem_5y);xx=xx+1

axis(2,las=2,cex.axis=1.3,at=seq(range[1],range[2],by=10),
     labels=c(rep("",length(seq(range[1],range[2],by=10)))),pos=1971,tck=tckk)

axis(2,las=2,cex.axis=1.3,at=c(seq(range[1],range[2],by=20)),pos=1971,tck=tckk)

axis(side=1,las=1,cex.axis=1.3,at=c(seq(1980,2020,by=10)),pos=range[1],
     labels = c(rep("",5)),tck=tckk)
lines(c(1971,2020),c(range[1],range[1]))
text(1974,range[2]-0.05*scal,"e",cex=2.2)

##### slope plot --------------
range<-range_c;scal=scal_c
par(mar=c(2.5,1.4,0,0))

plot(NULL,xlim=c(1971,2020),ylim=c(range[1],range[2]),xlab="",ylab=""
     ,las=1, axes=FALSE)
rect(1971,range[2],2020,range[1],col="grey90",border=NA)
x<-range[1]
repeat{
  lines(c(1971,2020),c(x,x),col="white",lwd=0.7)
  x=x+0.5;if(x>100){break}}

x<-1980
repeat{
  lines(c(x,x),c(range[1],range[2]),col="white",lwd=0.7)
  x=x+10;if(x>2020){break}}

xx=1
polyg_slope(post_melanoma_den_fem_1y);xx=xx+1
polyg_slope(post_melanoma_den_fem_cond);xx=xx+1
polyg_slope(post_melanoma_den_fem_5y)

xx=1
slope_fit(post_melanoma_den_fem_1y);xx=xx+1
slope_fit(post_melanoma_den_fem_cond);xx=xx+1
slope_fit(post_melanoma_den_fem_5y)

axis(2,las=2,cex.axis=1.3,at=c(seq(range[1],range[2],by=0.5)),
     labels=c(rep("",length(seq(range[1],range[2],by=0.5)))),pos=1971,tck=tckk*2)
axis(2,las=2,cex.axis=1.3,at=c(seq(-1,range[2],by=1)),pos=1971,tck=tckk*2)
lines(c(1971,2020),c(0,0),col="grey50")

axis(side=1,las=1,cex.axis=1.3,at=c(seq(1980,2020,by=10)),pos=range[1],labels=c(rep("",5)),tck=tckk*2)
text(c(seq(1977.3,2020,by=20)),c(rep(range[1]-0.1*scal,5)),c(seq(1980,2020,by=20)),xpd=TRUE,cex=1.4);
lines(c(1971,2020),c(range[1],range[1]))
title(xlab="Year", line=1.4, cex.lab=1.4,xpd=TRUE)



### FIN  % ------------------------------------------------------
#### males --------------
range<-range_b;scal=scal_b
par(mgp=c(1.6,0.4,0))
par(mar=c(0,1.4,0,0.3))

plot(NULL,xlim=c(1971,2020),ylim=c(range[1],range[2]),xlab="",ylab=""
     ,las=1, axes=FALSE)
rect(1971,range[2],2020,range[1],col="grey90",border=NA)
x<-10
repeat{
  lines(c(1971,2020),c(x,x),col="white",lwd=0.7)
  x=x+10;if(x>100){break}}
lines(c(1971,2020),c(50,50),col="white",lwd=1.7)


x<-1980
repeat{
  lines(c(x,x),c(range[1],range[2]),col="white",lwd=0.7)
  x=x+10;if(x>2020){break}}

# breaking points identification
xx=1
breakpo(post_melanoma_fin_mal_1y,0)
xx=xx+1
breakpo(post_melanoma_fin_mal_cond,-0.4)
xx=xx+1
breakpo(post_melanoma_fin_mal_5y,0.4)

# 95% credible interval for survival
xx=1
polyg_surv(post_melanoma_fin_mal_1y);xx=xx+1
polyg_surv(post_melanoma_fin_mal_cond);xx=xx+1
polyg_surv(post_melanoma_fin_mal_5y)

# data points of estimated survival
xx=1
points(melanoma[melanoma$shou=="fin_mal_",]$surv_1y~melanoma[melanoma$shou=="fin_mal_",]$year
       ,pch=17,col=colc[xx],cex=1);xx=xx+1
points(melanoma[melanoma$shou=="fin_mal_",]$surv_cond~melanoma[melanoma$shou=="fin_mal_",]$year
       ,pch=17,col=colc[xx],cex=1);xx=xx+1
points(melanoma[melanoma$shou=="fin_mal_",]$surv_5y~melanoma[melanoma$shou=="fin_mal_",]$year
       ,pch=17,col=colc[xx],cex=1)

tckk=-0.016
# fitted lines for survival
xx=1
surv_fit(post_melanoma_fin_mal_1y);xx=xx+1
surv_fit(post_melanoma_fin_mal_cond);xx=xx+1
surv_fit(post_melanoma_fin_mal_5y);xx=xx+1

axis(2,las=2,cex.axis=1.3,at=seq(range[1],range[2],by=10),
     labels=c(rep("",length(seq(range[1],range[2],by=10)))),pos=1971,tck=tckk)

axis(2,las=2,cex.axis=1.3,at=c(seq(range[1],range[2],by=20)),pos=1971,tck=tckk)

axis(side=1,las=1,cex.axis=1.3,at=c(seq(1980,2020,by=10)),pos=range[1],
     labels = c(rep("",5)),tck=tckk)
lines(c(1971,2020),c(range[1],range[1]))

text(1974,range[2]-0.05*scal,"b",cex=2.2)

##### slope plot --------------
range<-range_c;scal=scal_c
par(mar=c(2.5,1.4,0,0))

plot(NULL,xlim=c(1971,2020),ylim=c(range[1],range[2]),xlab="",ylab=""
     ,las=1, axes=FALSE)
rect(1971,range[2],2020,range[1],col="grey90",border=NA)
x<-range[1]
repeat{
  lines(c(1971,2020),c(x,x),col="white",lwd=0.7)
  x=x+0.5;if(x>100){break}}

x<-1980
repeat{
  lines(c(x,x),c(range[1],range[2]),col="white",lwd=0.7)
  x=x+10;if(x>2020){break}}

xx=1
polyg_slope(post_melanoma_fin_mal_1y);xx=xx+1
polyg_slope(post_melanoma_fin_mal_cond);xx=xx+1
polyg_slope(post_melanoma_fin_mal_5y)

xx=1
slope_fit(post_melanoma_fin_mal_1y);xx=xx+1
slope_fit(post_melanoma_fin_mal_cond);xx=xx+1
slope_fit(post_melanoma_fin_mal_5y)

axis(2,las=2,cex.axis=1.3,at=c(seq(range[1],range[2],by=0.5)),
     labels=c(rep("",length(seq(range[1],range[2],by=0.5)))),pos=1971,tck=tckk*2)
axis(2,las=2,cex.axis=1.3,at=c(seq(-1,range[2],by=1)),pos=1971,tck=tckk*2)
lines(c(1971,2020),c(0,0),col="grey50")

axis(side=1,las=1,cex.axis=1.3,at=c(seq(1980,2020,by=10)),pos=range[1],labels=c(rep("",5)),tck=tckk*2)
text(c(seq(1977.3,2020,by=20)),c(rep(range[1]-0.1*scal,5)),c(seq(1980,2020,by=20)),xpd=TRUE,cex=1.4);
lines(c(1971,2020),c(range[1],range[1]))
title(xlab="Year", line=1.4, cex.lab=1.4,xpd=TRUE)



#### females --------------
range<-range_b;scal=scal_b
par(mgp=c(1.6,0.4,0))
par(mar=c(0,1.4,0,0.3))

plot(NULL,xlim=c(1971,2020),ylim=c(range[1],range[2]),xlab="",ylab=""
     ,las=1, axes=FALSE)
rect(1971,range[2],2020,range[1],col="grey90",border=NA)
x<-10
repeat{
  lines(c(1971,2020),c(x,x),col="white",lwd=0.7)
  x=x+10;if(x>100){break}}
lines(c(1971,2020),c(50,50),col="white",lwd=1.7)


x<-1980
repeat{
  lines(c(x,x),c(range[1],range[2]),col="white",lwd=0.7)
  x=x+10;if(x>2020){break}}

# breaking points identification
xx=1
breakpo(post_melanoma_fin_fem_1y,0)
xx=xx+1
breakpo(post_melanoma_fin_fem_cond,0.7)
xx=xx+1
breakpo(post_melanoma_fin_fem_5y,-0.7)

# 95% credible interval for survival
xx=1
polyg_surv(post_melanoma_fin_fem_1y);xx=xx+1
polyg_surv(post_melanoma_fin_fem_cond);xx=xx+1
polyg_surv(post_melanoma_fin_fem_5y)

# data points of estimated survival
xx=1
points(melanoma[melanoma$shou=="fin_fem_",]$surv_1y~melanoma[melanoma$shou=="fin_fem_",]$year
       ,pch=17,col=colc[xx],cex=1);xx=xx+1
points(melanoma[melanoma$shou=="fin_fem_",]$surv_cond~melanoma[melanoma$shou=="fin_fem_",]$year
       ,pch=17,col=colc[xx],cex=1);xx=xx+1
points(melanoma[melanoma$shou=="fin_fem_",]$surv_5y~melanoma[melanoma$shou=="fin_fem_",]$year
       ,pch=17,col=colc[xx],cex=1)

tckk=-0.016
# fitted lines for survival
xx=1
surv_fit(post_melanoma_fin_fem_1y);xx=xx+1
surv_fit(post_melanoma_fin_fem_cond);xx=xx+1
surv_fit(post_melanoma_fin_fem_5y);xx=xx+1

axis(2,las=2,cex.axis=1.3,at=seq(range[1],range[2],by=10),
     labels=c(rep("",length(seq(range[1],range[2],by=10)))),pos=1971,tck=tckk)

axis(2,las=2,cex.axis=1.3,at=c(seq(range[1],range[2],by=20)),pos=1971,tck=tckk)

axis(side=1,las=1,cex.axis=1.3,at=c(seq(1980,2020,by=10)),pos=range[1],
     labels = c(rep("",5)),tck=tckk)
lines(c(1971,2020),c(range[1],range[1]))
text(1974,range[2]-0.05*scal,"f",cex=2.2)

##### slope plot --------------
range<-range_c;scal=scal_c
par(mar=c(2.5,1.4,0,0))

plot(NULL,xlim=c(1971,2020),ylim=c(range[1],range[2]),xlab="",ylab=""
     ,las=1, axes=FALSE)
rect(1971,range[2],2020,range[1],col="grey90",border=NA)
x<-range[1]
repeat{
  lines(c(1971,2020),c(x,x),col="white",lwd=0.7)
  x=x+0.5;if(x>100){break}}

x<-1980
repeat{
  lines(c(x,x),c(range[1],range[2]),col="white",lwd=0.7)
  x=x+10;if(x>2020){break}}

xx=1
polyg_slope(post_melanoma_fin_fem_1y);xx=xx+1
polyg_slope(post_melanoma_fin_fem_cond);xx=xx+1
polyg_slope(post_melanoma_fin_fem_5y)

xx=1
slope_fit(post_melanoma_fin_fem_1y);xx=xx+1
slope_fit(post_melanoma_fin_fem_cond);xx=xx+1
slope_fit(post_melanoma_fin_fem_5y)

axis(2,las=2,cex.axis=1.3,at=c(seq(range[1],range[2],by=0.5)),
     labels=c(rep("",length(seq(range[1],range[2],by=0.5)))),pos=1971,tck=tckk*2)
axis(2,las=2,cex.axis=1.3,at=c(seq(-1,range[2],by=1)),pos=1971,tck=tckk*2)
lines(c(1971,2020),c(0,0),col="grey50")

axis(side=1,las=1,cex.axis=1.3,at=c(seq(1980,2020,by=10)),pos=range[1],labels=c(rep("",5)),tck=tckk*2)
text(c(seq(1977.3,2020,by=20)),c(rep(range[1]-0.1*scal,5)),c(seq(1980,2020,by=20)),xpd=TRUE,cex=1.4);
lines(c(1971,2020),c(range[1],range[1]))
title(xlab="Year", line=1.4, cex.lab=1.4,xpd=TRUE)


### NOR  % ------------------------------------------------------
#### males --------------
range<-range_b;scal=scal_b
par(mgp=c(1.6,0.4,0))
par(mar=c(0,1.4,0,0.3))

plot(NULL,xlim=c(1971,2020),ylim=c(range[1],range[2]),xlab="",ylab=""
     ,las=1, axes=FALSE)
rect(1971,range[2],2020,range[1],col="grey90",border=NA)
x<-10
repeat{
  lines(c(1971,2020),c(x,x),col="white",lwd=0.7)
  x=x+10;if(x>100){break}}
lines(c(1971,2020),c(50,50),col="white",lwd=1.7)


x<-1980
repeat{
  lines(c(x,x),c(range[1],range[2]),col="white",lwd=0.7)
  x=x+10;if(x>2020){break}}

# breaking points identification
xx=1
breakpo(post_melanoma_nor_mal_1y,0)
xx=xx+1
breakpo(post_melanoma_nor_mal_cond,-0.6)
xx=xx+1
breakpo(post_melanoma_nor_mal_5y,1)

# 95% credible interval for survival
xx=1
polyg_surv(post_melanoma_nor_mal_1y);xx=xx+1
polyg_surv(post_melanoma_nor_mal_cond);xx=xx+1
polyg_surv(post_melanoma_nor_mal_5y)

# data points of estimated survival
xx=1
points(melanoma[melanoma$shou=="nor_mal_",]$surv_1y~melanoma[melanoma$shou=="nor_mal_",]$year
       ,pch=17,col=colc[xx],cex=1);xx=xx+1
points(melanoma[melanoma$shou=="nor_mal_",]$surv_cond~melanoma[melanoma$shou=="nor_mal_",]$year
       ,pch=17,col=colc[xx],cex=1);xx=xx+1
points(melanoma[melanoma$shou=="nor_mal_",]$surv_5y~melanoma[melanoma$shou=="nor_mal_",]$year
       ,pch=17,col=colc[xx],cex=1)

tckk=-0.016
# fitted lines for survival
xx=1
surv_fit(post_melanoma_nor_mal_1y);xx=xx+1
surv_fit(post_melanoma_nor_mal_cond);xx=xx+1
surv_fit(post_melanoma_nor_mal_5y);xx=xx+1

axis(2,las=2,cex.axis=1.3,at=seq(range[1],range[2],by=10),
     labels=c(rep("",length(seq(range[1],range[2],by=10)))),pos=1971,tck=tckk)

axis(2,las=2,cex.axis=1.3,at=c(seq(range[1],range[2],by=20)),pos=1971,tck=tckk)

axis(side=1,las=1,cex.axis=1.3,at=c(seq(1980,2020,by=10)),pos=range[1],
     labels = c(rep("",5)),tck=tckk)
lines(c(1971,2020),c(range[1],range[1]))

text(1974,range[2]-0.05*scal,"c",cex=2.2)

##### slope plot --------------
range<-range_c;scal=scal_c
par(mar=c(2.5,1.4,0,0))

plot(NULL,xlim=c(1971,2020),ylim=c(range[1],range[2]),xlab="",ylab=""
     ,las=1, axes=FALSE)
rect(1971,range[2],2020,range[1],col="grey90",border=NA)
x<-range[1]
repeat{
  lines(c(1971,2020),c(x,x),col="white",lwd=0.7)
  x=x+0.5;if(x>100){break}}

x<-1980
repeat{
  lines(c(x,x),c(range[1],range[2]),col="white",lwd=0.7)
  x=x+10;if(x>2020){break}}

xx=1
polyg_slope(post_melanoma_nor_mal_1y);xx=xx+1
polyg_slope(post_melanoma_nor_mal_cond);xx=xx+1
polyg_slope(post_melanoma_nor_mal_5y)

xx=1
slope_fit(post_melanoma_nor_mal_1y);xx=xx+1
slope_fit(post_melanoma_nor_mal_cond);xx=xx+1
slope_fit(post_melanoma_nor_mal_5y)

axis(2,las=2,cex.axis=1.3,at=c(seq(range[1],range[2],by=0.5)),
     labels=c(rep("",length(seq(range[1],range[2],by=0.5)))),pos=1971,tck=tckk*2)
axis(2,las=2,cex.axis=1.3,at=c(seq(-1,range[2],by=1)),pos=1971,tck=tckk*2)
lines(c(1971,2020),c(0,0),col="grey50")

axis(side=1,las=1,cex.axis=1.3,at=c(seq(1980,2020,by=10)),pos=range[1],labels=c(rep("",5)),tck=tckk*2)
text(c(seq(1977.3,2020,by=20)),c(rep(range[1]-0.1*scal,5)),c(seq(1980,2020,by=20)),xpd=TRUE,cex=1.4);
lines(c(1971,2020),c(range[1],range[1]))
title(xlab="Year", line=1.4, cex.lab=1.4,xpd=TRUE)



#### females --------------
range<-range_b;scal=scal_b
par(mgp=c(1.6,0.4,0))
par(mar=c(0,1.4,0,0.3))

plot(NULL,xlim=c(1971,2020),ylim=c(range[1],range[2]),xlab="",ylab=""
     ,las=1, axes=FALSE)
rect(1971,range[2],2020,range[1],col="grey90",border=NA)
x<-10
repeat{
  lines(c(1971,2020),c(x,x),col="white",lwd=0.7)
  x=x+10;if(x>100){break}}
lines(c(1971,2020),c(50,50),col="white",lwd=1.7)


x<-1980
repeat{
  lines(c(x,x),c(range[1],range[2]),col="white",lwd=0.7)
  x=x+10;if(x>2020){break}}

# breaking points identification
xx=1
breakpo(post_melanoma_nor_fem_1y,0)
xx=xx+1
breakpo(post_melanoma_nor_fem_cond,0)
xx=xx+1
breakpo(post_melanoma_nor_fem_5y,1)

# 95% credible interval for survival
xx=1
polyg_surv(post_melanoma_nor_fem_1y);xx=xx+1
polyg_surv(post_melanoma_nor_fem_cond);xx=xx+1
polyg_surv(post_melanoma_nor_fem_5y)

# data points of estimated survival
xx=1
points(melanoma[melanoma$shou=="nor_fem_",]$surv_1y~melanoma[melanoma$shou=="nor_fem_",]$year
       ,pch=17,col=colc[xx],cex=1);xx=xx+1
points(melanoma[melanoma$shou=="nor_fem_",]$surv_cond~melanoma[melanoma$shou=="nor_fem_",]$year
       ,pch=17,col=colc[xx],cex=1);xx=xx+1
points(melanoma[melanoma$shou=="nor_fem_",]$surv_5y~melanoma[melanoma$shou=="nor_fem_",]$year
       ,pch=17,col=colc[xx],cex=1)

tckk=-0.016
# fitted lines for survival
xx=1
surv_fit(post_melanoma_nor_fem_1y);xx=xx+1
surv_fit(post_melanoma_nor_fem_cond);xx=xx+1
surv_fit(post_melanoma_nor_fem_5y);xx=xx+1

axis(2,las=2,cex.axis=1.3,at=seq(range[1],range[2],by=10),
     labels=c(rep("",length(seq(range[1],range[2],by=10)))),pos=1971,tck=tckk)

axis(2,las=2,cex.axis=1.3,at=c(seq(range[1],range[2],by=20)),pos=1971,tck=tckk)

axis(side=1,las=1,cex.axis=1.3,at=c(seq(1980,2020,by=10)),pos=range[1],
     labels = c(rep("",5)),tck=tckk)
lines(c(1971,2020),c(range[1],range[1]))
text(1974,range[2]-0.05*scal,"g",cex=2.2)

##### slope plot --------------
range<-range_c;scal=scal_c
par(mar=c(2.5,1.4,0,0))

plot(NULL,xlim=c(1971,2020),ylim=c(range[1],range[2]),xlab="",ylab=""
     ,las=1, axes=FALSE)
rect(1971,range[2],2020,range[1],col="grey90",border=NA)
x<-range[1]
repeat{
  lines(c(1971,2020),c(x,x),col="white",lwd=0.7)
  x=x+0.5;if(x>100){break}}

x<-1980
repeat{
  lines(c(x,x),c(range[1],range[2]),col="white",lwd=0.7)
  x=x+10;if(x>2020){break}}

xx=1
polyg_slope(post_melanoma_nor_fem_1y);xx=xx+1
polyg_slope(post_melanoma_nor_fem_cond);xx=xx+1
polyg_slope(post_melanoma_nor_fem_5y)

xx=1
slope_fit(post_melanoma_nor_fem_1y);xx=xx+1
slope_fit(post_melanoma_nor_fem_cond);xx=xx+1
slope_fit(post_melanoma_nor_fem_5y)

axis(2,las=2,cex.axis=1.3,at=c(seq(range[1],range[2],by=0.5)),
     labels=c(rep("",length(seq(range[1],range[2],by=0.5)))),pos=1971,tck=tckk*2)
axis(2,las=2,cex.axis=1.3,at=c(seq(-1,range[2],by=1)),pos=1971,tck=tckk*2)
lines(c(1971,2020),c(0,0),col="grey50")

axis(side=1,las=1,cex.axis=1.3,at=c(seq(1980,2020,by=10)),pos=range[1],labels=c(rep("",5)),tck=tckk*2)
text(c(seq(1977.3,2020,by=20)),c(rep(range[1]-0.1*scal,5)),c(seq(1980,2020,by=20)),xpd=TRUE,cex=1.4);
lines(c(1971,2020),c(range[1],range[1]))
title(xlab="Year", line=1.4, cex.lab=1.4,xpd=TRUE)


### SWE  % ------------------------------------------------------
#### males --------------
range<-range_b;scal=scal_b
par(mgp=c(1.6,0.4,0))
par(mar=c(0,1.4,0,0.3))

plot(NULL,xlim=c(1971,2020),ylim=c(range[1],range[2]),xlab="",ylab=""
     ,las=1, axes=FALSE)
rect(1971,range[2],2020,range[1],col="grey90",border=NA)
x<-10
repeat{
  lines(c(1971,2020),c(x,x),col="white",lwd=0.7)
  x=x+10;if(x>100){break}}
lines(c(1971,2020),c(50,50),col="white",lwd=1.7)


x<-1980
repeat{
  lines(c(x,x),c(range[1],range[2]),col="white",lwd=0.7)
  x=x+10;if(x>2020){break}}

# breaking points identification
xx=1
breakpo(post_melanoma_swe_mal_1y,-0.2)
xx=xx+1
breakpo(post_melanoma_swe_mal_cond,-0.2)
xx=xx+1
breakpo(post_melanoma_swe_mal_5y,1.2)

# 95% credible interval for survival
xx=1
polyg_surv(post_melanoma_swe_mal_1y);xx=xx+1
polyg_surv(post_melanoma_swe_mal_cond);xx=xx+1
polyg_surv(post_melanoma_swe_mal_5y)

# data points of estimated survival
xx=1
points(melanoma[melanoma$shou=="swe_mal_",]$surv_1y~melanoma[melanoma$shou=="swe_mal_",]$year
       ,pch=17,col=colc[xx],cex=1);xx=xx+1
points(melanoma[melanoma$shou=="swe_mal_",]$surv_cond~melanoma[melanoma$shou=="swe_mal_",]$year
       ,pch=17,col=colc[xx],cex=1);xx=xx+1
points(melanoma[melanoma$shou=="swe_mal_",]$surv_5y~melanoma[melanoma$shou=="swe_mal_",]$year
       ,pch=17,col=colc[xx],cex=1)

tckk=-0.016
# fitted lines for survival
xx=1
surv_fit(post_melanoma_swe_mal_1y);xx=xx+1
surv_fit(post_melanoma_swe_mal_cond);xx=xx+1
surv_fit(post_melanoma_swe_mal_5y);xx=xx+1

axis(2,las=2,cex.axis=1.3,at=seq(range[1],range[2],by=10),
     labels=c(rep("",length(seq(range[1],range[2],by=10)))),pos=1971,tck=tckk)

axis(2,las=2,cex.axis=1.3,at=c(seq(range[1],range[2],by=20)),pos=1971,tck=tckk)

axis(side=1,las=1,cex.axis=1.3,at=c(seq(1980,2020,by=10)),pos=range[1],
     labels = c(rep("",5)),tck=tckk)
lines(c(1971,2020),c(range[1],range[1]))

text(1974,range[2]-0.05*scal,"d",cex=2.2)

##### slope plot --------------
range<-range_c;scal=scal_c
par(mar=c(2.5,1.4,0,0))

plot(NULL,xlim=c(1971,2020),ylim=c(range[1],range[2]),xlab="",ylab=""
     ,las=1, axes=FALSE)
rect(1971,range[2],2020,range[1],col="grey90",border=NA)
x<-range[1]
repeat{
  lines(c(1971,2020),c(x,x),col="white",lwd=0.7)
  x=x+0.5;if(x>100){break}}

x<-1980
repeat{
  lines(c(x,x),c(range[1],range[2]),col="white",lwd=0.7)
  x=x+10;if(x>2020){break}}

xx=1
polyg_slope(post_melanoma_swe_mal_1y);xx=xx+1
polyg_slope(post_melanoma_swe_mal_cond);xx=xx+1
polyg_slope(post_melanoma_swe_mal_5y)

xx=1
slope_fit(post_melanoma_swe_mal_1y);xx=xx+1
slope_fit(post_melanoma_swe_mal_cond);xx=xx+1
slope_fit(post_melanoma_swe_mal_5y)

axis(2,las=2,cex.axis=1.3,at=c(seq(range[1],range[2],by=0.5)),
     labels=c(rep("",length(seq(range[1],range[2],by=0.5)))),pos=1971,tck=tckk*2)
axis(2,las=2,cex.axis=1.3,at=c(seq(-1,range[2],by=1)),pos=1971,tck=tckk*2)
lines(c(1971,2020),c(0,0),col="grey50")

axis(side=1,las=1,cex.axis=1.3,at=c(seq(1980,2020,by=10)),pos=range[1],labels=c(rep("",5)),tck=tckk*2)
text(c(seq(1977.3,2020,by=20)),c(rep(range[1]-0.1*scal,5)),c(seq(1980,2020,by=20)),xpd=TRUE,cex=1.4);
lines(c(1971,2020),c(range[1],range[1]))
title(xlab="Year", line=1.4, cex.lab=1.4,xpd=TRUE)



#### females --------------
range<-range_b;scal=scal_b
par(mgp=c(1.6,0.4,0))
par(mar=c(0,1.4,0,0.3))

plot(NULL,xlim=c(1971,2020),ylim=c(range[1],range[2]),xlab="",ylab=""
     ,las=1, axes=FALSE)
rect(1971,range[2],2020,range[1],col="grey90",border=NA)
x<-10
repeat{
  lines(c(1971,2020),c(x,x),col="white",lwd=0.7)
  x=x+10;if(x>100){break}}
lines(c(1971,2020),c(50,50),col="white",lwd=1.7)


x<-1980
repeat{
  lines(c(x,x),c(range[1],range[2]),col="white",lwd=0.7)
  x=x+10;if(x>2020){break}}

# breaking points identification
xx=1
breakpo(post_melanoma_swe_fem_1y,0)
xx=xx+1
breakpo(post_melanoma_swe_fem_cond,0.6)
xx=xx+1
breakpo(post_melanoma_swe_fem_5y,-0.6)

# 95% credible interval for survival
xx=1
polyg_surv(post_melanoma_swe_fem_1y);xx=xx+1
polyg_surv(post_melanoma_swe_fem_cond);xx=xx+1
polyg_surv(post_melanoma_swe_fem_5y)

# data points of estimated survival
xx=1
points(melanoma[melanoma$shou=="swe_fem_",]$surv_1y~melanoma[melanoma$shou=="swe_fem_",]$year
       ,pch=17,col=colc[xx],cex=1);xx=xx+1
points(melanoma[melanoma$shou=="swe_fem_",]$surv_cond~melanoma[melanoma$shou=="swe_fem_",]$year
       ,pch=17,col=colc[xx],cex=1);xx=xx+1
points(melanoma[melanoma$shou=="swe_fem_",]$surv_5y~melanoma[melanoma$shou=="swe_fem_",]$year
       ,pch=17,col=colc[xx],cex=1)

tckk=-0.016
# fitted lines for survival
xx=1
surv_fit(post_melanoma_swe_fem_1y);xx=xx+1
surv_fit(post_melanoma_swe_fem_cond);xx=xx+1
surv_fit(post_melanoma_swe_fem_5y);xx=xx+1

axis(2,las=2,cex.axis=1.3,at=seq(range[1],range[2],by=10),
     labels=c(rep("",length(seq(range[1],range[2],by=10)))),pos=1971,tck=tckk)

axis(2,las=2,cex.axis=1.3,at=c(seq(range[1],range[2],by=20)),pos=1971,tck=tckk)

axis(side=1,las=1,cex.axis=1.3,at=c(seq(1980,2020,by=10)),pos=range[1],
     labels = c(rep("",5)),tck=tckk)
lines(c(1971,2020),c(range[1],range[1]))
text(1974,range[2]-0.05*scal,"h",cex=2.2)

##### slope plot --------------
range<-range_c;scal=scal_c
par(mar=c(2.5,1.4,0,0))

plot(NULL,xlim=c(1971,2020),ylim=c(range[1],range[2]),xlab="",ylab=""
     ,las=1, axes=FALSE)
rect(1971,range[2],2020,range[1],col="grey90",border=NA)
x<-range[1]
repeat{
  lines(c(1971,2020),c(x,x),col="white",lwd=0.7)
  x=x+0.5;if(x>100){break}}

x<-1980
repeat{
  lines(c(x,x),c(range[1],range[2]),col="white",lwd=0.7)
  x=x+10;if(x>2020){break}}

xx=1
polyg_slope(post_melanoma_swe_fem_1y);xx=xx+1
polyg_slope(post_melanoma_swe_fem_cond);xx=xx+1
polyg_slope(post_melanoma_swe_fem_5y)

xx=1
slope_fit(post_melanoma_swe_fem_1y);xx=xx+1
slope_fit(post_melanoma_swe_fem_cond);xx=xx+1
slope_fit(post_melanoma_swe_fem_5y)

axis(2,las=2,cex.axis=1.3,at=c(seq(range[1],range[2],by=0.5)),
     labels=c(rep("",length(seq(range[1],range[2],by=0.5)))),pos=1971,tck=tckk*2)
axis(2,las=2,cex.axis=1.3,at=c(seq(-1,range[2],by=1)),pos=1971,tck=tckk*2)
lines(c(1971,2020),c(0,0),col="grey50")

axis(side=1,las=1,cex.axis=1.3,at=c(seq(1980,2020,by=10)),pos=range[1],labels=c(rep("",5)),tck=tckk*2)
text(c(seq(1977.3,2020,by=20)),c(rep(range[1]-0.1*scal,5)),c(seq(1980,2020,by=20)),xpd=TRUE,cex=1.4);
lines(c(1971,2020),c(range[1],range[1]))
title(xlab="Year", line=1.4, cex.lab=1.4,xpd=TRUE)
