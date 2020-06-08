# Authors: Ken + Anis
# v0.2 Change log:
# added Anis' analysis
# v0.3 Change log:
# changed Ken's analysis to use the 'data' instead of 'Data'
# if dbsource = 'both', then assumed patient is followed for 4 years for mortality
# added log-rank test for Kaplan Meier curves
# removed age as a covariate
# v0.4 Change log:
# v0.5 Change log:
# changed color of censor cross on fig 3,4 to black
# changed legend label on all survival figures
# v0.6 Change log:
# Fixed ORs in the fitted LR models. Previous versions reported the log-odd, not odd ratio
# v0.7 Change log:
# Ben is messing everything up.
# changed survival plots to white background, fixed missing x/y axes
# unround to give exact p-values for the fitted LR results

library(pROC);library(survival);library(GGally);library(dplyr);library(ggplot2);library(grid);
library(caret);library(survey)

rm(list=ls())

#--------------------------------------------------------------***CAUTION***------------------------------------------------------------------
# Set the working directory to the source file location by doing the following steps: 1. Session 2. Set Working Directory 3. To Source File Location
# If you are planning to make changes to the code, PLEASE do it in a newer version of code.

dirname=(paste(getwd(),"/","Data", sep=""))
Data <- read.table(file = paste(dirname,"/",'Data_v1.1.csv', sep=""), header = TRUE, sep = ",",quote = "'", dec = ".",fill = TRUE)

# Section One: Data Preprocessing --------------------------------------------------------------------------

# Outcome Vectors 
OutcomeM <- c('hospital_expire_flag','icu_expire_flag', 'thirty_day_mort', 'one_year_mortality') # Mortality-related outcomes
OutcomeL <- c('los_hospital','los_icu') # LOS-related outcomes

# Exclusion Criteria
Ind1=(Data[,'age']>100)
cat('!',sum(Ind1),'patients with age>100')
Ind2=(Data[,'los_hospital']<=0) 
cat('!',sum(Ind2),'patients with los_hospital<=0')
Ind3=(Data[,'AG1']<0) 
Ind4=(Data[,'AG1']>100)
cat('!',sum(Ind3|Ind4),'patients with out of bound AG')
Ind5=(Data[,'survival']<=0) 
cat('!',sum(Ind5,na.rm = TRUE),'patients with survival<=0')
Ind5[is.na(Ind5)]=FALSE
Ind=Ind1|Ind2|Ind3|Ind4|Ind5
data=data.frame(Data[!Ind,])

# Section Two: Anis's Analyses ==============================================================================

## Descriptive Statistical Analysis
cat('*****',nrow(data),'patients included from CCU CSRU MICU SICU TSICU')
cat('Average age:',mean(data[,'age']))
cat((sum((data$gender)=='M')/nrow(data))*100,'male')
cat((sum((data$first_careunit)=='CCU')/nrow(data))*100,'CCU')
cat((sum((data$first_careunit)=='CSRU')/nrow(data))*100,'CSRU')
cat((sum((data$first_careunit)=='MICU')/nrow(data))*100,'MICU')
cat((sum((data$first_careunit)=='SICU')/nrow(data))*100,'SICU')
cat((sum((data$first_careunit)=='TSICU')/nrow(data))*100,'TSICU')
cat('RDW Characteristics:',mean(data$RDW1),'+/-',sd(data$RDW1))
cat('AG Characteristics:',mean(data$AG1),'+/-',sd(data$AG1))
summary(data)
sd(data$age)

## (RDW & AG) and outcomes
# 95% CI
outcome=c(OutcomeL,OutcomeM)
for (i in 1:(length(outcome))){
  x=data[,outcome[i]]
  print(outcome[i])
  error <- qt(0.975,df=length(x)-1)*sd(x)/sqrt(length(x))
  left <- round((mean(x)-error)*1,digits=2)
  right <- round((mean(x)+error)*1,digits=2)
  cat(round(mean(x)*1,digits=2),'(95% CI,',left,'to',right,')\n')
}
compar <- function(x,dec,i) {
  if (x>dec[i] && x<=dec[i+1]) {return(i)}} 

# Deciles
feat=c('RDW1','AG1')
y=matrix(NA,nrow=nrow(data),ncol=2)
levs=c('01','02','03','04','05','06','07','08','09','10')
for (k in c(1:2)){
  x=as.matrix(data[,feat[k]])
  deciles=unname(quantile(x, prob = seq(0, 1, length = 11), type = 5))
  # deciles=c(-Inf, deciles, Inf)
  for (j in c(1:10)){
    for (d in c(1:nrow(data))){
      
      if (j==1) {if (x[d]>=deciles[j] && x[d]<=deciles[j+1]) { y[d,k]=levs[(j)]}}
      else  {if (x[d]>deciles[j] && x[d]<=deciles[j+1]) { y[d,k]=levs[(j)]}}
      
    }
  }
}

colnames(y)=c('RDW_dec','AG_dec')
feat2=c('RDW_dec','AG_dec')
data_new=data.frame(data,y)
# OutcomeM
OutcomeM2 <- c('In-Hospital Mortality','In-ICU Mortality',
               '30-Days Mortality', 'One-Year Mortality') # Mortality-related outcomes
OutcomeL <- c('los_hospital','los_icu') # LOS-related outcomes


tiff("Figure1a.tiff",height = 8.3,width=4,units='in',res=600)
grid.newpage()
pushViewport(viewport(layout = grid.layout(length(OutcomeM), 1)))

for (out in c(1:length(OutcomeM))){
  Y=(data_new[data_new[,OutcomeM[out]]==1,])
  p=ggplot(Y, aes(RDW_dec))+ geom_bar()+ labs(title=OutcomeM2[out],x='RDW Deciles',y='Number of Patients')+ 
    theme_bw()
  print(p, vp = viewport(layout.pos.row = out, layout.pos.col = 1))
}
dev.off()

tiff("Figure1b.tiff",height = 8.3,width=4,units='in',res=600)
grid.newpage()
pushViewport(viewport(layout = grid.layout(length(OutcomeM), 1)))

for (out in c(1:length(OutcomeM))){
  Y=(data_new[data_new[,OutcomeM[out]]==1,])
  p=ggplot(Y, aes(AG_dec))+ geom_bar()+ labs(title=OutcomeM2[out],x='AG Deciles',y='Number of Patients')+
    theme_bw()
  print(p, vp = viewport(layout.pos.row = out, layout.pos.col = 1))
}

dev.off()
# Effect on discrimination and reclassification
# for (i in c(1:length(OutcomeM))){
i=4
O <- factor(data[,OutcomeM[i]], levels=0:1, labels=c('N', 'Y'))
Stats5 = function(...) c(twoClassSummary(...),defaultSummary(...))
train_control<- trainControl(method="repeatedcv", number=10,repeats=20,classProbs = TRUE)
# train_control$sampling <- "up"
train_control$summaryFunction <- Stats5
mydata<-cbind(data[,11],data$"AG1",data$"RDW1")
mydata<-data.frame(mydata,O)
colnames(mydata)=c('saps','ag','rdw','O')
set.seed(110)
model1<- train(O~saps, data=mydata, trControl=train_control, method="glm" ,metric='ROC')
set.seed(110)
model2<- train(O~saps+ag, data=mydata, trControl=train_control, method="glm" ,metric='ROC')
set.seed(110)
model3<- train(O~saps+rdw, data=mydata, trControl=train_control, method="glm" ,metric='ROC')
set.seed(110)
model4<- train(O~., data=mydata, trControl=train_control, method="glm" ,metric='ROC')
resamps <- resamples(list(SAPSii = model1,
                          SAPSii_AG = model2,
                          SAPSii_RDW = model3,
                          SAPSS_RDW_AG=model4))
# save(list = resamps, file = paste("Outcome",i,".RData",sep=""), envir = .GlobalEnv)
# bwplot(resamps, layout = c(3, 1))
# dotplot(resamps, metric = "ROC",auto.key=list(space="top", column=4, cex=.8, title="Sample size", 
                                               # cex.title=1, lines=TRUE, points=FALSE))
# }

dfrm <- data.frame(Input=gl(4, 3, 3*4, labels =(c('SAPSII','SAPSII+AG','SAPSII+RDW','SAPSSII+RDW+AG'))),
                   ROC=runif(12)*0+0.7,min=as.factor(runif(12)*0)
)
M1=t.test(resamps$values$`SAPSii~ROC`)
dfrm[c(1:3),2]=c(M1$conf.int,M1$estimate)
M2=t.test(resamps$values$`SAPSii_AG~ROC`)
dfrm[c(4:6),2]=c(M2$conf.int,M2$estimate)
M3=t.test(resamps$values$`SAPSii_RDW~ROC`)
dfrm[c(7:9),2]=c(M3$conf.int,M3$estimate)
M4=t.test(resamps$values$`SAPSS_RDW_AG~ROC`)
dfrm[c(10:12),2]=c(M4$conf.int,M4$estimate)
scale_color_brewer(type = 'seq', palette = 1)
tiff("Figure2_1_year.tiff",height = 4,width=7,units='in',res=600)
qplot(x=ROC, y=Input,data=dfrm,color=Input,shape=min)+
  # facet_grid(algo~day,scales = "free", space = "free")+
  geom_line(size=0.5)+
  geom_point(size=2)+
  scale_shape_manual(values=c(124,124,124))+
  theme(panel.grid.minor = element_blank())+
  theme_set(theme_gray(base_size = 9))+
  # scale_color_manual(values=c( '#000000','#000000','#000000'))+
  xlab("AUC\nConfidence Level:0.95")+
  ylab('Input')+theme_bw()
dev.off()
# Section Three: Ken's Analyses ############################################################################
Outcome <- c('hospital_expire_flag','icu_expire_flag', 'thirty_day_mort', 'one_year_mortality')

# Loop thru each outcome type and fit LR models with the covariates below
i = 1
for (i in 1:length(Outcome))
{
  
  col <- which(colnames(data) == Outcome[i])
  
  
  cat("Model of", colnames(data)[col], "\n")
  Model <- glm(as.factor(data[,col]) ~ data$AG1 
               + data$RDW1 
               + data$gender
               + data$sapsii
               + data$first_careunit
               , family = binomial())
  
  prob <- predict(Model,type=c("response"))
  data$prob=prob
  
  ROC <- roc(data[,col] ~ prob, data = data)
  
  plot(ROC)  
  
  #print(summary(Model))
  summary <- cbind(exp(Model$coefficients), summary(Model)$coefficients[,4])
  colnames(summary) <- c('OR', 'p-value')
  print(summary)
  print(auc(ROC))
  print("******************************")
  
  write.csv(summary, file = paste0(Outcome[i],"vsAG+RDW", ".csv"))
  
  i <- i+1
  
}


# RDW ONLY
print("RDW ONLY")

i = 1
for (i in 1:length(Outcome))
{
  
  col <- which(colnames(data) == Outcome[i])
  
  
  cat("Model of", colnames(data)[col], "\n")
  Model <- glm(data[,col] ~ data$RDW1 
               + data$gender
               + data$sapsii
               + data$first_careunit
               , family = binomial())
  
  prob <- predict(Model,type=c("response"))
  data$prob=prob
  
  ROC <- roc(data[,col] ~ prob, data = data)
  
  plot(ROC)  
  
  #print(summary(Model))
  summary <- cbind(exp(Model$coefficients), summary(Model)$coefficients[,4])
  colnames(summary) <- c('OR', 'p-value')
  print(summary)
  print(auc(ROC))
  print("******************************")
  
  write.csv(summary, file = paste0(Outcome[i],"vsRDW", ".csv"))
  
  i <- i+1
  
}




# ANION GAP ONLY
print("AG ONLY")

i = 1
for (i in 1:length(Outcome))
{
  
  col = which(colnames(data) == Outcome[i])
  
  
  cat("Model of", colnames(data)[col], "\n")
  Model = glm(data[,col] ~ data$AG1
              + data$gender
              + data$sapsii
              + data$first_careunit
              , family = binomial())
  
  prob=predict(Model,type=c("response"))
  data$prob=prob
  
  ROC <- roc(data[,col] ~ prob, data = data)
  
  plot(ROC)  
  
  #print(summary(Model))
  summary <- cbind(exp(Model$coefficients), summary(Model)$coefficients[,4])
  colnames(summary) <- c('OR', 'p-value')
  print(summary)
  print(auc(ROC))
  print("******************************")
  
  write.csv(summary, file = paste0(Outcome[i],"vsAG", ".csv"))
  
  i <- i+1
  
}

#SAPS-II ONlY
print("SAPS-2 ONLY")

i = 1
for (i in 1:length(Outcome))
{
  
  col <- which(colnames(data) == Outcome[i])
  
  
  cat("Model of", colnames(data)[col], "\n")
  Model <- glm(data[,col] ~ data$gender
               + data$sapsii
               + data$first_careunit
               , family=binomial())
  
  prob <- predict(Model,type=c("response"))
  data$prob=prob
  
  ROC <- roc(data[,col] ~ prob, data = data)
  
  plot(ROC)  
  
  #print(summary(Model))
  summary <- cbind(exp(Model$coefficients), summary(Model)$coefficients[,4])
  colnames(summary) <- c('OR', 'p-value')
  print(summary)
  print(auc(ROC))
  print("******************************")
  
  write.csv(summary, file = paste0(Outcome[i],"vsSAPS2only", ".csv"))
  
  i <- i+1
  
}



# RDW quartiles
data$RDW1Quart <- ntile(data$RDW1, 4) 

# AG quartiles
data$AG1Quart <- ntile(data$AG1, 4)


# fill in the NA with censored survival lengths
i <- 1
for (i in 1:nrow(data))
{
  if (is.na(data$survival[i]) &  data$dbsource[i] == 'metavision' )
  {
    data$survival[i] <- 90
  }
  
  if (is.na(data$survival[i]) & (data$dbsource[i] == 'carevue' | data$dbsource[i] == 'both'))
  {
    data$survival[i] <- 1460
  }
  
  i <- i +1
  
}

sum(is.na(data$survival) == FALSE)
sum(is.na(data$survival)) # Should be zero
sum(data$survival > 30)
## WARNING OF REMOVED N ROWS CONTAINING MISSING VALUES ARE NORMAL BECAUSE SOME PATIENTS ARE FOLLOWED FOR A VERY LONG TIME

# Color function
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
######


# 30 Days Mortality
surv.obj <- survfit(Surv(time = survival, event = thirty_day_mort) ~ RDW1Quart, data=data)
Fig1 <- ggsurv(surv.obj, CI=T) + 
  xlab("Time (days)") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12)) + 
  scale_colour_manual(name="RDW Quartiles", 
                      values = gg_color_hue(4)) +
  guides(linetype = F) +
  scale_x_continuous(limits = c(0,30))
  
Fig1
survdiff(Surv(time = survival, event = thirty_day_mort) ~ RDW1Quart, data=data)



surv.obj <- survfit(Surv(time = survival, event = thirty_day_mort) ~ AG1Quart, data=data)
Fig2 <- ggsurv(surv.obj, CI=T) + 
  xlab("Time (days)") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12)) + 
  scale_colour_manual(name="AG Quartiles", 
                      values = gg_color_hue(4)) +
  guides(linetype = F) +
  scale_x_continuous(limits = c(0,30))
Fig2
survdiff(Surv(time = survival, event = thirty_day_mort) ~ AG1Quart, data=data)


# 1YR Mortality
surv.obj <- survfit(Surv(time = survival, event = one_year_mortality) ~ RDW1Quart, data=data)
Fig3 <- ggsurv(surv.obj, CI=T, cens.col = 1) + 
  xlab("Time (days)") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12)) +
  scale_colour_manual(name="RDW Quartiles", 
                      values = gg_color_hue(4)) +
  guides(linetype = F) +
  scale_x_continuous(limits = c(0,365))
Fig3
survdiff(Surv(time = survival, event = one_year_mortality) ~ RDW1Quart, data=data)




surv.obj <- survfit(Surv(time = survival, event = one_year_mortality) ~ AG1Quart, data=data)
Fig4 <- ggsurv(surv.obj, CI=T,cens.col = 1) + 
  xlab("Time (days)") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=12), 
        legend.text=element_text(size=12)) +
  scale_colour_manual(name="AG Quartiles", 
                      values = gg_color_hue(4)) +
  guides(linetype = F) +
  scale_x_continuous(limits = c(0,365))
Fig4
survdiff(Surv(time = survival, event = one_year_mortality) ~ AG1Quart, data=data)

# Save via Export
grid_plot <- gridExtra::grid.arrange(Fig1, Fig2)
grid_plot <- gridExtra::grid.arrange(Fig3, Fig4)








