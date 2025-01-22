setwd("/Users/mkhan/Documents/fifth year")
# Read in the survival patient data for 112 patients ------------------------------------------------------
survivalpatient<-read.table("/Users/mkhan/Documents/ third year/cervix data/cervix info files/cervix 112 new.csv",sep=",",quote="",header=T,stringsAsFactors = FALSE)
#read in the lncRNA exon array data for the patients
lncRNAsurvivalexonarrayannotatedfiles<-get(load("cervixsurvival_RMA.RData"))
#remove in files of patients which had sequencing from two seperate FFPE blocks
##V618 and V362
#change the patient file names to match the survival data file name
FilterlncRNAsurvival<-lncRNAsurvivalexonarrayannotatedfiles[,-c(49,108,1:11)]
rownames(FilterlncRNAsurvival)<-lncRNAsurvivalexonarrayannotatedfiles$NONCODE3.ID
nameslncrnasurvival<-sapply(colnames(FilterlncRNAsurvival),function(x){
  NAME<-strsplit(x, fixed=T, "_")[[1]][[5]]
  PATIENT<-strsplit(NAME, fixed=T, ".")[[1]][[1]]}
)
nameslncrnasurvival<-toupper(nameslncrnasurvival)
colnames(FilterlncRNAsurvival)<-nameslncrnasurvival

rownames(survivalpatient)<-survivalpatient$vno
nolocal<-rownames(survivalpatient[survivalpatient$local_recur==0,])
local<-rownames(survivalpatient[survivalpatient$local_recur==1,])
#use limma model analysis to do a differential expression between patients with local recurrence and no local recurrence
target<-c(local,nolocal)
u<-FilterlncRNAsurvival[,match(target, colnames(FilterlncRNAsurvival))]
library("limma")
Group <- rep(c("local","Nolocal"),c(41,71))
design.matrix <- model.matrix(~ 0+~ factor(Group)) 
colnames(design.matrix) <- c("local","nolocal")
cont_matrix <- makeContrasts(localvsnolocal = local-nolocal, levels=design.matrix)
fit <- lmFit(as.matrix(u), design.matrix)
fit_contrast <- contrasts.fit(fit, cont_matrix)
fit_contrast <- eBayes(fit_contrast)
display<-topTable(fit_contrast,adjust="BH",n=11534)
###MDS plot showing seperation
df2 <- data.frame(apply(T24hypoxia0.2arrayAdeltaCt, 2, function(x) as.numeric(as.character(x))))
rownames(df2)<-rownames(T24hypoxia0.2arrayAdeltaCt)
colnames(df2)<-colnames(T24hypoxia0.2arrayAdeltaCt)
hypoxia.pca <- prcomp(t(df2),center = TRUE, scale=TRUE)
library(ggbiplot)
condition<-c("normoxia","normoxia","normoxia","normoxia","normoxia","normoxia","hypoxia","hypoxia","hypoxia","hypoxia","hypoxia","hypoxia")
ggbiplot(hypoxia.pca,obs.scale = 1, var.scale = 1,var.axes=FALSE,groups=condition,labels=colnames(df2) )+
  scale_colour_manual(name="Origin", values= c("forest green", "red3", "dark blue"))+
  ggtitle(" T24 hypoxia 0.2 percent and normoxia array A")+
  theme_minimal()+
  theme(legend.position = "bottom")+xlim(-30,30)







gene<-c("n340217","n346223,n385225,n387638","n337764","n326470,n326468,n326467,n326466,n326469","n342607","n339076","n339447")
#heatmap for the seven differentially expressed genes
u2<-u[gene,]
u2<-t(scale(t(u2)))
library(ComplexHeatmap)
library(circlize)
df = data.frame(type = c(rep("No local recurrence", 71), rep("Local recurrence", 41)))
ha = HeatmapAnnotation(df = df,c("black","white"))
Heatmap(u2,cluster_rows = function(m) hclust(dist(m),method="complete"),
        cluster_columns = function(m) hclust(dist(m),method="complete"), name = "Expression",column_title = "Cervix carcinoma",  
        column_title_side = "bottom",bottom_annotation = ha)
#Kaplan Maier for training on the seven genes
gene5<-c("n340217","n346223,n385225,n387638","n337764","n326470,n326468,n326467,n326466,n326469","n342607","n339076","n339447")
survival5<-FilterlncRNAsurvival[gene5,]
#censor the survival dataset
censorship<-60
tr_survival_event_censorship  <- ifelse( survivalpatient$local <= censorship & survivalpatient$local_recur == 1, 1 ,0)
tr_survival_time_censorship   <- ifelse(  tr_survival_event_censorship == 0 & survivalpatient$local>= censorship, censorship , survivalpatient$local)   
survivalpatient     <- cbind( survivalpatient, tr_survival_time_censorship, tr_survival_event_censorship )
colnames( survivalpatient )[ (ncol(survivalpatient)-1):ncol(survivalpatient ) ] <- c( "censored_lpfstime", "censored_lpfsstatus" )
censorship<-60
tr_survival_event_censorship  <- ifelse(survivalpatient$DFS <= censorship & survivalpatient$DFS.status == 1, 1 ,0)
tr_survival_time_censorship   <- ifelse(  tr_survival_event_censorship == 0 & survivalpatient$DFS>= censorship, censorship ,survivalpatient$DFS)   
survivalpatient      <- cbind( survivalpatient, tr_survival_time_censorship, tr_survival_event_censorship )
colnames( survivalpatient )[ (ncol(survivalpatient)-1):ncol(survivalpatient ) ] <- c( "censored_dfstime", "censored_dfsstatus" )
censorship<-60
tr_survival_event_censorship  <- ifelse(survivalpatient$MFS <= censorship & survivalpatient$distant_recur == 1, 1 ,0)
tr_survival_time_censorship   <- ifelse(  tr_survival_event_censorship == 0 & survivalpatient$MFS>= censorship, censorship ,survivalpatient$MFS)   
survivalpatient      <- cbind( survivalpatient, tr_survival_time_censorship, tr_survival_event_censorship )
colnames( survivalpatient )[ (ncol(survivalpatient)-1):ncol(survivalpatient ) ] <- c( "censored_mfstime", "censored_mfsstatus" )
#training Kaplan Maier
###n339076 was the 
genenegative<-c("n340217","n346223,n385225,n387638","n337764")
genepositive<-c("n326470,n326468,n326467,n326466,n326469",
                "n342607","n339076")
meanexpressiontestp<-sapply(colnames(survival5),function(x){
  u<-sum(as.numeric(survival5[genepositive,x]))
  p<-sum(as.numeric(survival5[genenegative,x]))
  p<-mean(u-p)})
rownames(survivalpatient)<-survivalpatient$vno
combinationcox<-merge(survivalpatient,meanexpressiontestp,by="row.names")
colnames(combinationcox)[ncol(combinationcox)]<-"meanexpressiontestp"
med<-median(combinationcox$meanexpressiontestp)
category<-sapply(combinationcox$meanexpressiontestp,function(x){
  if (x>=med){y="High score"}else if (x<med){y="Low score"}})
combinationcox$category<-category
library("survival")
library("survminer")
res.coxsurvival <- coxph(Surv(censored_lpfstime, censored_lpfsstatus) ~ category, data = combinationcox)
summary(res.coxsurvival)
#combining all differentially expressed genes with the seven genes using summation..
fit<-survfit(Surv(censored_lpfstime, censored_lpfsstatus) ~ category, data = combinationcox)
par(pty="s")
gplot5<-ggsurvplot(fit,size=2,censor=TRUE,pval="P<0.0001",legend.title="",linetype=c(1,3),legend.labs=c("Radioresistant","Radiosensitive"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=10,pval.coord=c(0,0.1),
                         risk.table = TRUE,palette=c("red","blue"),ylab="Local control",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",25), font.y = c("bold",25), font.tickslab = c("bold",20),
                         font.legend.labs= c(20,"bold"),risk.table.fontsize = 7.0,risk.table.height = 0.3,risk.table.col = "strata",risk.table.y.text = FALSE)
gplot5$plot <- gplot5$plot + labs(
  title    = "Cervix training"        
)
gplot5 <- ggpar(
  gplot5,
  font.title    = c(25, "bold"))
gplot5$plot<-gplot5$plot + theme(plot.title = element_text(hjust = 0.5))
 
gplot5$table <- ggpar(gplot5$table,
                            font.title = list(size = 16))


res.coxsurvival<- coxph(Surv(censored_mfstime, censored_mfsstatus) ~ Age.of.FFPE.block, data = combinationcox)
summary(res.coxsurvival)
res.coxsurvival<- coxph(Surv(censored_mfstime, censored_mfsstatus) ~ category, data = combinationcox)
summary(res.coxsurvival)
fit<-survfit(Surv(censored_mfstime, censored_mfsstatus) ~ category, data = combinationcox)
par(pty="s")
gplot5m<-ggsurvplot(fit,size=2,censor=TRUE,pval="P=0.2",legend.title="",linetype=c(1,3),legend.labs=c("Radioresistant","Radiosensitive"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=10,pval.coord=c(0,0.1),
                         risk.table = TRUE,palette=c("red","blue"),ylab="Metastasis free survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",25), font.y = c("bold",25), font.tickslab = c("bold",20),
                         font.legend.labs= c(20,"bold"),risk.table.fontsize = 7.0,risk.table.height = 0.3,risk.table.col = "strata",risk.table.y.text = FALSE)
gplot5m <- ggpar(
  gplot5m,
  font.title    = c(25, "bold"))
gplot5m$plot<-gplot5m$plot + theme(plot.title = element_text(hjust = 0.5))
 
gplot5m$table <- ggpar(gplot5m$table,
                            font.title = list(size = 16))

res.coxsurvival<- coxph(Surv(censored_dfstime, censored_dfsstatus ) ~ Age.of.FFPE.block, data = combinationcox)
summary(res.coxsurvival)
res.coxsurvival<- coxph(Surv(censored_dfstime, censored_dfsstatus ) ~ category, data = combinationcox)
summary(res.coxsurvival)
fit<-survfit(Surv(censored_dfstime, censored_dfsstatus ) ~ category, data = combinationcox)
par(pty="s")
gplotmd<-ggsurvplot(fit,size=2,censor=TRUE,pval="P=0.003",legend.title="",linetype=c(1,3),legend.labs=c("Radioresistant","Radiosensitive"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=10,pval.coord=c(0,0.1),
                         risk.table = TRUE,palette=c("red","blue"),ylab="Disease free survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",25), font.y = c("bold",25), font.tickslab = c("bold",20),
                         font.legend.labs= c(20,"bold"),risk.table.fontsize = 7.0,risk.table.height = 0.3,risk.table.col = "strata",risk.table.y.text = FALSE)
gplotmd <- ggpar(
  gplotmd,
  font.title    = c(25, "bold"))
gplotmd$plot<-gplotmd$plot + theme(plot.title = element_text(hjust = 0.5))
 
gplotmd$table <- ggpar(gplotmd$table,
                            font.title = list(size = 16))

library("gridExtra")
library("ggplot2")
splots<-list(gplot5,gplot5m,gplotmd)
s<-arrange_ggsurvplots(splots,ncol=1,nrow=3)
ggsave("cervixoutcomemean7genes.pdf",plot=s,height=12,width=12)

#Look at the validation dataset
#read the annotation files 
exonarrayannotatedfiles<-get(load("cervixSF2_RMA.RData"))
combination<-(exonarrayannotatedfiles[,12:59])
rownames(combination)<-exonarrayannotatedfiles$NONCODE3.ID
# what i need to do is to associate each of the values in the patient sample of the combination with an SF2 value and create a vector of it
patientnames<-sapply(colnames(combination)[1:47],function(x){
 NAME<-strsplit(x, fixed=T, "_")[[1]][[5]]
  PATIENT<-strsplit(NAME, fixed=T, ".")[[1]][[1]]}
)
patientnames<-toupper(patientnames)
patientnames_v<-sapply(colnames(combination)[48],function(x){
  PATIENT<-strsplit(x, fixed=T, ".")[[1]][[1]]}
)
patientnames<-c(patientnames,patientnames_v)
colnames(combination)<-patientnames
#now I need to read the SF2 survival files
SF2calculation<-read.table("/Users/mkhan/Documents/ third year/cervix data/cervix info files/cervixSF2.csv",sep=",",header=T,row.names=1,stringsAsFactors = FALSE)
SF2calculationfilter<-SF2calculation[SF2calculation$treatment_type==1,]
rownames(SF2calculationfilter)<-SF2calculationfilter$vno
SF2calculationfilter1<-SF2calculationfilter[-c(3),]
rownames(SF2calculationfilter1)<-SF2calculationfilter1$vno
combination5<-combination[,rownames(SF2calculationfilter1)]
#look at the genes as median
expressionsf2<-sapply(colnames(combination5),function(x){
  u<-sum(as.numeric(combination5[genepositive,x]))
  p<-sum(as.numeric(combination5[genenegative,x]))
  p<-mean(u-p)})
#censor the sf2 data
censorship<-60
tr_survival_event_censorship  <- ifelse( SF2calculationfilter1$local <= censorship & SF2calculationfilter1$local_recur == 1, 1 ,0)
tr_survival_time_censorship   <- ifelse(  tr_survival_event_censorship == 0 & SF2calculationfilter1$local>= censorship, censorship , SF2calculationfilter1$local)   
SF2calculationfilter1   <- cbind( SF2calculationfilter1, tr_survival_time_censorship, tr_survival_event_censorship )
colnames( SF2calculationfilter1 )[ (ncol(SF2calculationfilter1)-1):ncol(SF2calculationfilter1 ) ] <- c( "censored_lpfstime", "censored_lpfsstatus" )
censorship<-60
tr_survival_event_censorship  <- ifelse(SF2calculationfilter1$DFS <= censorship & SF2calculationfilter1$dfs_status == 1, 1 ,0)
tr_survival_time_censorship   <- ifelse(  tr_survival_event_censorship == 0 & SF2calculationfilter1$DFS>= censorship, censorship ,SF2calculationfilter1$DFS)   
SF2calculationfilter1   <- cbind( SF2calculationfilter1, tr_survival_time_censorship, tr_survival_event_censorship )
colnames( SF2calculationfilter1)[ (ncol(SF2calculationfilter1)-1):ncol(SF2calculationfilter1) ] <- c( "censored_dfstime", "censored_dfsstatus" )
censorship<-60
tr_survival_event_censorship  <- ifelse(SF2calculationfilter1$MFS <= censorship & SF2calculationfilter1$distant_recur == 1, 1 ,0)
tr_survival_time_censorship   <- ifelse(  tr_survival_event_censorship == 0 & SF2calculationfilter1$MFS>= censorship, censorship ,SF2calculationfilter1$MFS)   
SF2calculationfilter1  <- cbind( SF2calculationfilter1, tr_survival_time_censorship, tr_survival_event_censorship )
colnames( SF2calculationfilter1 )[ (ncol(SF2calculationfilter1)-1):ncol(SF2calculationfilter1 ) ] <- c( "censored_mfstime", "censored_mfsstatus" )
combinationcox<-merge(SF2calculationfilter1,expressionsf2,by="row.names")
colnames(combinationcox)[ncol(combinationcox)]<-"expressionsf2"
med<-median(combinationcox$expressionsf2)
category<-sapply(combinationcox$expressionsf2,function(x){
  if (x>=med){y="High score"}else if (x<med){y="Low score"}})
combinationcox$category<-category
library("survival")
library("survminer")
res.coxsurvival <- coxph(Surv(censored_lpfstime, censored_lpfsstatus) ~ category, data = combinationcox)
summary(res.coxsurvival)
fit<-survfit(Surv(censored_lpfstime, censored_lpfsstatus) ~ category, data = combinationcox)
par(pty="s")
gplot5all<-ggsurvplot(fit,size=2,censor=TRUE,pval="P=0.74",legend.title="",linetype=c(1,3),legend.labs=c("Radioresistant","Radiosensitive"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=10,pval.coord=c(0,0.1),
                         risk.table = TRUE,palette=c("red","blue"),ylab="Local control",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",25), font.y = c("bold",25), font.tickslab = c("bold",20),
                         font.legend.labs= c(20,"bold"),risk.table.fontsize = 7.0,risk.table.height = 0.3,risk.table.col = "strata",risk.table.y.text = FALSE)
gplot5all$plot <- gplot5all$plot + labs(
  title    = "Cervix validation"        
)
gplot5all <- ggpar(
  gplot5all,
  font.title    = c(25, "bold"))
gplot5all$plot<-gplot5all$plot + theme(plot.title = element_text(hjust = 0.5))
 
gplot5all$table <- ggpar(gplot5all$table,
                            font.title = list(size = 16))


res.coxsurvival<- coxph(Surv(censored_mfstime, censored_mfsstatus) ~ category, data = combinationcox)
summary(res.coxsurvival)
fit<-survfit(Surv(censored_mfstime, censored_mfsstatus) ~ category, data = combinationcox)
par(pty="s")
gplot5mall<-ggsurvplot(fit,size=2,censor=TRUE,pval="P=0.23",legend.title="",linetype=c(1,3),legend.labs=c("Radioresistant","Radiosensitive"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=10,pval.coord=c(0,0.1),
                         risk.table = TRUE,palette=c("red","blue"),ylab="Metastasis free survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",25), font.y = c("bold",25), font.tickslab = c("bold",20),
                         font.legend.labs= c( 20,"bold"),risk.table.fontsize = 7.0,risk.table.height = 0.3,risk.table.col = "strata",risk.table.y.text = FALSE)
gplot5mall <- ggpar(
  gplot5mall,
  font.title    = c(25, "bold"))
gplot5mall$plot<-gplot5mall$plot + theme(plot.title = element_text(hjust = 0.5))
 
gplot5mall$table <- ggpar(gplot5mall$table,
                            font.title = list(size = 16))

res.coxsurvival<- coxph(Surv(censored_dfstime, censored_dfsstatus) ~ category, data = combinationcox)
summary(res.coxsurvival)
fit<-survfit(Surv(censored_dfstime, censored_dfsstatus) ~ category, data = combinationcox)
par(pty="s")
gplotmdall<-ggsurvplot(fit,size=2,censor=TRUE,pval="P=0.19",legend.title="",linetype=c(1,3),legend.labs=c("Radioresistant","Radiosensitive"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=10,pval.coord=c(0,0.1),
                         risk.table = TRUE,palette=c("red","blue"),ylab="Disease free survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",25), font.y = c("bold",25), font.tickslab = c("bold",20),
                         font.legend.labs= c(20,"bold"),risk.table.fontsize = 7.0,risk.table.height = 0.3,risk.table.col = "strata",risk.table.y.text = FALSE)
gplotmdall <- ggpar(
  gplotmdall,
  font.title    = c(25, "bold"))
gplotmdall$plot<-gplotmdall$plot + theme(plot.title = element_text(hjust = 0.5))
 
gplotmdall$table <- ggpar(gplotmdall$table,
                            font.title = list(size = 16))


splots<-list(gplot5,gplot5m,gplotmd,gplot5all,gplot5mall,gplotmdall)
s<-arrange_ggsurvplots(splots,ncol=2,nrow=3)
ggsave("cervixSF27genesignaturescoretrainingvalidation.pdf",plot=s,height=18,width=18)
