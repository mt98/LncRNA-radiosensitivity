#####Use spearman correlations in a network approach. SF2 distribution is non-linear, hence spearman will be used


module load apps/gcc/R/3.6.1 
qrsh -l short -V -cwd -pe smp.pe 8 R --vanilla --interactive

setwd("C:/Users/mqbpkmk3/Documents/fifth year")
#load libraries
library("tidyr")
library("dplyr")
library("Hmisc")
 library( reshape2 )
library( igraph )
library( survival )
library( ggplot2 )
library( gplots )
library( glmnet )
library( cluster )
library( plyr )
library(survival)
library(survminer)
library( pamr )
library( packHV )
library( fpc )

#read the annotation files 

exonarrayannotatedfiles<-get(load("cervixSF2_RMA.RData"))

#now I need to read the SF2 files
SF2calculation<-read.table("cervixSF2.csv",sep=",",header=T,row.names=1,stringsAsFactors = FALSE)
#SF2calculation<-read.table("C:/Users/mqbpkmk3/Documents/third year/cervix data/cervix info files/cervixSF2.csv",sep=",",header=T,row.names=1,stringsAsFactors = FALSE)
#develop a correlation model between the SF2 values and the lncRNA expression values 
#log transform microarray values

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

SF2values<-lapply(patientnames,function(x){
  p<-SF2calculation[SF2calculation$vno==x,]
  p$sf2})
SF2values<-as.vector(unlist(SF2values))
names(SF2values)<-patientnames
SF2values<-na.omit(SF2values)
colnames(combination)=patientnames
combineSF2expression<-merge(as.data.frame(t(combination)),SF2values,by="row.names")
colnames(combineSF2expression)[ncol(combineSF2expression)] = "SF2"

#cervix correlation based on SF2

cervixcorrelation<- sapply(colnames(combineSF2expression)[2:9967],function(x){
  cancer<-combineSF2expression[, c(x,"SF2")]
  significance <- cor.test(cancer[,x],cancer$SF2,use="everything", method=c("spearman") )
  pvalue<-as.numeric(unlist(significance[3]))
  cor<-as.numeric(unlist(significance[4]))
cor<-cbind(pvalue,cor)
})
pvalue<-cervixcorrelation[1,]
FDRcorrection<- p.adjust(pvalue, method="fdr" )
pearsoncervix<- merge(as.data.frame(t(cervixcorrelation)),FDRcorrection, by="row.names")
rownames(pearsoncervix)<-pearsoncervix$Row.names

sig_genes<-rownames(pearsoncervix[pearsoncervix$V1<0.05&pearsoncervix$V2>0,])

#survival data
SF2calculationfilter<-SF2calculation[SF2calculation$treatment_type==1,]
SF2calculationfilter1<-SF2calculationfilter[-c(3),]
rownames(SF2calculationfilter1)<-SF2calculationfilter1$vno
censorship<-60
tr_survival_event_censorship  <- ifelse(  SF2calculationfilter1$local <= censorship & SF2calculationfilter1$local_recur == 1, 1 ,0)
tr_survival_time_censorship   <- ifelse(  tr_survival_event_censorship == 0 & SF2calculationfilter1$local >= censorship, censorship , SF2calculationfilter1$local )   
SF2calculationfilter1        <- cbind( SF2calculationfilter1, tr_survival_time_censorship, tr_survival_event_censorship )
colnames( SF2calculationfilter1 )[ (ncol(SF2calculationfilter1)-1):ncol(SF2calculationfilter1) ] <- c( "censored_localtime", "censored_localstatus" )

censorship<-60
tr_survival_event_censorship  <- ifelse(  SF2calculationfilter1$MFS <= censorship & SF2calculationfilter1$distant_recur == 1, 1 ,0)
tr_survival_time_censorship   <- ifelse(  tr_survival_event_censorship == 0 & SF2calculationfilter1$MFS >= censorship, censorship , SF2calculationfilter1$MFS)   
SF2calculationfilter1        <- cbind( SF2calculationfilter1, tr_survival_time_censorship, tr_survival_event_censorship )
colnames( SF2calculationfilter1 )[ (ncol(SF2calculationfilter1)-1):ncol(SF2calculationfilter1) ] <- c( "censored_mtime", "censored_mstatus" )

censorship<-60
tr_survival_event_censorship  <- ifelse(  SF2calculationfilter1$DFS <= censorship & SF2calculationfilter1$dfs_status== 1, 1 ,0)
tr_survival_time_censorship   <- ifelse(  tr_survival_event_censorship == 0 & SF2calculationfilter1$DFS >= censorship, censorship , SF2calculationfilter1$DFS)   
SF2calculationfilter1        <- cbind( SF2calculationfilter1, tr_survival_time_censorship, tr_survival_event_censorship )
colnames( SF2calculationfilter1 )[ (ncol(SF2calculationfilter1)-1):ncol(SF2calculationfilter1) ] <- c( "censored_dfstime", "censored_dfsstatus" )

#test with those genes which are involved in the network without looking at MAD                         
                          
cervix_network_dataset<-combination[sig_genes,]

#make a co-expression network of 197 lncRNAs in SF2 cohort of patients
#### Compute the correlation among the literature genes

which_correlation <- "spearman"
correlation_threshold <- 0.5
#correlation

tr_cor          <- cor( t(cervix_network_dataset), method= which_correlation )
#### Remove the interactions with weak coexpression strength
tr_cor[upper.tri(tr_cor, diag=TRUE)] <- NA
tr_cor_list <- melt(tr_cor)
tr_cor_list <- na.omit(tr_cor_list)
rm(tr_cor)
tr_cor_list_strong <- tr_cor_list[ tr_cor_list$value >= correlation_threshold, ] 
tr_cor_list_strong$abs_value <- abs( tr_cor_list_strong$value )
#### Keep only the most connected component of the network
tr_network      <- graph.data.frame( tr_cor_list_strong, directed=F ) 
#write_graph(tr_network,"wholenetwork.graphml",format=c("graphml"))
check_isolate   <- clusters( tr_network ) 
major_component <- which( check_isolate$csize == max( check_isolate$csize ) ) 
kept_genes      <- names( check_isolate$membership[ which( check_isolate$membership %in% major_component ) ] ) 
#add the 25 genes one by one based on their correlation strength and use the number where the correlation stops

selectlncRNAs<-pearsoncervix[kept_genes,]
selectlncRNAs<-selectlncRNAs[order(-selectlncRNAs[,3]),]

combinationtotal<-cbind(as.data.frame(t(combineSF2expression[,rownames(selectlncRNAs)])))
#order the lncRNAs by rank
combinationtotal<-combineSF2expression[,c(rownames(selectlncRNAs),"SF2")]

#mean expression calculated 
performance <- NULL
for( current_size in 2:length( rownames(selectlncRNAs) ) ){
         current_signature     <- rownames( selectlncRNAs)[ 1:current_size]   
         current_score      <- apply( combinationtotal[ ,current_signature ] ,1 , mean )  
         current_cor_score <- cor( current_score, combinationtotal[ ,"SF2"], method="spearman" )
         print( c( current_size, current_cor_score) )
         performance <- c( performance, current_cor_score )
}

siggenes<-rownames(selectlncRNAs)[1:12]
#use the expression of the genes
sig_expression<-combination[siggenes,]

#use the mean expression of the genes and then put in the development of the signature
#mean expression of the genes 

meanexpressiontest<-sapply(colnames(combination),function(x){
  p<-mean(as.numeric(sig_expression[,x]))})


#develop radiosensitivity categories
combinationcoxtest<-merge(SF2calculationfilter1,meanexpressiontest,by="row.names")
#try different stratification methods
med<-quantile(combinationcoxtest$y,0.75)
category<-sapply(combinationcoxtest$y,function(x){
  if (x>=med){y="High score"}else if (x<med){y="Low score"}})
combinationcoxtest$category<-category
library("survminer")
library("survival")
res.cox <- coxph(Surv( censored_localtime, censored_localstatus) ~ category, data = combinationcoxtest)
summary(res.cox)
fit<-survfit(Surv(censored_localtime, censored_localstatus) ~ category, data = combinationcoxtest)
par(pty="s")
gplotmedian<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend.title="",linetype=c(1,3),legend.labs=c("Radioresistant","Radiosensitive"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=10,pval.coord=c(0,0.1),
                         risk.table = TRUE,palette=c("black","gray"),ylab="Local control",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",25), font.y = c("bold",25), font.tickslab = c("bold",20),
                         font.legend.labs= c( "bold",20),risk.table.fontsize = 7.0,risk.table.height = 0.3,risk.table.col = "strata",risk.table.y.text = FALSE
)
gplotmedian$plot <- gplotmedian$plot + labs(
  title    = "Cervix SF2"        
)
gplotmedian<- ggpar(
  gplotmedian,
  font.title    = c(25, "bold"))
gplotmedian$plot<-gplotmedian$plot + theme(plot.title = element_text(hjust = 0.5))
 
gplotmedian$table <- ggpar(gplotmedian$table,
                            font.title = list(size = 18))


res.cox <- coxph(Surv(censored_mtime, censored_mstatus) ~ category, data = combinationcoxtest)
summary(res.cox)
fit<-survfit(Surv(censored_mtime, censored_mstatus) ~ category, data = combinationcoxtest)
par(pty="s")
gplotmedian1<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend="none",linetype=c(1,3),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=10,pval.coord=c(0,0.1),
                         risk.table = TRUE,palette=c("black","gray"),ylab="Metastasis free survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",25), font.y = c("bold",25), font.tickslab = c("bold",20),
                         font.legend.labs= c( "bold",20),risk.table.fontsize = 7.0,risk.table.height = 0.3,risk.table.col = "strata",risk.table.y.text = FALSE)
gplotmedian1 <- ggpar(
  gplotmedian1,
  font.title    = c(25, "bold"))
gplotmedian1$plot<-gplotmedian1$plot + theme(plot.title = element_text(hjust = 0.5))
 
gplotmedian1$table <- ggpar(gplotmedian1$table,
                            font.title = list(size = 18))
res.cox <- coxph(Surv(censored_dfstime, censored_dfsstatus) ~ category, data = combinationcoxtest)
summary(res.cox)
fit<-survfit(Surv(censored_dfstime, censored_dfsstatus) ~ category, data = combinationcoxtest)
par(pty="s")
gplot2median<-ggsurvplot(fit,size=2,censor=TRUE,pval=T,legend="none",linetype=c(1,3),legend.labs=c("Radioresistant","Radiosensitive"),
                         risk.table.x.text = FALSE,tables.theme = clean_theme(),pval.size=10,pval.coord=c(0,0.1),
                         risk.table = TRUE,palette=c("black","gray"),ylab="Disease free survival",break.time.by = 20,
                         xlab="Time in months", font.x = c("bold",25), font.y = c("bold",25), font.tickslab = c("bold",20),
                         font.legend.labs= c( "bold",20),risk.table.fontsize = 7.0,risk.table.height = 0.3,risk.table.col = "strata",risk.table.y.text = FALSE)
gplot2median <- ggpar(
  gplot2median,
  font.title    = c(25, "bold"))
gplot2median$plot<-gplot2median$plot + theme(plot.title = element_text(hjust = 0.5))
 
gplot2median$table <- ggpar(gplot2median$table,
                            font.title = list(size = 18))
                            
splots<-list(gplotmedian,gplotmedian1, gplot2median)
s<-arrange_ggsurvplots(splots,ncol=2,nrow=3)
ggsave("/mnt/iusers01/cw01/mqbpkmk3/SupplementaryFigure4.pdf",plot=s,height=16,width=16)                           

scp mqbpkmk3@csf3.itservices.manchester.ac.uk:~/SupplementaryFigure4.pdf "/Users/mkhan/Documents/LncRNA paper rad res figures" 





