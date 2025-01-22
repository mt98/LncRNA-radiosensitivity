#The purpose of this code is to calculate correlation between the radiosensitivity measurements (AUC) and the lncRNAs in a pan cancer setting
setwd("/Users/mkhan/Documents/ third year")
#load the AUC radiosensitivity measurements for the CCLE cell lines
pancancer<-read.csv(file="/Users/mkhan/Documents/ third year/CCLE data/info files//anotated.csv", header=T, sep=",",stringsAsFactors = FALSE)
cellines<-pancancer$Name
#load the R object for the RNA-seq expression data for the CCLE
load("CCLE.RData")
#In order to identify the lncRNAs from the CCLE object, the names from the TCGA lncRNAs were identified
files <- list.files(path="/Users/mkhan/Documents/ third year/TCGA data", pattern="*rnaexpr.tsv",full.names=T, recursive=TRUE)
Data1 <- lapply(files
                ,function(x){
                  t <- read.table(x, header=T, sep="\t", row.names=1)
                }
)
Cancers <-lapply(Data1
                 ,function(x){
                   Patient <-lapply(colnames(x), function(y){strsplit(y, fixed=T, ".")[[1]][[2]]})
                   Patient <- unlist(Patient)
                   Tumour <- Patient == "Tumor"
                   Cancer <- x[, Tumour]
                 }  
)
CombinedCancers<-do.call('cbind', Cancers)
#acquire the ENSEMBL ids of lncRNAs from the TCGA
lncRNAs<-rownames(CombinedCancers)
genes<-rownames(Data)
#identify the lncRNAs which are present in CCLE
intersection<-intersect(lncRNAs,genes)
increasedvector<-sapply(cellines, function(x){
  p<-which(grepl(x, colnames(Data)))
  })
  #only get the lncRNAs from the CCLE data
allCancers<-Data[intersection,unlist(increasedvector)]
#log transform the CCLE data
allCancers<-log2(allCancers+1)
#filter out the lncRNAs which are zero across all the cell lines
allCancers<- allCancers[rowSums(allCancers) > 0, ]#12432
#from the AUC measurement file only take the AUC value column
pancancerradiation<-sapply(cellines, function(x)
{
  measurement<-pancancer[pancancer$Name==x,"AUC"]
  
  
})
names(pancancerradiation)=cellines
#merge the gene expression data of the cell lines and their AUC values
combine <- merge(t(allCancers),pancancerradiation,by="row.names")
colnames(combine)[ncol(combine)] = "AUC"
#calculate correlation between AUC and expression values 
pancancercorrelation<- sapply( colnames(combine)[2:12433],function(i){
  cancer<-combine[, c(i,"AUC")]
  significance <- cor.test(cancer[,i],cancer$AUC, method=c("pearson") )
 cor<-significance[c(3,4)]
})
pancancercorrelation<-t(pancancercorrelation)
pvalue<-pancancercorrelation[,1]
FDRcorrection<- p.adjust(pvalue, method="fdr" )
pearson<- merge(as.data.frame(pancancercorrelation),FDRcorrection, by="row.names")
#plot for the lncRNAs which are significant based on pvalue< 0.05
write.table(as.matrix(pearson),file="pancancercorrelation.csv",quote=F,sep=",")
# Pan-cancer correlations based on  using SF2 measurements from Torres paper rather than the AUC measurements---------------------------------------
torres<-read.table(file="/Users/mkhan/Documents/ third year/CCLE data/info files//SF2measurementsNCI-60paneltorres.csv", header=T, sep=",",stringsAsFactors=F)
rownames(torres)<-torres$Cancer
#reduce the CCLE data set to the cell lines which have SF2 measurements according to the Torres paper
increasedvectort<-sapply(torres$Cancer, function(x){
  p<-which(grepl(x, colnames(Data)))
})
Cancerst<-Data[intersection,unlist(increasedvectort)]

cancerradiation<-sapply(torres$Cancer, function(x)
{
  measurement<-torres[torres$Cancer==x,"SF2"]
  
  
})
names(cancerradiation)=torres$Cancer
Cancerst<-t(Cancerst)
combinet <- merge(as.data.frame(Cancerst),cancerradiation,by="row.names")
colnames(combinet)[ncol(combinet)] = "SF2"
#pan cancer correlations filter out all zero
pancancercorrelationt<- sapply(colnames(combinet)[2:12557],function(x){
  cancer<-combinet[, c(x,"SF2")]
  significance <- cor.test(cancer[,x],cancer$SF2,use="everything", method=c("pearson") )
  cor<-significance[c(3,4)]
})
pancancercorrelationt<-t(pancancercorrelationt)
pvaluet<-pancancercorrelationt[,1]
FDRcorrectiont<- p.adjust(pvaluet, method="fdr" )
pearsont<- merge(as.data.frame(pancancercorrelationt),FDRcorrectiont, by="row.names")
pearsont<-as.matrix(pearsont)#need to unlist the variables in pearson
#write.table(pearsont,file="pearsontorres.csv",sep=",")

#Pan cancer correlations based on SF2 measurements from Alberts paper
albertsf2<-read.table(file="/Users/mkhan/Documents/ third year/CCLE data/info files/SF2 measurements alberts.csv", header=T, sep=",",stringsAsFactors=F)
rownames(albertsf2)<-albertsf2$Cancer
increasedvectora<-sapply(albertsf2$Cancer, function(x){
  p<-which(grepl(x, colnames(Data)))
})
Cancersa<-Data[intersection,unlist(increasedvectora)]
cancerradiationa<-sapply(albertsf2$Cancer, function(x)
{
  measurement<-albertsf2[albertsf2$Cancer==x,"SF2"]
  
  
})
names(cancerradiationa)=albertsf2$Cancer
Cancersa<-t(Cancersa)
combinea <- merge(as.data.frame(Cancersa),cancerradiationa,by="row.names")
colnames(combinea)[ncol(combinea)] = "SF2"
pancancercorrelationa<- sapply(colnames(combinea)[2:12557],function(x){
  cancer<-combinea[, c(x,"SF2")]
  significance <- cor.test(cancer[,x],cancer$SF2,use="everything", method=c("pearson") )
  cor<-significance[c(3,4)]
})
pancancercorrelationa<-t(pancancercorrelationa)
pvaluea<-pancancercorrelationa[,1]
FDRcorrectiona<- p.adjust(pvaluea, method="fdr" )
pearsona<- merge(as.data.frame(pancancercorrelationa),FDRcorrectiona, by="row.names")
#write.table(pearsona,file="pearsonalberts.csv",sep=",")
# To see the correlations of the SF2 of Torres paper with AUC metric
AUC<-read.table(file="/Users/mkhan/Documents/ third year/CCLE data/info files/AUC measurements.csv",sep=",",header=T,stringsAsFactors=F)
rownames(AUC)<-AUC$Name
intersectionaucsf2<-intersect(torres$Cancer, AUC$Name)
cancerradiationin<-sapply(intersectionaucsf2, function(x)
{measurement<-torres[torres$Cancer==x,"SF2"]})
names(cancerradiationin)<-intersectionaucsf2
pancancerradiationauc<-sapply(intersectionaucsf2, function(x)
{measurement<-AUC[AUC$Name==x,"AUC"]})
names(pancancerradiationauc)<-intersectionaucsf2
aucsf2<-merge(cancerradiationin,pancancerradiationauc,by="row.names")
correlationaucsf2<-cor(aucsf2$x,aucsf2$y)
#0.376
# To see the correlations of the SF2 of alberts  paper with AUC metric
intersectionaucsf2a<-intersect(albertsf2$Cancer, AUC$Name)
cancerradiationina<-sapply(intersectionaucsf2a, function(x)
{measurement<-albertsf2[albertsf2$Cancer==x,"SF2"]})
names(cancerradiationina)<-intersectionaucsf2a
pancancerradiationauca<-sapply(intersectionaucsf2a, function(x)
{measurement<-AUC[AUC$Name==x,"AUC"]})
names(pancancerradiationauca)<-intersectionaucsf2a
aucsf2a<-merge(cancerradiationina,pancancerradiationauca,by="row.names")
correlationaucsf2<-cor(aucsf2a$x,aucsf2a$y)
#look at the pan-cancer using cell lines from SF2 measurements to see if smaller values effect the correlations
increasedvectorA<-sapply(torres$Cancer, function(x){
  p<-which(grepl(x, colnames(Data)))
})
allCancersauc<-Data[intersection,unlist(increasedvectorA)]
inter<-intersect(torres$Cancer,cellines)
pancancerradiationA<-sapply(inter, function(x)
{
  measurement<-pancancer[pancancer$Name==x,"AUC"]})
names(pancancerradiationA)=inter
allCancersauc<-t(allCancersauc)
combineauc <- merge(as.data.frame(allCancersauc),pancancerradiationA,by="row.names")
colnames(combineauc)[ncol(combineauc)] = "AUC"

pancancercorrelationA<- sapply(colnames(combineauc)[2:12557],function(x){
  cancer<-combineauc[, c(x,"AUC")]
  significance <- cor.test(cancer[,x],cancer$AUC,use="everything", method=c("pearson") )
  cor<-significance[c(3,4)]
})
pancancercorrelationA<-t(pancancercorrelationA)
pvalue<-pancancercorrelationA[,1]
FDRcorrection<- p.adjust(pvalue, method="fdr" )
pearsonA<- merge(as.data.frame(pancancercorrelationA),FDRcorrection, by="row.names")
library("ggplot2")
cols <- c("red","black")
merge<-ggplot() + 
  geom_point(data=aucsf2, aes(x=aucsf2[,2], y=aucsf2[,3], colour="Torres-Roca"),size=2) + 
  geom_point(data=aucsf2a, aes(x=aucsf2a[,2], y=aucsf2a[,3], colour="Alberts"),size=2)+
  theme_classic()+ labs(title="", x="SF2", y = "AUC")+
  scale_colour_manual(name="SF2 datasets", values = c("black", "grey"))+theme(legend.position="top")+
  annotate("text", x = 0.30, y = 6, label ="R=0.38, n=22 (Torres-Roca)",size=5)+
  annotate("text", x = 0.30, y = 5.8, label ="R=0.45, n=27 (Alberts)",size=5)+
  theme(axis.title=element_text(size=14,face="bold"),aspect.ratio=1,axis.text=element_text(size=14,face="bold"),legend.text=element_text(size=14),legend.title=element_text(size=14))

# Use brewer color palettes
p<-ggplot(aucsf2, aes(aucsf2[,2],aucsf2[,3]) )+
  geom_point(size=2)+
  theme_classic()+ 
  labs(title="",
       x="SF2 (Torres)", y = "AUC")+annotate("text", x = 0.25, y = 6, label = "R=0.38, n=22",size=7)+
  theme(plot.title = element_text(hjust = 0.5),title =element_text(size=14, face='bold'))+
  theme(axis.title=element_text(size=14,face="bold"),aspect.ratio=1,axis.text=element_text(size=14,face="bold"))
# Use brewer color palettes
p+scale_color_brewer(palette="Dark2")
# Use grey scale
p + scale_color_grey()
p1<-ggplot(aucsf2a, aes(aucsf2a[,2],aucsf2a[,3]))+
  geom_point(size=2)+
  theme_classic()+ 
  labs(title="",
       x="Albert SF2", y = "AUC")+annotate("text", x = 0.25, y = 6, label = "R=0.45, n=27",size=7)+
  theme(axis.title=element_text(size=14,face="bold"),aspect.ratio=1,axis.text=element_text(size=12,face="bold"))+
  theme(plot.title = element_text(hjust = 0.5),title =element_text(size=14, face='bold'))
# Use brewer color palettes
p1+scale_color_brewer(palette="Dark2")
# Use grey scale
p1+ scale_color_grey()
#correlation alberts and Torres
intersectionat<-intersect(torres$Cancer, albertsf2$Cancer)
cancerradiationinat<-sapply(intersectionat, function(x)
{measurement<-albertsf2[albertsf2$Cancer==x,"SF2"]})
names(cancerradiationinat)<-intersectionat

cancerradiationint<-sapply(intersectionat, function(x)
{measurement<-torres[torres$Cancer==x,"SF2"]})
names(cancerradiationint)<-intersectionat

atm<-merge(cancerradiationinat,cancerradiationint,by="row.names")
correlationat<-cor(atm$x,atm$y) #0.51

p2<-ggplot(atm, aes(atm[,2],atm[,3]))+
  geom_point(size=2)+
  theme_classic()+ 
  labs(title="",
       x="SF2 (Alberts)", y = "SF2 (Torres)")+annotate("text", x = 0.20, y = 0.77, label = "R=0.51, n=47",size=5)+
  theme(axis.title=element_text(size=14,face="bold"),aspect.ratio=1,axis.text=element_text(size=12,face="bold"))
# Use brewer color palettes
p2+scale_color_brewer(palette="Dark2")
# Use grey scale
p2+ scale_color_grey()

surv<-grid.arrange(p,p1,p2)
ggsave("myfile.pdf", surv,width=15,height=15)


library("gridExtra")

surv<-grid.arrange(merge,p2,nrow=2,ncol=2)
ggsave("AUC.pdf", surv,width=10,height=10)






































