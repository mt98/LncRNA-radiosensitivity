#correlations between the expression values of the bladder and the head and neck cell lines with the  AUC measurements calculated
setwd("/Users/mkhan/Documents/ third year")
#load CCLE data
load("CCLE.RData")
#load TCGA data to identify lncRNAs which are in the CCLE list
files <- list.files(path="/Users/mkhan/Documents/ third year/TCGA data", pattern="*rnaexpr.tsv",full.names=T, recursive=TRUE)
Data1 <- lapply(files
                ,function(x){
                  t <- read.table(x, header=T, sep="\t", row.names=1)
                }
)# reading the TCGA files
Cancers <-lapply(Data1
                 ,function(x){
                   Patient <-lapply(colnames(x), function(y){strsplit(y, fixed=T, ".")[[1]][[2]]})
                   Patient <- unlist(Patient)
                   Tumour <- Patient == "Tumor"
                   Cancer <- x[, Tumour]
                 }  
)# remove the samples which are normal
CombinedCancers<-do.call('cbind', Cancers)
ENSEMBL<-rownames(CombinedCancers)# ensembl ids for the TCGA TANRIC dataset
LncRNAs<-rownames(Data)# ensembl ids for the lncRNAs
intersectionlncRNAs<-intersect(ENSEMBL,LncRNAs)# intersection of the ensembl IDs to give only lncRNA expression levels
FilteredData<-Data[intersectionlncRNAs,]# lncrna expression levels
FilteredData<-log2(FilteredData+1)#log transform

#Read in the radiosensitivity measurement data file and then look at the AUC values
radiosensitivity<-read.table("/Users/mkhan/Documents/ third year/CCLE data/info files//anotated.csv", header=T, sep=",",stringsAsFactors = FALSE)
Names=radiosensitivity$Name
radiation<-sapply (Names, function(x){
measurement<-radiosensitivity[radiosensitivity$Name==x,"AUC"]})
names(radiation)=Names
library("plyr")

#frequency of cell lines for each tumour type
frequency<-count(radiosensitivity,"Cancer") 
select<-subset(frequency,freq>=15) #select only those cell lines which have cell lines greater than 15
sites<-select$Cancer 

#correlations for all the cancer types individually 
individualcancersprimary<-lapply(sites,function(x){
    #for each individual cancer type
    primarycancer <- subset(radiosensitivity, Cancer== x)
    rownames(primarycancer)<-primarycancer$Name
    #find the lncRNAs in the CCLE data 
    intersection<-intersect(rownames(primarycancer),colnames(FilteredData))
    individualtype<-FilteredData[,intersection]
   #remove the lncRNAs which are zero across all samples
    individualtype<- individualtype[rowSums(individualtype) > 0, ]
    combine = merge(as.data.frame(t(individualtype)), radiation, by="row.names" )
    rownames(combine)=combine$Row.names
                    colnames(combine)[ncol(combine)] = "AUC"
        #calculate the correlations between the expression levels and AUC           
                  correlations<- sapply(colnames(combine)[2:ncol(individualtype)],function(z){
          cancertype<-combine[, c(z,"AUC")]
          colnames(cancertype)[1]=z 
          significance <- cor.test(cancertype[,z],cancertype$AUC,use="everything", method=c("pearson") )
          cor<-significance[c(3,4)]
         
          })
                  correlations<-t(correlations)
                  pvalue<-correlations[,1]
                  FDRcorrection<- p.adjust(pvalue, method="fdr" )
                  pearson<- merge(as.data.frame(correlations),FDRcorrection, by="row.names")
                  
                  })

write.table( as.matrix(individualcancersprimary[1][[1]]),file="breastcancers.csv",quote=F,sep=",")
write.table( as.matrix(individualcancersprimary[2][[1]]),file="endometriumcancers.csv",quote=F,sep=",")
write.table( as.matrix(individualcancersprimary[3][[1]]),file="glioblastoma.csv",quote=F,sep=",")
write.table( as.matrix(individualcancersprimary[15][[1]]),file="bladder.csv",quote=F,sep=",")
write.table( as.matrix(individualcancersprimary[14][[1]]),file="hnsc.csv",quote=F,sep=",")
write.table( as.matrix(individualcancersprimary[5][[1]]),file="colorectal.csv",quote=F,sep=",")

#look at the significant correlations based on p-value
significantindividualcancersprimary<-lapply(sites,function(x){
  primarycancer <- subset(radiosensitivity, Cancer== x)
  rownames(primarycancer)<-primarycancer$Name
  intersection<-intersect(rownames(primarycancer),colnames(FilteredData))
  individualtype<-FilteredData[,intersection]
  individualtype<- individualtype[rowSums(individualtype) > 0, ]
  combine = merge(as.data.frame(t(individualtype)), radiation, by="row.names" )
  rownames(combine)=combine$Row.names
  colnames(combine)[ncol(combine)] = "AUC"
  correlations<- sapply(colnames(combine)[2:ncol(individualtype)],function(z){
    cancertype<-combine[, c(z,"AUC")]
    colnames(cancertype)[1]=z 
    significance <- cor.test(cancertype[,z],cancertype$AUC,use="everything", method=c("pearson") )
    cor<-significance[c(3,4)]})
  correlations<-t(correlations)
  pvalue<-correlations[,1]
  FDRcorrection<- p.adjust(pvalue, method="fdr" )
  pearson<- merge(as.data.frame(correlations),FDRcorrection, by="row.names")
  filtration<-subset(pearson,p.value<0.05)
})


write.table( as.matrix(significantindividualcancersprimary[1][[1]]),file="pavaluebreastcancerspvalue.csv",quote=F,sep=",")
write.table( as.matrix(significantindividualcancersprimary[2][[1]]),file="pvalueendometriumcancers.csv",quote=F,sep=",")
write.table( as.matrix(significantindividualcancersprimary[14][[1]]),file="pvaluehead-neck.csv",quote=F,sep=",")
write.table( as.matrix(significantindividualcancersprimary[15][[1]]),file="pvaluebladder.csv",quote=F,sep=",")
write.table( as.matrix(significantindividualcancersprimary[5][[1]]),file="pvaluecolo.csv",quote=F,sep=",")
#FILTER OUT LNCRNAS WITH THE LOW EXPRESSION using MAD and calculate correlations only for head and neck and bladder 
madindividualcancersprimary<-sapply(sites[c(14,15)],function(x){
  primarycancer <- subset(radiosensitivity, Cancer== x)
  rownames(primarycancer)<-primarycancer$Name
  intersection<-intersect(rownames(primarycancer),colnames(FilteredData))
  #we study the MAD (median absolute deviation)
  filteredindividualtype<-FilteredData[,intersection]
  combine = merge(as.data.frame(t(filteredindividualtype)), radiation, by="row.names" )
  rownames(combine)=combine$Row.names
  colnames(combine)[ncol(combine)] = "AUC"
medianabsolutedeviation<-sapply(colnames(combine)[2:12557],function(z){
median<-mad(combine[,z], center = median(combine[,z]), constant = 1.4826, na.rm = FALSE,
              low = FALSE, high = FALSE)
                                  })
})
#write.table(madindividualcancersprimary,file="Medvalues.csv",quote=F,sep=",")
colnames(madindividualcancersprimary)<-c("headneck","bladder")
madindividualcancersprimary<-as.data.frame(madindividualcancersprimary)
madindividualcancersprimary$headneckrank<-rank(-madindividualcancersprimary$headneck)
madindividualcancersprimary$bladderrank<-rank(-madindividualcancersprimary$bladder)
#write.table(madindividualcancersprimary,file="Medvalues.csv",quote=F,sep=",")
HNSC<-madindividualcancersprimary[madindividualcancersprimary$headneckrank<=2000,]
Bladder<-madindividualcancersprimary[madindividualcancersprimary$bladderrank<=2000,]
#HNSC correlations and significance
primarycancerHNSC <- subset(radiosensitivity, Cancer== sites[14])
  rownames(primarycancerHNSC)<-primarycancerHNSC$Name
  intersectionHNSC<-intersect(rownames(primarycancerHNSC),colnames(FilteredData))
  individualtypeHNSC<-FilteredData[rownames(HNSC),intersectionHNSC]
  combineHNSC = merge(as.data.frame(t(individualtypeHNSC)), radiation, by="row.names" )
  rownames(combineHNSC)=combineHNSC$Row.names
  colnames(combineHNSC)[ncol(combineHNSC)] = "AUC"
  correlationsHNSC<- sapply(colnames(combineHNSC)[2:2001],function(z){
    cancertype<-combineHNSC[, c(z,"AUC")]
    colnames(cancertype)[1]=z 
    significance <- cor.test(cancertype[,z],cancertype$AUC,use="everything", method=c("pearson") )
    cor<-significance[c(3,4)]
    
  })
  correlationsHNSC<-t(correlationsHNSC)
  pvalue<-correlationsHNSC[,1]
  FDRcorrection<- p.adjust(pvalue, method="fdr" )
  pearsonHNSC<- merge(as.data.frame(correlationsHNSC),FDRcorrection, by="row.names")
  rownames(pearsonHNSC)<-colnames(combineHNSC)[2:2001]
  #write.table( as.matrix(pearsonHNSC),file="headandneckaftermediancorrectiontop200070118.csv",quote=F,sep=",")  
#bladder cancer correlation
  primarycancerbladder <- subset(radiosensitivity, Cancer== sites[15])
  rownames(primarycancerbladder)<-primarycancerbladder$Name
  intersectionbladder<-intersect(rownames(primarycancerbladder),colnames(FilteredData))
  individualtypebladder<-FilteredData[rownames(Bladder),intersectionbladder]
  combinebladder = merge(as.data.frame(t(individualtypebladder)), radiation, by="row.names" )
  rownames(combinebladder)=combinebladder$Row.names
  colnames(combinebladder)[ncol(combinebladder)] = "AUC"
  correlationsbladder<- sapply(colnames(combinebladder)[2:2001],function(z){
    cancertype<-combinebladder[, c(z,"AUC")]
    colnames(cancertype)[1]=z 
    significance <- cor.test(cancertype[,z],cancertype$AUC,use="everything", method=c("pearson") )
    cor<-significance[c(3,4)]
    
  })
  correlationsbladder<-t(correlationsbladder)
  pvaluebladder<-correlationsbladder[,1]
  FDRcorrectionbladder<- p.adjust(pvaluebladder, method="fdr" )
  pearsonbladder<- merge(as.data.frame(correlationsbladder),FDRcorrectionbladder, by="row.names")
  rownames(pearsonbladder)<-colnames(combinebladder)[2:2001] 
  #write.table( as.matrix(pearsonbladder),file="bladderaftermediancorrectionnew.csv",quote=F,sep=",")  
# calculate the range of the correlations to tabulate in a table  
  individualcancersprimarytest<-lapply(sites,function(x){
    primarycancer <- subset(radiosensitivity, Cancer== x)
    rownames(primarycancer)<-primarycancer$Name
    intersection<-intersect(rownames(primarycancer),colnames(FilteredData))
    individualtype<-FilteredData[,intersection]
    individualtype<- individualtype[rowSums(individualtype) > 0, ]
    combine = merge(as.data.frame(t(individualtype)), radiation, by="row.names" )
    rownames(combine)=combine$Row.names
    colnames(combine)[ncol(combine)] = "AUC"
    correlations<- sapply(colnames(combine)[2:ncol(individualtype)],function(z){
      cancertype<-combine[, c(z,"AUC")]
      colnames(cancertype)[1]=z 
      significance <- cor.test(cancertype[,z],cancertype$AUC,use="everything", method=c("pearson") )
      cor<-significance[c(3,4)]  })
    correlations<-t(correlations)
    pvalue<-correlations[,1]
    FDRcorrection<- p.adjust(pvalue, method="fdr" )
    pearson<- merge(as.data.frame(correlations),FDRcorrection, by="row.names")
    pearson1<-subset(pearson,pearson$p.value<0.05) })
  par(mfrow = c(5, 3))
  range<-sapply(individualcancersprimarytest,function(x){
    P<-range(x$estimate,na.rm=FALSE)
  })
  dim(range)
  range<-t(range)
  rownames(range)<-sites
  rangeFDR<-sapply(individualcancersprimarytest,function(x){
    P<-range(x$y,na.rm=FALSE)
  })
  rangeFDR<-t(rangeFDR)
  rownames(rangeFDR)<-sites
  p<-cbind(range,rangeFDR)
  
range(individualcancersprimarytest[2][[1]]$estimate,na.rm=FALSE)

###look at the differences in AUC values between using the mean of the 24 vs mean of the ten
significant<-pearsonbladder[pearsonbladder$estimate>=0.6,]
#expression score calculation using the 24 genes
primarycancer <- subset(radiosensitivity, Cancer== sites[15])
  rownames(primarycancer)<-primarycancer$Name
  intersection<-intersect(rownames(primarycancer),colnames(FilteredData))
  filteredindividualtype<-FilteredData[rownames(significant),intersection]
mean_24<-apply(filteredindividualtype,2,mean)  
#merge with AUC  
combine = merge(((mean_24)), radiation, by="row.names" )
  rownames(combine)=combine$Row.names
  colnames(combine)[ncol(combine)-1] = "mean24"
  colnames(combine)[ncol(combine)] = "AUC"
med<-median(combine$mean24)
strat<-ifelse(combine$mean24>=med,"high","low")
combine$strat<-strat
highAUC<-combine[combine$strat=="high","AUC"]
lowAUC<-combine[combine$strat=="low","AUC"]
t.test(highAUC,lowAUC)

##using the mean of ten lncRNAs
ten<-c( "ENSG00000179743.2", "ENSG00000186056.5", "ENSG00000204623.4" , "ENSG00000223573.2", "ENSG00000227674.1", 
 "ENSG00000234684.2", "ENSG00000235560.3", "ENSG00000237836.1" , "ENSG00000245711.2", "ENSG00000261716.1" )

filtered_ten<-filteredindividualtype[ten,]
###
mean_10<-apply(filtered_ten,2,mean)  
#merge with AUC  
combine_10 = merge(((mean_10)), radiation, by="row.names" )
  rownames(combine_10)=combine_10$Row.names

  colnames(combine_10)[ncol(combine_10)] = "AUC"
med_10<-median(combine_10$x)
strat<-ifelse(combine_10$x>=med_10,"high","low")
combine_10$strat<-strat
highAUC<-combine_10[combine_10$strat=="high","AUC"]
lowAUC<-combine_10[combine_10$strat=="low","AUC"]
t.test(highAUC,lowAUC)
cor.test(combine_10$x,combine_10$AUC)
 
 
  
 