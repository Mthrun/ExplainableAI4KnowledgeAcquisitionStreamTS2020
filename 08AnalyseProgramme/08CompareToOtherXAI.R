#08CompareToOtherXAI.R
Disk="E"
path=ReDi("ExplainableAI4TimeSeries2020/03eUD35",Disk)
source(paste0(path,'/Clustering_Output_eUD35.R'))
setwd(path)

#Compare to eUD3.5 ----
#Apply Source Code of
#https://github.com/miguelmedinaperez/eUD3.5 
#with default parameters to
# ExplainableAI4TimeSeries2020/03eUD35/HydrologyToExplain.lrn
# Outputs of algorithm are stored in
# ExplainableAI4TimeSeries2020/03eUD35
V=Clustering_Output_eUD35("Hydrology_Clusters.txt","Hydrology_Patterns.txt")
Cls=V$Cls+1
Patterns=V$Patterns

Patterns2Rules=function(Patterns){
  Patterns=gsub("\\[",",c(",x = Patterns)
  Patterns=gsub("\\] ","\\)",x = Patterns)
  Patterns=gsub("\\]","\\)",x = Patterns)
  Patterns1List=strsplit(Patterns,",")
  ClassNames=lapply(Patterns1List, function(x) gsub(" ","\\,",x[[2]]))
  ClassPops=lapply(Patterns1List, function(x) gsub(" ","\\,",x[[3]]))
  
  ClassNames_e <- new.env()
  ClassNames=lapply(ClassNames,function(x) eval(parse(text=x),envir = ClassNames_e))
  
  ClassPops_e <- new.env()
  ClassPopss=lapply(ClassPops,function(x) eval(parse(text=x),envir =ClassPops_e))
  
  ClassNamesM=do.call(rbind,ClassNames)
  ClassPopsM=do.call(rbind,ClassPopss)
  Rules=lapply(Patterns1List, "[[",1)
  Names=paste0("R",1:length(Rules))
  names(Rules)=Names
  rownames(ClassNamesM)=Names
  rownames(ClassPopsM)=Names
  
  return(list(ClassPopsM=ClassPopsM,ClassNamesM=ClassNamesM,Rules=Rules))
}

RulesV=Patterns2Rules(Patterns)

NoOfRules=length(RulesV$Rules)
NoOfRules

Disk="E"
setwd(ReDi("ExplainableAI4TimeSeries2020/04DBS",Disk))
load('HydrologieTaeglich_hellinger3Clusters.rda')
ind=which(ClstTrue<4)

table(ClstTrue[ind],Cls)#no overlap
ClusteringAccuracy(Cls,ClstTrue[ind])

#procedure requires https://github.com/aultsch/DataIO
V=ReadLRN(FileName = "HydrologyToExplain.lrn",ReDi("ExplainableAI4TimeSeries2020/03eUD35",Disk))
Data=V$Data
ind=which(ClstTrue<4) #ignore outliers

names=c('C1','C2','C3')
#no clear distinction
ClassMDplot(Data[ind,1],Cls,MinimalAmoutOfData = 10,PlotLegend=F,ClassNames = names,ColorSequence = GeneralizedUmatrix::DefaultColorSequence)$ggobject+theme_bw()+ggtitle('Class MDplot of Nitrate')+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+ylab('PDE')
#c2 and c3 lower than c1
ClassMDplot(Data[ind,7],Cls,MinimalAmoutOfData = 10,PlotLegend=F,ClassNames = names,ColorSequence = GeneralizedUmatrix::DefaultColorSequence)$ggobject+theme_bw()+ggtitle('Class MDplot of Electric Conductivity')+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+ylab('PDE')

## Explainable K-means ----
#procedure requires https://github.com/aultsch/DataIO
V=ReadLRN(FileName = "HydrologyToExplain.lrn",ReDi("ExplainableAI4TimeSeries2020/03eUD35",Disk))
Data=V$Data
ind=which(ClstTrue<4) #ignore outliers

#source code for [Dasgupta et al., 2020]  is not available
#[Dasgupta et al., 2020]  Dasgupta, S., Frost, N., Moshkovitz, M., & Rashtchian, C.: Explainable $ k $-Means and $ k $-Medians Clustering, 37th International Conference on Machine Learning,, Vienna, Austria, 2020.

#But we can try a clustering with k-means of data (except of outliers)
library(FCPS)#on cran, requires version 1.2.6 for kcentroids, otherwise chose another k-means type
ClsV=FCPS::kmeansClustering(Data[ind,],3,Type = "kcentroids")
ClsVeri=ClsV$Cls
require(data.table)#on cran
require(FeatureImpCluster)#on cran
#transforming data to measure feature importance
DF=as.data.frame(Data[ind,])
DF=data.table::as.data.table(DF)
f <- FeatureImpCluster::FeatureImpCluster(ClsVeri$Object,DF)
f$featureImp
table(ClstTrue[ind],ClsVeri)#no overlap
ClusteringAccuracy(ClsVeri,ClstTrue[ind])

#hypothesis
#if features are seperable, they threshold decision tree will work
# names=c("Dry days with warm water of hgl",
#         "Duality",
#         "Dry days with cold water",
#         "Outliers")

names=c('C1','C2','C3')
ClassMDplot(Data[ind,1],ClsVeri,MinimalAmoutOfData = 10,PlotLegend=F,ClassNames = names,ColorSequence = GeneralizedUmatrix::DefaultColorSequence)$ggobject+theme_bw()+ggtitle('Class MDplot of Nitrate')+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+ylab('PDE')
ClassMDplot(Data[ind,7],ClsVeri,MinimalAmoutOfData = 10,PlotLegend=F,ClassNames = names,ColorSequence = GeneralizedUmatrix::DefaultColorSequence)$ggobject+theme_bw()+ggtitle('Class MDplot of Electric Conductivity')+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+ylab('PDE')
#we see that features are equally distributed in the three classes
# therefore it can be assumed that the features cannot be explained by k-means