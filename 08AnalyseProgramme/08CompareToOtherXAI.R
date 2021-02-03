#08CompareToOtherXAI.R
Disk="E"
path=ReDi("ExplainableAI4TimeSeries2020/03eUD35",Disk)
source(paste0(path,'/Clustering_Output_eUD35.R'))
setwd(path)
#function used to prepare eUD3.5 patterns
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

#Compare to eUD3.5 ----
#Apply Source Code of
#https://github.com/miguelmedinaperez/eUD3.5 
#with default parameters to
# ExplainableAI4TimeSeries2020/03eUD35/HydrologyToExplain.lrn
# Outputs of algorithm are stored in
# ExplainableAI4TimeSeries2020/03eUD35

#2013-2014
V=Clustering_Output_eUD35("Hydrology_Clusters.txt","Hydrology_Patterns.txt")
Cls=V$Cls+1
Patterns=V$Patterns

RulesV=Patterns2Rules(Patterns)

NoOfRules=length(RulesV$Rules)
NoOfRules #541

V=Clustering_Output_eUD35("Hydrology_Clusters2015.txt","Hydrology_Paterns2015.txt")
Cls2015=V$Cls+1#234
Patterns=V$Patterns

RulesV=Patterns2Rules(Patterns)

NoOfRules=length(RulesV$Rules)
NoOfRules #503

V=Clustering_Output_eUD35("Hydrology_Clusters2016.txt","Hydrology_Paterns2016.txt")
Cls2016=V$Cls+1#291
Patterns=V$Patterns

RulesV=Patterns2Rules(Patterns)

NoOfRules=length(RulesV$Rules)
NoOfRules #552

#procedure requires https://github.com/aultsch/DataIO
V=ReadLRN(FileName = "HydrologyToExplain.lrn",ReDi("ExplainableAI4TimeSeries2020/03eUD35",Disk))
Data=V$Data
ind=which(ClstTrue<4) #ignore outliers

Names=c(1:3)
names(Names)=c('C 1','C 2','C 3')
#no clear distinction
ClassMDplot(Data[ind,1],Cls,MinimalAmoutOfData = 10,PlotLegend=F,ClassNames = Names,ColorSequence = GeneralizedUmatrix::DefaultColorSequence,Ordering = "Alphabetical")$ggobject+theme_bw()+ggtitle('Class MDplot of Nitrate (NO3) for 2013/2014')+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+ylab('PDE of NO3 in mg/L')+xlab("eUD3.5 Classes")

ClassMDplot(Data[ind,7],Cls,MinimalAmoutOfData = 10,PlotLegend=F,ClassNames = Names,ColorSequence = GeneralizedUmatrix::DefaultColorSequence,Ordering = "Alphabetical")$ggobject+theme_bw()+ggtitle('Class MDplot of Electric Conductivity (EC) for 2013/2014')+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+ylab('PDE of EC in	mS/m')+xlab("eUD3.5 Classes")


V=ReadDates("20201226Hydrologie2016",ReDi("Hydrologie2019/09Originale",Disk))
Data2016=V$Data

ClassMDplot(Data2016[,1],Cls2016,MinimalAmoutOfData = 10,PlotLegend=F,ClassNames = Names,ColorSequence = GeneralizedUmatrix::DefaultColorSequence,Ordering = "Alphabetical")$ggobject+theme_bw()+ggtitle('Class MDplot of Nitrate (NO3) for 2016')+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+ylab('PDE of NO3 in mg/L')+xlab("eUD3.5 Classes")

ClassMDplot(Data2016[,7],Cls2016,MinimalAmoutOfData = 10,PlotLegend=F,ClassNames = Names,ColorSequence = GeneralizedUmatrix::DefaultColorSequence,Ordering = "Alphabetical")$ggobject+theme_bw()+ggtitle('Class MDplot of Electric Conductivity (EC) for 2016')+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+ylab('PDE of EC in	mS/m')+xlab("eUD3.5 Classes")

V=ReadDates("20201226Hydrologie2015",ReDi("Hydrologie2019/09Originale",Disk))
Data2015=V$Data

ClassMDplot(Data2015[,1],Cls2015,MinimalAmoutOfData = 10,PlotLegend=F,ClassNames = Names,ColorSequence = GeneralizedUmatrix::DefaultColorSequence,Ordering = "Alphabetical")$ggobject+theme_bw()+ggtitle('Class MDplot of Nitrate (NO3) for 2015')+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+ylab('PDE of NO3 in mg/L')+xlab("eUD3.5 Classes")

ClassMDplot(Data2015[,7],Cls2015,MinimalAmoutOfData = 10,PlotLegend=F,ClassNames = Names,ColorSequence = GeneralizedUmatrix::DefaultColorSequence,Ordering = "Alphabetical")$ggobject+theme_bw()+ggtitle('Class MDplot of Electric Conductivity (EC) for 2015')+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+ylab('PDE of EC in	mS/m')+xlab("eUD3.5 Classes")

## IMM ----
Disk="E"
setwd(ReDi("ExplainableAI4TimeSeries2020/04PBC",Disk))
load('HydrologieTaeglich_hellinger3Clusters.rda')
ind=which(ClstTrue<4)

table(ClstTrue[ind],Cls)#no overlap
ClusteringAccuracy(Cls,ClstTrue[ind])

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

Names=c(1:3)
names(Names)=c('C 1','C 2','C 3')
ClassMDplot(Data[ind,1],ClsVeri,MinimalAmoutOfData = 10,PlotLegend=F,ClassNames = Names,ColorSequence = GeneralizedUmatrix::DefaultColorSequence,Ordering = "Alphabetical")$ggobject+theme_bw()+ggtitle('Class MDplot of Nitrate (NO3) for 2013/2014')+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+ylab('PDE of NO3 in mg/L')+xlab("IMM Classes")
ClassMDplot(Data[ind,7],ClsVeri,MinimalAmoutOfData = 10,PlotLegend=F,ClassNames = Names,ColorSequence = GeneralizedUmatrix::DefaultColorSequence,Ordering = "Alphabetical")$ggobject+theme_bw()+ggtitle('Class MDplot of Electric Conductivity (EC) for 2013/2014')+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+ylab('PDE of EC in	mS/m')+xlab("IMM Classes")
#we see that features are equally distributed in the three classes
# therefore it can be assumed that the features cannot be explained by k-means


## IMM on 2016 ----
Disk="E"
setwd(ReDi("ExplainableAI4TimeSeries2020/04PBC",Disk))

load(file="Hydro2016.rda")#,Data2016,DF2016,Trans2016,DecisionTree,Rules,Distance,gmm,outModel,Cls2016,ScriptPfad)

ClsV=FCPS::kmeansClustering(Data2016,2,Type = "kcentroids")

ClsVeri=ClsV$Cls
require(data.table)#on cran
require(FeatureImpCluster)#on cran
#transforming data to measure feature importance
DF=as.data.frame(Data2016)
DF=data.table::as.data.table(DF)
f <- FeatureImpCluster::FeatureImpCluster(ClsV$Object,DF)
f$featureImp
table(Cls2016,ClsVeri)#no overlap
ClusteringAccuracy(ClsVeri,Cls2016)

Names=c(1:2)
names(Names)=c('C 1','C 2')
ClassMDplot(Data2016[,1],ClsVeri,MinimalAmoutOfData = 100,PlotLegend=F,ClassNames = Names,ColorSequence = GeneralizedUmatrix::DefaultColorSequence,Ordering = "Alphabetical")$ggobject+theme_bw()+ggtitle('Class MDplot of Nitrate (NO3) for 2016')+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+ylab('PDE of NO3 in mg/L')+xlab("IMM Classes")
ClassMDplot(Data2016[,7],ClsVeri,MinimalAmoutOfData = 100,PlotLegend=F,ClassNames = Names,ColorSequence = GeneralizedUmatrix::DefaultColorSequence,Ordering = "Alphabetical")$ggobject+theme_bw()+ggtitle('Class MDplot of Electric Conductivity (EC) for 2016')+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+ylab('PDE of EC in	mS/m')+xlab("IMM Classes")

# 
# sol71        At47        Wt13        Wt18     
# 0.428865979 0.002405498 0.001718213 0.001030928
## IMM on 2015 ----
load(file="Hydro2015.rda")#,Data2016,DF2016,Trans2016,DecisionTree,Rules,Distance,gmm,outModel,Cls2016,ScriptPfad)
ind=which(Cls<7)
ClsV=FCPS::kmeansClustering(Data[ind,],6,Type = "kcentroids")

ClsVeri=ClsV$Cls
require(data.table)#on cran
require(FeatureImpCluster)#on cran
#transforming data to measure feature importance
DF=as.data.frame(Data[ind,])
DF=data.table::as.data.table(DF)
f <- FeatureImpCluster::FeatureImpCluster(ClsV$Object,DF)
f$featureImp

# 
# sol71         q13         q18        At47      
# 0.714410480 0.090829694 0.067248908 0.003056769 
table(Cls[ind],ClsVeri)#no overlap
ClusteringAccuracy(ClsVeri,Cls[ind])


Names=c(1:6)
names(Names)=c('C 1','C 2','C 3','C 4','C 5','C 6')
ClassMDplot(Data[ind,1],ClsVeri,MinimalAmoutOfData = 100,PlotLegend=F,ClassNames = Names,ColorSequence = GeneralizedUmatrix::DefaultColorSequence,Ordering = "Alphabetical")$ggobject+theme_bw()+ggtitle('Class MDplot of Nitrate (NO3) for 2016')+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+ylab('PDE of NO3 in mg/L')+xlab("IMM Classes")
ClassMDplot(Data[ind,7],ClsVeri,MinimalAmoutOfData = 100,PlotLegend=F,ClassNames = Names,ColorSequence = GeneralizedUmatrix::DefaultColorSequence,Ordering = "Alphabetical")$ggobject+theme_bw()+ggtitle('Class MDplot of Electric Conductivity (EC) for 2016')+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+ylab('PDE of EC in	mS/m')+xlab("IMM Classes")

