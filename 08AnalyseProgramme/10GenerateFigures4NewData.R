#10GenerateFigures4NewData.R
## 2015 ----
Disk="E"
setwd(ReDi("ExplainableAI4TimeSeries2020/04DBS",Disk))

load(file="Hydro2015.rda")#,Data2016,DF2016,Trans2016,DecisionTree,Rules,Distance,gmm,outModel,Cls2016,ScriptPfad)

sum(is.finite(intradist2015[intradist2015[,2]>BayesBoundaries2015[2],2]))/sum(is.finite(intradist2015[,2]))
sum(is.finite(intradist2015[intradist2015[,3]>BayesBoundaries2015[2],3]))/sum(is.finite(intradist2015[,3]))

sum(is.finite(intradist2015[intradist2015[,3]>BayesBoundaries2015[1],3]))/sum(is.finite(intradist2015[,3]))

SignedLog(BayesBoundaries2015)
SignedLog(apply(intradist2015, 2, mean,na.rm=TRUE))
SignedLog(apply(intradist2015, 2, max,na.rm=TRUE))

Silhouetteplot(Trans,Cls,main="2015 Data")

Distvec2015=Distance[upper.tri(Distance,diag = F)]
AdaptGauss::PlotMixturesAndBoundaries(Distvec2015,gmm$Means,gmm$SDs,gmm$Weights,xlab='Range of Distances',ylab='PDE',main='GMM of Distribution of Distances df',SingleGausses=T)
legend("topright",c("Distances df","GMM","Single Gaussians","Bayes Boundary"),col=c("black","red","blue","magenta"), lwd=2, lty=c(1,1,1,1),bty = "n")
abline(v = BayesBoundaries2015[2],col="magenta",lwd=5)

par(pty="s")
AdaptGauss::QQplotGMM(Distvec2015,gmm$Means,gmm$SDs,gmm$Weights,ylab='Data = GMM of Distribution of Distances df')

deT2015=XAI::DecisionTree(Data[,c(13,3,6)],Cls = Cls,PlotIt=TRUE,CartPKG = F)
Rules2015=XAI::DecisionTree2Rules(deT2015,digits = 1)
Rules2015$RuleSet

(55+39+27+37+30+21)/nrow(Trans)
Cls=RenameDescendingClassSize(Cls)
Heatmap(Distance,Cls)



TopviewTopographicMap(genU$Umatrix,genU$Bestmatches,Cls = Cls,Imx = imx2015,BmSize = 8)

ClusterCount(Cls)
24/nrow(Data)

Names2015_=c("WarmWater
            
            Without HeavyRain",
"LightRain&MildWater

AtHighLevel", 
"ColdWater

AtHighLevel",
"CoolerWater

AtLowLevel",
"CoolerWater

AtIntermediateLevel",
"HighRainIntensity

AtHighLevel","Outliers")

Names2015=1:7
names(Names2015)=Names2015_
unique(Cls)

ClassMDplot(Data[,"nnit13"],Cls,ClassNames = Names2015[unique(Cls)],MinimalAmoutOfData = 100,Ordering = "Alphabetical")$ggobject+theme_bw()+ggtitle('Class MDplot of Nitrate (NO3)')+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+ylab('NO3	in mg/L')
#ClassMDplot(Data[,"nnit13"],Cls,MinimalAmoutOfData = 100,Ordering = "Average")$ggobject+theme_bw()+ggtitle('Class MDplot of Electric Conductivity')+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+ylab('PDE')

ClassMDplot(Data[,"con47"],Cls,ClassNames = Names2015[unique(Cls)],MinimalAmoutOfData = 100,Ordering = "Alphabetical")$ggobject+theme_bw()+ggtitle('Class MDplot of Electric Conductivity (EC)')+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+ylab('EC in	mS/m')
##2016 ----
Disk="E"
setwd(ReDi("ExplainableAI4TimeSeries2020/04DBS",Disk))

load(file="Hydro2016.rda")#,Data2016,DF2016,Trans2016,DecisionTree,Rules,Distance,gmm,outModel,Cls2016,ScriptPfad)
plotTopographicMap(outModel$Umatrix,outModel$Bestmatches,Cls2016,Imx = imx2016)

TopviewTopographicMap(outModel$Umatrix,outModel$Bestmatches,Cls = Cls2016,Imx = imx2016)

Rules$RuleSet

sum(is.finite(intradist2016[intradist2016[,2]>BayesBoundaries2016[2],2]))/sum(is.finite(intradist2016[,2]))
sum(is.finite(intradist2016[intradist2016[,3]>BayesBoundaries2016[2],3]))/sum(is.finite(intradist2016[,3]))
Heatmap(Trans2016,Cls2016)

ClassMDplot(Data2016[,"nnit13"],Cls2016,MinimalAmoutOfData = 100,Ordering = "Alphabetical")
ClassMDplot(Data2016[,"con47"],Cls2016,MinimalAmoutOfData = 100,Ordering = "Alphabetical")

Distvec=Distance[upper.tri(Distance,diag = F)]
AdaptGauss::PlotMixturesAndBoundaries(Distvec,gmm$Means,gmm$SDs,gmm$Weights,gmm$SDs,gmm$Weights,xlab='Range of Distances',ylab='PDE',main='GMM of Distribution of Distances df',SingleColor = "blue",SingleGausses = T,MixtureColor = "magenta",DataColor = "black",lwd=2)
legend("topright",c("Distances df","GMM","Single Gaussians","Bayes Boundary"),col=c("black","red","blue","magenta"), lwd=2, lty=c(1,1,1,1),bty = "n")
abline(v = BayesBoundaries2016[2],col="magenta",lwd=5)

par(pty="s")
AdaptGauss::QQplotGMM(Distvec,gmm$Means,gmm$SDs,gmm$Weights,ylab='Data = GMM of Distribution of Distances df')
Heatmap(Trans2016,Cls2016)
Silhouetteplot(Trans2016,Cls2016,main="2016 Data")

deT2016=XAI::DecisionTree(Data2016[,c(13,3,6)],Cls = Cls2016,PlotIt=TRUE,CartPKG = F)

EvolutionTree <- function(Data, Cls){
  Cls <- factor(Cls)
  fullData = cbind(data.frame(Data) , Cls)
  fit <- evtree::evtree(Cls ~ ., fullData)
  plot(fit)
  return(fit)
}
EvolutionTree(Data2016[,c(13,3,6)],Cls = Cls2016)
ind=which(Data2016[,3]<1.137)
ind=which(Data2016[,3]>=1.137&Data2016[,c(6)]>8.706)
Cls3=Cls2016
Cls3[ind]=3

Rules2016=XAI::DecisionTree2Rules(deT2016,digits = 1)

Rules2016$RuleSet

ClassMDplot(Data2016[,"rain"],Cls2016,MinimalAmoutOfData = 100,Ordering = "Alphabetical")

out=ClassPDEplotMaxLikeli(Data2016[,'rain'],Cls2016,lwd=1.5)

out$ggobject+ggtitle('Class PDEplot of rain')+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+
  ylab('PDE')+xlab('Wt18 in ?C')+
  geom_vline(xintercept = 12.5, color = "red",size=1.5,linetype='dashed')


out=ClassPDEplotMaxLikeli(Data2016[,'Wt18'],Cls2016,lwd=1.5)

out$ggobject+ggtitle('Class PDEplot of Water Temperature')+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+
  ylab('PDE')+xlab('Wt18 in ?C')+
  geom_vline(xintercept = 9.1, color = "red",size=1.5,linetype='dashed')

ClusterCount(Cls2016)


32/nrow(Data2016)


Names2016_=c("ColdWater
            
            AtLowerLevel","WarmWater
            
            AtHigherLevel")
Names2016=c(1:2)
names(Names2016)=Names2016_
ClassMDplot(Data2016[,"nnit13"],Cls2016,ClassNames = Names2016,MinimalAmoutOfData = 100,Ordering = "Alphabetical")$ggobject+theme_bw()+ggtitle('Class MDplot of Nitrate (NO3)')+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+ylab('PDE of NO3	in mg/L')
ClassMDplot(Data2016[,"con47"],Cls2016,ClassNames = Names2016,MinimalAmoutOfData = 100,Ordering = "Alphabetical")$ggobject+theme_bw()+ggtitle('Class MDplot of Electric Conductivity (EC)')+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+ylab('PDE of EC in	mS/m')
