#05ClusterAnalysis.R
require(DatabionicSwarm)
require(GeneralizedUmatrix)
require(ProjectionBasedClustering)
library(DataVisualizations)
setwd(ReDi("ExplainableAI4TimeSeries2020/01Transformierte"))
load(file='hellinger_DistancesHydrologie.rda')#,pvalueChitest,distvec,gmm,InputDistances,Trans3,Trans,backtrafo)
#Clustering ----

#Module 1: Projection
res=DatabionicSwarm::Pswarm(InputDistances)

#Transformation of Distances to Data
Trans4=ProjectionBasedClustering::MDS(InputDistances,OutputDimension = ncol(Trans3))$ProjectedPoints

#Model 2: Compute U-Matrix for Visualization
resUmatrix=DatabionicSwarm::GeneratePswarmVisualization(Trans4,res$ProjectedPoints,res$LC)
#Evaluate Number of Clusters in Visualization
GeneralizedUmatrix::plotTopographicMap(resUmatrix$Umatrix,resUmatrix$Bestmatches,NoLevels=10)
#Normalize U-Matrix with Abstract U-Matrix
norm=GeneralizedUmatrix::NormalizeUmatrix(Trans4,resUmatrix$Umatrix,resUmatrix$Bestmatches)
#Evaluate Number of Clusters in Visualization
GeneralizedUmatrix::plotTopographicMap(norm,resUmatrix$Bestmatches,NoLevels=10,Cls)

#3. Module: Clustering, No. of Clusters is No. of Valleys
Cls=DatabionicSwarm:DBSclustering(3,Trans4,resUmatrix$Bestmatches,res$LC,PlotIt = T,StructureType = T)

#Detect Outliers Interactively
Cls2=ProjectionBasedClustering::interactiveClustering(resUmatrix$Umatrix,resUmatrix$Bestmatches,Cls = Cls)$Cls

#Manually Cut-out Island
imx=ProjectionBasedClustering::interactiveGeneralizedUmatrixIsland(resUmatrix$Umatrix,resUmatrix$Bestmatches,ClstTrue2)

#Rename Cluster Labels
ClstTrue=FCPS::ClusterRenameDescendingSize(Cls2)

method='hellinger'
setwd(ReDi("ExplainableAI4TimeSeries2020/04PBC"))
save(file='HydrologieTaeglich_hellinger3Clusters.rda',Trans3,InputDistances,Trans4,method,imx,res,resUmatrix,norm,backtrafo,Trans,Time)

## Ploting ----
setwd(ReDi("ExplainableAI4TimeSeries2020/04PBC"))
load(file='HydrologieTaeglich_hellinger3Clusters.rda')

library(rgl)
r3dDefaults$windowRect = c(0,0,1200,1200) 
cc=ClstTrue;cc[cc>4]=4;
unique(cc)
plotTopographicMap(resUmatrix$Umatrix,resUmatrix$Bestmatches,cc,BmSize = 0.8,Imx = imx,ClsColors=c("black", "red", "yellow", "magenta"),
                   Names=c("Cluster 1","Cluster 2","Cluster 3","Outliers"),NamesPosition="bottomright",NamesCex=2)

setwd(ReDi("HydrologieSchwarmClustering2016/07Dokumentationen"))
rgl.snapshot('hellinger3cluster_HydrologyStructures.png')

##Evaluation ----
require(FCPS)
library(DataVisualizations)
DataVisualizations::Silhouetteplot(Trans4,ClstTrue,main="2013/2014 data")
DataVisualizations::Heatmap(InputDistances,ClstTrue)+ggtitle('Distances sorted by Clustering (Cls)')
cc=ClstTrue;cc[cc>4]=4;Heatmap(InputDistances,cc)+ggtitle('Distances sorted by clustering (Cls)')
intracluster=FCPS::ClusterDistances(InputDistances,Cls = ClstTrue)
DataVisualizations::MDplot(intracluster)
apply(intracluster,2,median,na.rm=T)

DataVisualizations::Heatmap(InputDistances,ClstTrue)+ggtitle('Distances sorted by Clustering (Cls)')
require('fpc')
fpc::cluster.stats(d = InputDistances, clustering=ClstTrue,silhouette=T)


setwd(ReDi("ExplainableAI4TimeSeries2020/04PBC"))
save(file='HydrologieTaeglich_hellinger3Clusters.rda',ClstTrue,Trans3,InputDistances,Trans4,method,imx,res,resUmatrix,norm,backtrafo,Trans,Time)