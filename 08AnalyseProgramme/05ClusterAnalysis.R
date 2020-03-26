#05ClusterAnalysis.R
library(DatabionicSwarm)
library(GeneralizedUmatrix)
library(ProjectionBasedClustering)
library(DataVisualizations)
setwd(ReDi("ExplainableAI4KnowledgeAcquisitionStreamTS2020/01Transformierte"))
load(file='hellinger_DistancesHydrologie.rda')#,pvalueChitest,distvec,gmm,InputDistances,Trans3,Trans,backtrafo)
#Clustering ----

res=DatabionicSwarm::Pswarm(InputDistances)

Trans4=ProjectionBasedClustering::MDS(InputDistances,OutputDimension = ncol(Trans3))$ProjectedPoints
resUmatrix=DatabionicSwarm::GeneratePswarmVisualization(Trans4,res$ProjectedPoints,res$LC)
GeneralizedUmatrix::plotTopographicMap(resUmatrix$Umatrix,resUmatrix$Bestmatches,NoLevels=10)
norm=GeneralizedUmatrix::NormalizeUmatrix(Trans4,resUmatrix$Umatrix,resUmatrix$Bestmatches)
GeneralizedUmatrix::plotTopographicMap(norm,resUmatrix$Bestmatches,NoLevels=10,Cls)

Cls=DatabionicSwarm:DBSclustering(3,Trans4,resUmatrix$Bestmatches,res$LC,PlotIt = T,StructureType = T)

Cls2=ProjectionBasedClustering::interactiveClustering(resUmatrix$Umatrix,resUmatrix$Bestmatches,Cls = Cls)$Cls

imx=ProjectionBasedClustering::interactiveGeneralizedUmatrixIsland(resUmatrix$Umatrix,resUmatrix$Bestmatches,ClstTrue2)

ClstTrue=RenameDescendingClassSize(Cls2)

plotTopographicMap(resUmatrix$Umatrix,resUmatrix$Bestmatches,ClstTrue,BmSize = 0.8,Imx = imx)
library(rgl)
setwd(ReDi("HydrologieSchwarmClustering2016/07Dokumentationen"))
rgl.snapshot('hellinger3cluster_HydrologyStructures.png')

Silhouette(Trans4,ClstTrue)
Heatmap(InputDistances,ClstTrue)+ggtitle('Distances sorted by Clustering (Cls)')
cc=ClstTrue;cc[cc>4]=4;Heatmap(InputDistances,cc)+ggtitle('Distances sorted by clustering (Cls)')
intracluster=ClusterDistances(InputDistances,Cls = ClstTrue)
MDplot(intracluster)
apply(intracluster,2,median,na.rm=T)

Heatmap(InputDistances,ClstTrue)+ggtitle('Distances sorted by Clustering (Cls)')
requireRpackage('fpc')
fpc::cluster.stats(d = InputDistances, clustering=ClstTrue,silhouette=T)

method='hellinger'
setwd(ReDi("ExplainableAI4KnowledgeAcquisitionStreamTS2020/04DBS"))
save(file='HydrologieTaeglich_hellinger3Clusters.rda',ClstTrue,Trans3,InputDistances,Trans4,method,imx,res,resUmatrix,norm,backtrafo,Trans,Time)