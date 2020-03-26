#07OckhamsRazor.R
##
library(PPCI)#cran
library(rgl)#cran
library(GeneralizedUmatrix)#cran
library(ProjectionBasedClustering)#cran
##nonlinear model proof ----
setwd(ReDi("ExplainableAI4KnowledgeAcquisitionStreamTS2020/04DBS"))

load('HydrologieTaeglich_hellinger3Clusters.rda')

proj=ProjectionBasedClustering::ICA(Trans4)
genmodellinear=GeneralizedUmatrix(Trans3,proj$ProjectedPoints,F,Cls=ClstTrue)  
plotTopographicMap(genmodellinear$Umatrix,genmodellinear$Bestmatches,Cls=ClstTrue)

setwd(ReDi("ExplainableAI4KnowledgeAcquisitionStreamTS2020/07Dokumentationen"))
rgl.snapshot('HydrologyStructuresLinearhellinger3Clusters.png')


proj=PPCI::mcdr(Trans4,p = 2,maxit=1000,ftol=1e-7)

plot(proj$fitted,col=ClstTrue)

genmodellinear=GeneralizedUmatrix(Trans3,proj$fitted,F,Cls=ClstTrue)  
plotTopographicMap(genmodellinear$Umatrix,genmodellinear$Bestmatches,Cls=ClstTrue,BmSize = 1.5)

setwd(ReDi("ExplainableAI4KnowledgeAcquisitionStreamTS2020/07Dokumentationen"))
rgl.snapshot('HydrologyStructuresProjectionPursuit.png')



#Trans4=MDS(DistanceMatrix(Trans3,'hellinger'),OutputDimension = ncol(Trans3))$ProjectedPoints

#library(PPCI)
#proj=PPCI::mcdr(Trans4,p = 2,maxit=1000,ftol=1e-7)

#plot(proj$fitted,col=ClstTrue)
#Proj=proj$fitted
#Proj[,2]=Proj[,2]*10
#plot(Proj,col=ClstTrue)
#genmodellinear=GeneralizedUmatrix(Trans4,Proj,F,Cls=ClstTrue,Tiled = T)  
#plotTopographicMap(genmodellinear$Umatrix,genmodellinear$Bestmatches,Cls=ClstTrue,BmSize = 1.2)

