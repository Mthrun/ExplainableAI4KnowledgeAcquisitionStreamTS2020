#07OckhamsRazor.R
##
library(PPCI)#cran
library(rgl)#cran
library(GeneralizedUmatrix)#cran
library(ProjectionBasedClustering)#cran
##nonlinear model proof ----
setwd(ReDi("ExplainableAI4TimeSeries2020/04DBS"))

load('HydrologieTaeglich_hellinger3Clusters.rda')

proj=ProjectionBasedClustering::ICA(Trans4)
genmodellinear=GeneralizedUmatrix(Trans3,proj$ProjectedPoints,F,Cls=ClstTrue)  
plotTopographicMap(genmodellinear$Umatrix,genmodellinear$Bestmatches,Cls=ClstTrue)

setwd(ReDi("ExplainableAI4TimeSeries2020/07Dokumentationen"))
rgl.snapshot('HydrologyStructuresLinearhellinger3Clusters.png')

require(FCPS)
require(PPCI)
proj=PPCI::mcdr(Trans4,p = 2,maxit=1000,ftol=1e-7)

linear_cls=FCPS::ProjectionPursuitClustering(Trans4,ClusterNo = 3,Type="MaximumClusterbility")
plot(proj$fitted,col=ClstTrue)
plot(proj$fitted,col=linear_cls$Cls)
cc=ClstTrue;cc[cc>4]=4;
table(linear_cls$Cls,cc)

genmodellinear=GeneralizedUmatrix(Trans3,proj$fitted,F,Cls=ClstTrue)  
plotTopographicMap(genmodellinear$Umatrix,genmodellinear$Bestmatches,Cls=linear_cls$Cls,BmSize = 1.5)

plotTopographicMap(genmodellinear$Umatrix,genmodellinear$Bestmatches,Cls=ClstTrue,BmSize = 1.5)

setwd(ReDi("ExplainableAI4TimeSeries2020/07Dokumentationen"))
rgl.snapshot('HydrologyStructuresProjectionPursuit.png')

##2015 ----

library(PPCI)#cran
library(rgl)#cran
library(GeneralizedUmatrix)#cran
library(ProjectionBasedClustering)#cran
##nonlinear model proof ----
Disk=""
setwd(ReDi("ExplainableAI4TimeSeries2020/04DBS",Disk))

load('Hydro2015.rda')

require(FCPS)
require(PPCI)
proj=PPCI::mcdr(Trans,p = 2,maxit=1000,ftol=1e-7)

linear_cls=FCPS::ProjectionPursuitClustering(Trans,ClusterNo = 6,Type="MaximumClusterbility")
plot(proj$fitted,col=Cls)
plot(proj$fitted,col=linear_cls$Cls)

table(linear_cls$Cls,Cls)

genmodellinear=GeneralizedUmatrix(Trans,proj$fitted,F,Cls=Cls)  
TopviewTopographicMap(genmodellinear$Umatrix,genmodellinear$Bestmatches,Cls=linear_cls$Cls,BmSize = 1.5)

plotTopographicMap(genmodellinear$Umatrix,genmodellinear$Bestmatches,Cls=Cls,BmSize = 1.5)

setwd(ReDi("ExplainableAI4TimeSeries2020/07Dokumentationen"))
rgl.snapshot('HydrologyStructuresProjectionPursuit2015.png')


##2016 ----

library(PPCI)#cran
library(rgl)#cran
library(GeneralizedUmatrix)#cran
library(ProjectionBasedClustering)#cran
##nonlinear model proof ----
Disk=""
setwd(ReDi("ExplainableAI4TimeSeries2020/04DBS",Disk))

load('Hydro2016.rda')

require(FCPS)
require(PPCI)
proj=PPCI::mcdr(Trans2016,p = 2,maxit=1000,ftol=1e-7)

linear_cls=FCPS::ProjectionPursuitClustering(Trans2016,ClusterNo = 2,Type="MaximumClusterbility")
plot(proj$fitted,col=Cls)
plot(proj$fitted,col=linear_cls$Cls)

table(linear_cls$Cls,Cls2016)

genmodellinear=GeneralizedUmatrix(Trans2016,proj$fitted,F,Cls=Cls)  
TopviewTopographicMap(genmodellinear$Umatrix,genmodellinear$Bestmatches,Cls=linear_cls$Cls,BmSize = 1.5)

plotTopographicMap(genmodellinear$Umatrix,genmodellinear$Bestmatches,Cls=Cls2016,BmSize = 1.5)

setwd(ReDi("ExplainableAI4TimeSeries2020/07Dokumentationen"))
rgl.snapshot('HydrologyStructuresProjectionPursuit2016.png')