#07OckhamsRazor.R
## Section 2.2.2
library(PPCI)#cran
library(rgl)#cran
library(GeneralizedUmatrix)#cran
library(ProjectionBasedClustering)#cran
##search for a simpler model ----
setwd(ReDi("ExplainableAI4TimeSeries2020/04PBC"))

load('HydrologieTaeglich_hellinger3Clusters.rda')
#use linear projection method
proj=ProjectionBasedClustering::ICA(Trans4)
genmodellinear=GeneralizedUmatrix(Trans3,proj$ProjectedPoints,F,Cls=ClstTrue)  

plotTopographicMap(genmodellinear$Umatrix,genmodellinear$Bestmatches,Cls=ClstTrue)

setwd(ReDi("ExplainableAI4TimeSeries2020/07Dokumentationen"))
rgl.snapshot('HydrologyStructuresLinearhellinger3Clusters.png')

require(FCPS)
require(PPCI)
#use linear projection method
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

#=> no linear model could be found
##2015 data----

library(PPCI)#cran
library(rgl)#cran
library(GeneralizedUmatrix)#cran
library(ProjectionBasedClustering)#cran
##nonlinear model proof ----
Disk=""
setwd(ReDi("ExplainableAI4TimeSeries2020/04PBC",Disk))

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


##2016 data ----

library(PPCI)#cran
library(rgl)#cran
library(GeneralizedUmatrix)#cran
library(ProjectionBasedClustering)#cran
##nonlinear model proof ----
Disk=""
setwd(ReDi("ExplainableAI4TimeSeries2020/04PBC",Disk))

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