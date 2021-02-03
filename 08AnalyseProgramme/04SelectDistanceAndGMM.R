#04SelectDistanceAndGMM.R
library(DataVisualizations)#cran
library(parallelDist)#cran
library(FCPS) #cran, using  #ClusterDistances
setwd(ReDi('ExplainableAI4TimeSeries2020/01Transformierte'))
load(file='BackTransformationHydrologieAverageAgregation.rda')#,backtrafo,Trans3,Trans)

#Investigating Distance Distributions for multimodality

MDplot(ClusterDistances(as.matrix(parallelDist::parDist(Trans3,method = 'manhattan'))))#+xlab(names(FCPS)[i])+theme_bw()

MDplot(ClusterDistances(as.matrix(parallelDist::parDist(Trans3,method = 'maximum'))))#+xlab(names(FCPS)[i])+theme_bw()

MDplot(ClusterDistances(as.matrix(parallelDist::parDist(Trans3,method = 'minkowski',p = 3))))#+xlab(names(FCPS)[i])+theme_bw()

MDplot(ClusterDistances(DistanceMatrix(Trans3,method = 'pearson')))#+xlab(names(FCPS)[i])#+theme_bw()

#geht so
MDplot(ClusterDistances(as.matrix(parallelDist::parDist(Trans3,method = 'bhjattacharyya'),SampleSize=1e5)))#+xlab(names(FCPS)[i])+theme_bw()

MDplot(ClusterDistances(as.matrix(parallelDist::parDist(Trans3,method = 'bray'),SampleSize=1e5)))#+xlab(names(FCPS)[i])+theme_bw()

MDplot(ClusterDistances(as.matrix(parallelDist::parDist(Trans3,method = 'canberra'),SampleSize=1e5)))#+xlab(names(FCPS)[i])+theme_bw()

MDplot(ClusterDistances(as.matrix(parallelDist::parDist(Trans3,method = 'chord'))),SampleSize=1e5)#+xlab(names(FCPS)[i])+theme_bw()

MDplot(ClusterDistances(as.matrix(parallelDist::parDist(Trans3,method = 'divergence'))),SampleSize=1e5)#+xlab(names(FCPS)[i])+theme_bw()

MDplot(ClusterDistances(as.matrix(parallelDist::parDist(Trans3,method = 'geodesic'))),SampleSize=1e5)#+xlab(names(FCPS)[i])+theme_bw()

MDplot(ClusterDistances(as.matrix(parallelDist::parDist(Trans3,method = 'kullback'))),SampleSize=1e5)#+xlab(names(FCPS)[i])+theme_bw()

MDplot(ClusterDistances(as.matrix(parallelDist::parDist(Trans3,method = 'podani'))),SampleSize=1e5)#+xlab(names(FCPS)[i])+theme_bw()

MDplot(ClusterDistances(as.matrix(parallelDist::parDist(Trans3,method = 'soergel'))),SampleSize=1e5)#+xlab(names(FCPS)[i])+theme_bw()

MDplot(ClusterDistances(as.matrix(parallelDist::parDist(Trans3,method = 'wave'))),SampleSize=1e5)#+xlab(names(FCPS)[i])+theme_bw()

MDplot(ClusterDistances(as.matrix(parallelDist::parDist(Trans3,method = 'whittaker'))),Scaling = 'CompleteRobust',SampleSize=1e5)#+xlab(names(FCPS)[i])+theme_bw()

MDplot(ClusterDistances(as.matrix(parallelDist::parDist(Trans3,method = 'mahalanobis'))),SampleSize=1e5)#+xlab(names(FCPS)[i])+theme_bw()

MDplot(ClusterDistances(as.matrix(DistanceMatrix(Trans3,method = 'cosine'))),SampleSize=1e5)#+xlab(names(FCPS)[i])+theme_bw()

MDplot(ClusterDistances(as.matrix(DistanceMatrix(Trans3,method = 'sqEuclidean'))),SampleSize=1e5)#+xlab(names(FCPS)[i])+theme_bw()

MDplot(ClusterDistances(as.matrix(DistanceMatrix(Trans3))),SampleSize=1e5)#+xlab(names(FCPS)[i])+theme_bw()

#appropriate distance
MDplot(ClusterDistances(as.matrix(parallelDist::parDist(Trans3,method = 'hellinger'))),SampleSize=1e5)#+xlab(names(FCPS)[i])+theme_bw()

InputDistances=as.matrix(parallelDist::parDist(Trans3,method = 'hellinger'))

setwd(ReDi("ExplainableAI4TimeSeries2020/01Transformierte"))
save(file='hellinger_DistancesHydrologie.rda',InputDistances,Trans3,Trans,backtrafo,Time)
##
## gaussian mixuture model (gmm) if distances ----
library(diptest)#cran
library(AdaptGauss)#cran
library(DistributionOptimization)
setwd(ReDi("ExplainableAI4TimeSeries2020/01Transformierte"))
load(file='hellinger_DistancesHydrologie.rda')#,pvalueChitest,distvec,gmm,InputDistances,Trans3,Trans,backtrafo)

distvec=InputDistances[upper.tri(InputDistances,diag = F)]
dip=diptest::dip.test(distvec)
#algorithm so that gaussians overlap as less as possible
modelproposal=DistributionOptimization::DistributionOptimization(distvec,3)
gmm=list(Means=modelproposal$Means,SDs=modelproposal$SDs,Weights=modelproposal$Weights)
#interactive adjustments
gmm=AdaptGauss::AdaptGauss(distvec,gmm$Means,gmm$SDs,gmm$Weights)

dput(gmm)
gmm=list(Means = c(0.269584, 0.461048, 0.574664), 
         SDs = c(0.0846784, 0.060246, 0.1382311), 
         Weights = c(0.4353, 0.2355, 0.3277))


dec=AdaptGauss::BayesDecisionBoundaries(gmm$Means,gmm$SDs,gmm$Weights)

PlotMixtures(distvec,gmm$Means,gmm$SDs,gmm$Weights,lwd=5,xlab='Range of Distances',ylab='PDE',main='GMM of Distribution of Distances df',SingleGausses=T)
legend("topright",c("Distances df","GMM","Single Gaussians","Bayes Boundary"),col=c("black","red","blue","magenta"), lwd=2, lty=c(1,1,1,1),bty = "n")
abline(v = dec[1],col="magenta",lwd=5)

par(pty="s")
AdaptGauss::QQplotGMM(distvec,gmm$Means,gmm$SDs,gmm$Weights,ylab='Data = GMM of Distribution of Distances df')#,xlab='Gaussian Mixture Model (GMM)')

setwd(ReDi("ExplainableAI4TimeSeries2020/01Transformierte"))
save(file='hellinger_DistancesHydrologie.rda',distvec,gmm,InputDistances,Trans3,Trans,backtrafo,dip)
