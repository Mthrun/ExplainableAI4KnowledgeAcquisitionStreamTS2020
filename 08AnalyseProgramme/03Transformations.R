#03Transformations.R

#required packages available on cran
library(TSAT)
library(DataVisualizations)#CRAN
library(ABCanalysis)#CRAN
library(DatabionicSwarm)#CRAN
library(parallelDist)#CRAN
KNNimputation=function (Data, k = 3){
#imputes data using k nearest neighbors
  dist2All=function (x, data, defined = rep(1, length(x)), distance = "euclidean") 
  {
    requireNamespace("parallelDist")
    xData <- rbind(x, data)
    xData <- xData[, which(defined == 1)]
    distToAll <- as.vector(parallelDist::parallelDist(xData, 
                                                      method = distance))[1:nrow(data)]
    return(distToAll = distToAll)
  }
  
  if (is.vector(Data)) 
    stop("Wrong dimension of Data.")
  AnzData = nrow(Data)
  d = ncol(Data)
  k = max(1, k)
  ImputedData = Data
  AnzNaNs = matrix(0,AnzData, d)
  dummy <- apply(Data, 1, function(x) !any(is.na(x)))
  HasNaNind = which(dummy == FALSE)
  PureInd = which(dummy == TRUE)
  AnzPure = length(PureInd)
  if (AnzPure < 1) 
    stop("KNNimputation: there is no case in Data which is complete")
  PureData = Data[PureInd, ]
  AnzNaNdata = length(HasNaNind)
  for (i in 1:AnzNaNdata) {
    n = HasNaNind[i]
    x = Data[n, ]
    NoNanColum = !(is.nan(x))
    DistToAll = t(dist2All(x, PureData, NoNanColum, "euclid"))
    Sorted <- sort(na.last = NA, DistToAll, index.return = TRUE)
    DistToAll = matrix(Sorted$x)
    Sind = matrix(Sorted$ix)
    KNNmean = colMeans(PureData[Sind[1:k], ])
    NaNsInd = which(is.nan(x), arr.ind = TRUE)
    ImputedData[n, NaNsInd] = KNNmean[NaNsInd]
  }
  return(ImputedData = ImputedData)
}
#load data form path
setwd(ReDi('ExplainableAI4TimeSeries2020/01Transformierte'))
V=TSAT::ReadDates('HydrologieAggregatedByMean2013bis2014.csv',ReDi('ExplainableAI4TimeSeries2020/09Originale'))
Time=V$Time
Trans=V$Data
Trans[!is.finite(Trans)]=NaN
DataVisualizations::Pixelmatrix(Trans)
Header=colnames(Trans)
#impute missing values specified as nan
Trans2=KNNimputation(Trans,7)
DataVisualizations::Pixelmatrix(Trans2)
Trans=Trans2
#distribution analysis
for(i in 1:14)
  DataVisualizations::InspectVariable(Trans[,i],Header[i])

#outliers def
OutlierInd=ABCanalysis::ABCanalysis(Trans[,14],PlotIt = T)$Aind

min(Trans[OutlierInd,14])
#capping
Trans[,14]=Trans[,14]/min(Trans[OutlierInd,14])

#standardization of features
i=12
DataVisualizations::InspectVariable(Trans[,i],Header[i])
Trans[,c(12,13)]=DataVisualizations::SignedLog(Trans[,c(12,13)])
DataVisualizations::InspectVariable(Trans[,i],Header[i])

ind=which( Trans[,14]>1)
Trans[ind,14]=1.1
DataVisualizations::InspectVariable(Trans[,14])

backtrafo=DatabionicSwarm::RobustNormalization(Trans[,c(1:13)],WithBackTransformation = T)

TransformedData=backtrafo$TransformedData
mback2=TransformedData
colnames(TransformedData)
#making sure backtransformation is possible
for(i in 1:ncol(TransformedData)){
  mback2[,i]=TransformedData[,i]*backtrafo$Denom[i]+backtrafo$MinX[i]
  if(i==11){
    print(backtrafo$Denom[i])
    print(backtrafo$MinX[i])
  }
}
range(mback2[,'St24'],na.rm = T)
range(Trans[,'St24'],na.rm = T)
sum(Trans[,c(1:13)]-mback2,na.rm = T)

Trans[,c(1:13)]=RobustNormalization(Trans[,c(1:13)],WithBackTransformation = F)
DataVisualizations::PlotMissingvalues(Trans,Names = colnames(Trans))
MDplot(Trans)

#correlations
Header=Header[1:14]

DataVisualizations::Pixelmatrix(cor(Trans),XNames = Header,YNames = Header)

pairs(Trans)
pairs(Trans[,c(5,6,9)])
pairs(Trans[,c(2,12,13)])
pairs(Trans[,c(10,11)])
pairs(Trans[,c(3,4)])

pairs(Trans[,c(11,9)])

x=cbind(Trans[,c(11,9)])
x[!is.finite(x[,1]),1]=median(x[,1],na.rm = T)
x[!is.finite(x[,2]),2]=median(x[,2],na.rm = T)
plot(x)
cor(x = Trans[,5],y = Trans[,6],method ="spearman")
cor.test(x = Trans[,5],y = Trans[,6],method ="spearman")
cor(x = Trans[,12],y = Trans[,13],method ="spearman")
cor.test(x = Trans[,12],y = Trans[,13],method ="spearman")


Header=Header[-c(5)]
Trans2=Trans[,-c(5)]

# pairs(Trans2)
# pairs(Trans2[,c(2,5,8)])
Trans3=Trans2#Trans2[,-c(12)]
Trans3=Trans2[,-c(12)]
Header=colnames(Trans3)

DataVisualizations::Pixelmatrix(cor(Trans3),XNames = Header,YNames = Header)

pairs(Trans3)

# cor(x = Trans[,'Wt18'],y = Trans[,'At47'],method ="spearman")
# x=(Trans[,'Wt18']-Trans[,'At47'])*2
# x=x-min(x)
# y=(Trans[,'Wt18']+Trans[,'At47'])/2
# Trans3[,'At47']=x
# Trans3[,'Wt18']=y

Pixelmatrix(cor(Trans3),XNames = Header,YNames = Header)
pairs(Trans3)
DataVisualizations::MDplot(Trans3)
# pairs(Trans3)
Trans3[Trans3<0]=0
Trans3[Trans3>1.1]=1
DataVisualizations::MDplot(Trans3)
#save transformed data
setwd(ReDi('ExplainableAI4TimeSeries2020/01Transformierte'))
save(file='BackTransformationHydrologieAverageAgregation.rda',backtrafo,Trans3,Trans,Time)
