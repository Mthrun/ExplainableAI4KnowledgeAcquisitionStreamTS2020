#20NeuesFramework.R
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
DecisionTree = function(Data,Cls,AutomaticPruning=TRUE,CartPKG="partykit",SplittingIndex="gini",...){
  # use crossvalidation and pruning to get the best fitting  decision tree
  #
  # INPUT
  # Data(1:d,1:n)        data Array of d cases with n variable
  # Cls(1:d)             vector of classes, d integer numbers, number k indicates class k
  # OPTINAL
  # AutomaticPruning     selects the optimal pruning based on a complexity parameter, overrides any manual setting of Complexity
  # MinClassSize         a subtree's class must have at least this number of members
  #                      10 if omitted
  # Complexity            defual: 0.01, value defining the level of pruning: Any split that does not decrease the overall lack of fit by a factor of cp is not attempted
  # PruningLevel          defualt 20, maximal tree depth 
  #
  # NrOfCrossvalidations The number of cross-validations for pruning, default NrOfCrossvalidations=10 (AutomaticPruning=TRUE) or NrOfCrossvalidations=0 (AutomaticPruning=TRUE)
  #                      10 if omitted
  # PlotIt               zeichnen des baumes   default  TRUE
  #
  # OUTPUT
  # BestCartTree         a binary tree where each non-terminal node is a
  #                      condition on a variable (clolumn of Data)
  #                      left son, means condition satisfied
  #author: MT 06/2020

  if(missing(Data)){
    stop('Data not given')
  }
  
  if(missing(Cls)){
    stop('Cls not given')
  }
  if(!is.matrix(Data)){
    Data=as.matrix(Data) 
    warning('Data is not a matrix')
  }

  if(!is.vector(Cls)){
    Cls=as.vector(Cls) 
    warning('Cls is not a vector')
  }  

  dots=list(...)
  
  if(is.null(dots[["Names"]])){
    Names=colnames(Data)
    if(is.null(Names)){
      Names=paste0('C',1:ncol(Data))
    }
  }else{
    Names=dots$Names
  }
  
  if(!is.vector(Names)){
    Names=as.vector(Names) 
    warning('Names is not a vector')
  }  
  
  if(is.null(dots[["PlotIt"]]))
    PlotIt=FALSE
  else
    PlotIt=dots$PlotIt
  
  if(is.null(dots[["MinClassSize"]]))
    MinClassSize=10
  else
    MinClassSize=dots$MinClassSize
  
  if(is.null(dots[["NrOfCrossvalidations"]])){
    if(isTRUE(AutomaticPruning))
      NrOfCrossvalidations=10
    else
      NrOfCrossvalidations=0
    
  }else{
    NrOfCrossvalidations=dots$NrOfCrossvalidations
  }
    
  
  if(is.null(dots[["Complexity"]]))
    Complexity=0.01 #default
  else
    Complexity=dots$Complexity
  
  if(is.null(dots[["PruningLevel"]]))
    PruningLevel=ncol(Data) #nicht hoeher als anzahl variablen
  else
    PruningLevel=dots$PruningLevel
  
  
  requireNamespace('rpart')

  if(is.null(dots[["PlotComplexity"]]))
    PlotComplexity=FALSE
  else
    PlotComplexity=dots$PlotComplexity
  
  if(is.null(dots[["Minbucket"]]))
    Minbucket=round(MinClassSize/3)
  else
    Minbucket=dots$Minbucket
  
  Data=cbind(Data,Cls) 
  Names=append(Names,"Cls", after = length(Names))
  Data=data.frame(Data) 
  names(Data)=make.names(Names) 

  #CART tree generation 
  if(CartPKG=="partykit"){
    requireNamespace("partykit")
    CartTree=partykit::ctree(formula=Cls~.,data=Data,maxdepth=PruningLevel,minbucket=Minbucket,minsplit=MinClassSize)#,control=list(maxdepth=PruningLevel,minbucket=Minbucket)) 
  }else{
    CartTree=rpart::rpart(formula=Cls~.,data=Data,method= "class",minsplit=MinClassSize,parms=list(split=SplittingIndex),xval=NrOfCrossvalidations,maxdepth=PruningLevel,minbucket=Minbucket) 
  

  #gets a table of optimal prunings based on a complexity parameter.
    if(isTRUE(AutomaticPruning)==TRUE&NrOfCrossvalidations>1){
      results=CartTree$cptable #capture.output
      results=data.frame(results)
      xerror=sort(na.last=NA,results$xerror,index.return = TRUE)
      Complexity=results$CP[xerror$ix[1]]
      if (xerror$x[1]==xerror$x[2]) Complexity=results$CP[xerror$ix[2]]
    }
  
    CartTree=rpart::prune(CartTree,cp=Complexity)
    if(isTRUE(PlotComplexity)){
      rpart::plotcp(CartTree)
    }
  }
    #requireRpackage('rpart.plot')
 
  if(isTRUE(PlotIt)){
    if(CartPKG=="partykit"){
      plot(CartTree)
    }else{
      requireNamespace('rpart.plot')
      rpart.plot::prp(CartTree, main="No. of incorrect classifications/No. of observations",
                      type=1,
                      branch.lwd=2,        # gyp for user manual
                      extra=3,             # Class models: misclassification rate at the node, expressed as the number of incorrect classifications and the number of observations in the node
                      nn=TRUE,             # display the node numbers
                      fallen.leaves=TRUE,  # put the leaves on the bottom of the page
                      branch=.5,           # change angle of branch lines
                      faclen=0,            # do not abbreviate factor levels
                      trace=1,             # print the automatically calculated cex
                      shadow.col="gray",   # shadows under the leaves
                      branch.lty=3,        # draw branches using dotted lines
                      split.cex=1.2,       # make the split text larger than the node text
                      split.prefix="is ",  # put "is " before split text
                      split.suffix="?",    # put "?" after split text
                      #col=cols, border.col=cols,   # show proportions greater than .5 in green
                      split.box.col="lightgray",   # lightgray split boxes (default is white)
                      split.border.col="darkgray", # darkgray border on split boxes
                      digits=2#,
                      #split.round=.001 # round the split box corners somewhat
      )
    }

  }
  return(invisible(CartTree))
}
DecisionTree2Rules <- function(Tree, digits){
  
  # RuleSet = DecisionTree2Rules(Tree,VariableNames);
  # The Rules derived from the Trees path from root to leaf
  # 
  # INPUT:
  # Tree                a binary tree where each non-terminal node is a
  #                     condition on a variable (clolumn of Data)
  #                     left son, means condition satisfied
  #
  # OUTPUT:
  # RuleSet        Array of Rules derived from Tree.
  # NamesInRuleSet      character array containing all the unique Names that are used in the RuleSet
  #                     May contain empty strings!
  #
 
  #requireRpackage('rpart.utils')
  requireNamespace('rpart.utils')
  if (inherits(Tree, "rpart")){ #stop("Not a legitimate \"rpart\" object")
  
  isleaf = Tree$frame$var == "<leaf>"
  idl = which(isleaf)
  output = vector(mode = 'character',length = length(idl))
  Rules2Fun = vector(mode = 'character',length = length(idl))
  #LIST2Fun = vector(mode = 'character',length = length(idl))
  LIST2Fun=c()
  NamesInRuleSet = vector(mode='character')
  paths = rpart.utils::rpart.rules(Tree)
  rules = rpart.utils::rpart.lists(Tree)
  RuleClsString=c()
  for(num in 1:length(idl)){
    LIST2Fun = c(LIST2Fun,paste0("Ind",num))
    
    i = idl[num] # position in the list
    leaf = as.numeric(row.names(Tree$frame[i,])) # position in tree
    class = Tree$frame$yval[i]  # class of the leaf
    Rules2Fun[num] = paste0("Ind",num," = which(")
    
    output[num] = paste0('RuleNo ', class,'.',num,':   Case is_a ',class,'$if  ')
    RuleClsString=c(RuleClsString,paste0('RuleCls[Ind',num,"] = ",class))
    appendnext = ""
    appendnext2 =""
    path = paths[[leaf]]        # path to reach the leaf
    
    splitted = strsplit(x = (strsplit(path,',')[[1]]),split = "")
    for(r in splitted){
      output[num] = paste0(output[num],appendnext)
      Rules2Fun[num] = paste0(Rules2Fun[num],appendnext2)
      
      number = as.numeric(r[2])
      direction = r[1]
      varname = attributes(rules[[direction]][number])$names
      if(!any(NamesInRuleSet == varname)){
        NamesInRuleSet = append(NamesInRuleSet, varname)
      }
      varval = rules[[direction]][number][[1]]
      comparator = attributes(varval)$compare
      if(!missing(digits)){
        varval = round(varval, digits)
      }
      output[num] = paste0(output[num],varname,' ',comparator,' ', varval)
      Rules2Fun[num] = paste0(Rules2Fun[num],varname,' ',comparator,' ', varval)
      appendnext="  and $  "
      appendnext2="  &  "
    }
    output[num] = paste0(output[num]," .$")
    Rules2Fun[num] = paste0(Rules2Fun[num],")")
  }
  ind=order(output)
  LIST2Fun=LIST2Fun[ind]
  return(list(RuleSet=output[ind],NamesInRuleSet=NamesInRuleSet,Rules2Fun=Rules2Fun[ind],LIST2Fun=LIST2Fun,RuleClsString=RuleClsString))
  }
  if (inherits(Tree, "party")){ #stop("Not a legitimate \"rpart\" object")
    #rules=partykit:::.list.rules.party(Tree)
    GetRules=function (x, i = NULL, ...) 
    {
      leafs <- partykit::nodeids(x, terminal = TRUE)
      Data_Leafs <- partykit::data_party(x, leafs)
      VarNamesL=lapply(Data_Leafs, names)
      VarNames= unique(unlist(VarNamesL))
      ClsLeafL=lapply(Data_Leafs,function(x,y) x[,which(names(x)=="Cls")])
      ClsLeafCountL=lapply(ClsLeafL, ClusterCount)
      MajorVoteL=lapply(ClsLeafCountL, function(x) {
        
        return(x$UniqueClusters[which.max(x$ClusterPercentages)])
      })
      MajorVote=unlist(MajorVoteL)
      names(MajorVote)=leafs
      
      NamesInRuleSet=VarNames[!VarNames %in% c("Cls","(fitted)","(weights)","(response)")]
      #NamesInRuleSet=NamesInRuleSet[-which(NamesInRuleSet=="Cls")]
      #NamesInRuleSet=NamesInRuleSet[-which(NamesInRuleSet=="(fitted)")]
      #NamesInRuleSet=NamesInRuleSet[-which(NamesInRuleSet=="(weights)")]
      #NamesInRuleSet=NamesInRuleSet[-which(NamesInRuleSet=="(response)")]
      .list.rules.party=function (x, i = NULL, ...) {
        if (is.null(i)) 
          i <- nodeids(x, terminal = TRUE)
        if (length(i) > 1) {
          ret <- sapply(i, .list.rules.party, x = x)
          names(ret) <- if (is.character(i)) 
            i
          else names(x)[i]
          return(ret)
        }
        if (is.character(i) && !is.null(names(x))) 
          i <- which(names(x) %in% i)
        stopifnot(length(i) == 1 & is.numeric(i))
        stopifnot(i <= length(x) & i >= 1)
        i <- as.integer(i)
        dat <- partykit::data_party(x, i)
        if (!is.null(x$fitted)) {
          findx <- which("(fitted)" == names(dat))[1]
          fit <- dat[, findx:ncol(dat), drop = FALSE]
          dat <- dat[, -(findx:ncol(dat)), drop = FALSE]
          if (ncol(dat) == 0) 
            dat <- x$data
        }
        else {
          fit <- NULL
          dat <- x$data
        }
        rule <- c()
        recFun <- function(node) {
          if (partykit::id_node(node) == i) 
            return(NULL)
          kid <- sapply(partykit::kids_node(node), partykit::id_node)
          whichkid <- max(which(kid <= i))
          split <- partykit::split_node(node)
          ivar <- partykit::varid_split(split)
          svar <- names(dat)[ivar]
          index <- partykit::index_split(split)
          if (is.factor(dat[, svar])) {
            if (is.null(index)) 
              index <- ((1:nlevels(dat[, svar])) > partykit::breaks_split(split)) + 
                1
            slevels <- levels(dat[, svar])[index == whichkid]
            srule <- paste(svar, " %in% c(\"", paste(slevels, 
                                                     collapse = "\", \"", sep = ""), "\")", sep = "")
          }
          else {
            if (is.null(index)) 
              index <- 1:length(kid)
            breaks <- cbind(c(-Inf, partykit::breaks_split(split)), c(partykit::breaks_split(split), 
                                                            Inf))
            sbreak <- breaks[index == whichkid, ]
            right <- partykit::right_split(split)
            srule <- c()
            if (is.finite(sbreak[1])) 
              srule <- c(srule, paste(svar, ifelse(right, ">", 
                                                   ">="), sbreak[1]))
            if (is.finite(sbreak[2])) 
              srule <- c(srule, paste(svar, ifelse(right, "<=", 
                                                   "<"), sbreak[2]))
            srule <- paste(srule, collapse = " & ")
          }
          rule <<- c(rule, srule)
          return(recFun(node[[whichkid]]))
        }#endrecFun
     
        node <- recFun(partykit::node_party(x))
        paste(rule, collapse = " & ")
      }#end.list.rules.party
      
      if (is.null(i)) 
        i <- partykit::nodeids(x, terminal = TRUE)
      if (length(i) > 1) {
        ret <- sapply(i, .list.rules.party, x = x)

        names(ret) <- if (is.character(i)) 
          i
        else names(x)[i]
  
        
        return(list(Rules=ret,NamesInRuleSet=NamesInRuleSet,MajorVote=MajorVote))
      }
    }

    rulesV=GetRules(Tree)
    RuleSet=rulesV$Rules
    Class=rulesV$MajorVote
    LIST2Fun=c()
    Rules2Fun = vector(mode = 'character',length = length(RuleSet))
    RuleClsString = vector(mode = 'character',length = length(RuleSet))
    Output = vector(mode = 'character',length = length(RuleSet))
    for(num in 1:length(RuleSet)){
      Output[num] = paste0('RuleNo ',Class[num],'.',num,':   Case is_a ',Class[num],'$if  ',RuleSet[num])
      
      LIST2Fun = c(LIST2Fun,paste0("Ind",num))
      Rules2Fun[num] = paste0("Ind",num," = which(")
      Rules2Fun[num] =paste0( Rules2Fun[num],paste(RuleSet[[num]],collapse = " & "),")")
      
      RuleClsString[num]=paste0('RuleCls[Ind',num,"] = ",Class[num])
    }
    
    return(list(RuleSet=RuleSet,NamesInRuleSet=rulesV$NamesInRuleSet,Rules2Fun=Rules2Fun,LIST2Fun=LIST2Fun,RuleClsString=RuleClsString,RuleDesc=Output))
    
  }
}

library(readxl)
library(TSAT)
Disk="E"
setwd(ReDi("Hydrologie2019/90RawData",Disk))
raw=read_excel("ts_synchro2015.xls",col_names = T,sheet = 1)
colnames(raw)

nnit13=TSAT::aggregateTime2Days(raw$indexzoo,as.numeric(raw$`N-nitrate`),FUN = mean,na.rm=TRUE,Header = "nnit13")
GWl3=TSAT::aggregateTime2Days(raw$indexzoo,as.numeric(raw$GWlevel3),FUN = mean,na.rm=TRUE,Header = "GWl3")
GWl25=TSAT::aggregateTime2Days(raw$indexzoo,as.numeric(raw$GWlevel25),FUN = mean,na.rm=TRUE,Header = "GWl25")
GWl32=TSAT::aggregateTime2Days(raw$indexzoo,as.numeric(raw$GWlevel32),FUN = mean,na.rm=TRUE,Header = "GWl32")

Wt13=TSAT::aggregateTime2Days(raw$indexzoo,as.numeric(raw$Tstream13),FUN = mean,na.rm=TRUE,Header = "Wt13")
Wt18=TSAT::aggregateTime2Days(raw$indexzoo,as.numeric(raw$Tstream18),FUN = mean,na.rm=TRUE,Header = "Wt18")
sol71=TSAT::aggregateTime2Days(raw$indexzoo,as.numeric(raw$SolarRadiation),FUN = mean,na.rm=TRUE,Header = "sol71")

con47=TSAT::aggregateTime2Days(raw$indexzoo,as.numeric(raw$Conductivity),FUN = mean,na.rm=TRUE,Header = "con47")
At47=TSAT::aggregateTime2Days(raw$indexzoo,as.numeric(raw$Tair),FUN = mean,na.rm=TRUE,Header = "At47")
Smoist24=TSAT::aggregateTime2Days(raw$indexzoo,as.numeric(raw$SoilMoisture),FUN = mean,na.rm=TRUE,Header = "Smoist24")
St24=TSAT::aggregateTime2Days(raw$indexzoo,as.numeric(raw$Tsoil),FUN = mean,na.rm=TRUE,Header = "St24")
q13=TSAT::aggregateTime2Days(raw$indexzoo,as.numeric(raw$Discharge13),FUN = mean,na.rm=TRUE,Header = "q13")
q18=TSAT::aggregateTime2Days(raw$indexzoo,as.numeric(raw$Discharge18),FUN = mean,na.rm=TRUE,Header = "q18")
rain=TSAT::aggregateTime2Days(raw$indexzoo,as.numeric(raw$RainIntensity),FUN = mean,na.rm=TRUE,Header = "rain")


DF=merge(nnit13,GWl3,by.x = "Time",by.y = "Time",all.x = T,all.y = T)
DF=merge(DF,GWl25,by.x = "Time",by.y = "Time",all.x = T,all.y = T)
DF=merge(DF,GWl32,by.x = "Time",by.y = "Time",all.x = T,all.y = T)

DF=merge(DF,Wt13,by.x = "Time",by.y = "Time",all.x = T,all.y = T)
DF=merge(DF,Wt18,by.x = "Time",by.y = "Time",all.x = T,all.y = T)
DF=merge(DF,sol71,by.x = "Time",by.y = "Time",all.x = T,all.y = T)

DF=merge(DF,con47,by.x = "Time",by.y = "Time",all.x = T,all.y = T)
DF=merge(DF,At47,by.x = "Time",by.y = "Time",all.x = T,all.y = T)
DF=merge(DF,Smoist24,by.x = "Time",by.y = "Time",all.x = T,all.y = T)
DF=merge(DF,q13,by.x = "Time",by.y = "Time",all.x = T,all.y = T)
DF=merge(DF,q18,by.x = "Time",by.y = "Time",all.x = T,all.y = T)
DF=merge(DF,rain,by.x = "Time",by.y = "Time",all.x = T,all.y = T)


DF[which(is.na(DF),arr.ind = T)]=NaN


WriteDates("20201226Hydrologie2015",TSdata = DF,OutDirectory =ReDi("Hydrologie2019/09Originale",Disk) )

V=ReadDates("20201226Hydrologie2015",ReDi("Hydrologie2019/09Originale",Disk))
Data=V$Data
DataVisualizations::PlotMissingvalues(DF)

Data=as.matrix(DF[,c(2:14)])
Data=Data[is.finite(Data[,"nnit13"]),]
Data=KNNimputation(Data,7)
DataVisualizations::PlotMissingvalues(Data)

plot(Data[,"q13"],Data[,"q18"])
plot(Data[,"Wt13"],Data[,"Wt18"])
which(colnames(Data)=="q13")
which(colnames(Data)=="Wt13")
Data2=Data[,-c(5,11)]

Trans=Data2*NaN
Trans[,"q18"]=DataVisualizations::SignedLog(Data2[,"q18"])

which(colnames(Data)=="q18")
which(colnames(Data)=="rain")

Trans[,c(1:11)]=DatabionicSwarm::RobustNormalization(Data2[,c(1:11)],Capped = T)

boxplot(Data2[,"rain"])
ABCanalysis::ABCanalysis(Data2[,"rain"])$ABLimit

Trans[,"rain"]=Data2[,"rain"]/ABCanalysis::ABCanalysis(Data2[,"rain"])$ABLimit

Trans[Trans[,"rain"]>1,"rain"]=1.1


MDplot(Trans)

out=Distances::DistanceDistributions(Trans)
which(!is.finite(Trans))


Distance=as.matrix(parallelDist::parDist(x = Trans,method = 'minkowski',p=0.1))
Distvec=Distance[upper.tri(Distance,diag = F)]
library(AdaptGauss)
library(DataVisualizations)

gmm=AdaptGauss(Distvec)
QQplotGMM(Distvec,gmm$Means,gmm$SDs,gmm$Weights)
MDplot(RobustNormalization(Distance[upper.tri(Distance,diag = F)]))

Proj=Pswarm(Distance)
genU=GeneratePswarmVisualization(Data = Trans,Proj$ProjectedPoints,Proj$LC,PlotIt = T)
#normU=NormalizeUmatrix(Trans,genU$Umatrix,genU$Bestmatches)
Cls=interactiveClustering(genU$Umatrix,genU$Bestmatches)

ClusterCount(Cls)
DecisionTree=DecisionTree(Data,Cls = Cls,PlotIt=TRUE,CartPKG = F)

colnames(Data)
deT=DecisionTree(Data[,c(13,3,6)],Cls = Cls,PlotIt=TRUE,CartPKG = F)

Rules=DecisionTree2Rules(DecisionTree,digits = 1)

ClusterApply(Data[,"nnit13"],Cls,FUN = mean)

ClassMDplot(Data[,"nnit13"],Cls,MinimalAmoutOfData = 100,Ordering = "Average")
ClassMDplot(Data[,"con47"],Cls,MinimalAmoutOfData = 100,Ordering = "Average")


imx2015=ProjectionBasedClustering::interactiveGeneralizedUmatrixIsland(genU$Umatrix,genU$Bestmatches,Cls = Cls)

DistanceMethod="minkowskiP=0.1"
BayesBoundaries2015=AdaptGauss::BayesDecisionBoundaries(gmm$Means,gmm$SDs,gmm$Weights)
intradist2015=FCPS::ClusterIntraDistances(Distance,Cls = Cls,PlotIt = F)
MDplot(intradist2015)


ScriptPfad="/Hydrologie2019/08Analyse/20NeuesFramework.R"
setwd(ReDi("ExplainableAI4TimeSeries2020/04DBS",Disk))
#E:\Subversion\PRO\Research\ExplainableAI4TimeSeries2020\04PBC
save(file="Hydro2015.rda",Data,DF,Trans,DecisionTree,Rules,Distance,gmm,Proj,genU,Cls,ScriptPfad,imx2015,BayesBoundaries2015,intradist2015,DistanceMethod)

TopviewTopographicMap(genU$Umatrix,genU$Bestmatches,Cls = Cls,Imx = imx2015)
#save.image("E:/Subversion/PRO/Research/Hydrologie2019/01Transformierte/20202312Hydro2015XAI.RData")

Heatmap(Distance,Cls)

#2016 ----
library(readxl)
Disk="E"
setwd(ReDi("Hydrologie2019/90RawData",Disk))
raw=read_excel("ts_synchro2016.xls",col_names = T,sheet = 1)
colnames(raw)

nnit13=TSAT::aggregateTime2Days(raw$indexzoo,as.numeric(raw$`N-nitrate`),FUN = mean,na.rm=TRUE,Header = "nnit13")
GWl3=TSAT::aggregateTime2Days(raw$indexzoo,as.numeric(raw$GWlevel3),FUN = mean,na.rm=TRUE,Header = "GWl3")
GWl25=TSAT::aggregateTime2Days(raw$indexzoo,as.numeric(raw$GWlevel25),FUN = mean,na.rm=TRUE,Header = "GWl25")
GWl32=TSAT::aggregateTime2Days(raw$indexzoo,as.numeric(raw$GWlevel32),FUN = mean,na.rm=TRUE,Header = "GWl32")

Wt13=TSAT::aggregateTime2Days(raw$indexzoo,as.numeric(raw$Tstream13),FUN = mean,na.rm=TRUE,Header = "Wt13")
Wt18=TSAT::aggregateTime2Days(raw$indexzoo,as.numeric(raw$Tstream18),FUN = mean,na.rm=TRUE,Header = "Wt18")
sol71=TSAT::aggregateTime2Days(raw$indexzoo,as.numeric(raw$SolarRadiation),FUN = mean,na.rm=TRUE,Header = "sol71")

con47=TSAT::aggregateTime2Days(raw$indexzoo,as.numeric(raw$Conductivity),FUN = mean,na.rm=TRUE,Header = "con47")
At47=TSAT::aggregateTime2Days(raw$indexzoo,as.numeric(raw$Tair),FUN = mean,na.rm=TRUE,Header = "At47")
Smoist24=TSAT::aggregateTime2Days(raw$indexzoo,as.numeric(raw$SoilMoisture),FUN = mean,na.rm=TRUE,Header = "Smoist24")
St24=TSAT::aggregateTime2Days(raw$indexzoo,as.numeric(raw$Tsoil),FUN = mean,na.rm=TRUE,Header = "St24")
q13=TSAT::aggregateTime2Days(raw$indexzoo,as.numeric(raw$Discharge13),FUN = mean,na.rm=TRUE,Header = "q13")
q18=TSAT::aggregateTime2Days(raw$indexzoo,as.numeric(raw$Discharge18),FUN = mean,na.rm=TRUE,Header = "q18")
rain=TSAT::aggregateTime2Days(raw$indexzoo,as.numeric(raw$RainIntensity),FUN = mean,na.rm=TRUE,Header = "rain")


DF2016=merge(nnit13,GWl3,by.x = "Time",by.y = "Time",all.x = T,all.y = T)
DF2016=merge(DF2016,GWl25,by.x = "Time",by.y = "Time",all.x = T,all.y = T)
DF2016=merge(DF2016,GWl32,by.x = "Time",by.y = "Time",all.x = T,all.y = T)

DF2016=merge(DF2016,Wt13,by.x = "Time",by.y = "Time",all.x = T,all.y = T)
DF2016=merge(DF2016,Wt18,by.x = "Time",by.y = "Time",all.x = T,all.y = T)
DF2016=merge(DF2016,sol71,by.x = "Time",by.y = "Time",all.x = T,all.y = T)

DF2016=merge(DF2016,con47,by.x = "Time",by.y = "Time",all.x = T,all.y = T)
DF2016=merge(DF2016,At47,by.x = "Time",by.y = "Time",all.x = T,all.y = T)
DF2016=merge(DF2016,Smoist24,by.x = "Time",by.y = "Time",all.x = T,all.y = T)
DF2016=merge(DF2016,q13,by.x = "Time",by.y = "Time",all.x = T,all.y = T)
DF2016=merge(DF2016,q18,by.x = "Time",by.y = "Time",all.x = T,all.y = T)
DF2016=merge(DF2016,rain,by.x = "Time",by.y = "Time",all.x = T,all.y = T)

DF2016[which(is.na(DF2016),arr.ind = T)]=NaN

TSAT::WriteDates("20201226Hydrologie2016",TSdata = DF2016,OutDirectory =ReDi("Hydrologie2019/09Originale",Disk) )

V=TSAT::ReadDates("20201226Hydrologie2016",ReDi("Hydrologie2019/09Originale",Disk))
Data2016=V$Data

DataVisualizations::PlotMissingvalues(Data2016)

#Data2016=as.matrix(DF2016[,c(2:14)])
Data2016=Data2016[is.finite(Data2016[,"nnit13"]),]
Data2016=KNNimputation(Data2016,7)
DataVisualizations::PlotMissingvalues(Data2016)

plot(Data2016[,"q13"],Data2016[,"q18"])
plot(Data2016[,"Wt13"],Data2016[,"Wt18"])
which(colnames(Data2016)=="q18")
which(colnames(Data2016)=="Wt13")
Data20162=Data2016[,-c(5,12)]
colnames(Data20162)[10]="q18"
Trans2016=Data20162*NaN
Trans2016[,"q18"]=DataVisualizations::SignedLog(Data20162[,"q18"])

which(colnames(Data2016)=="q18")
which(colnames(Data2016)=="rain")

Trans2016[,c(1:11)]=DatabionicSwarm::RobustNormalization(Data20162[,c(1:11)],Capped = T)

boxplot(Data20162[,"rain"])
ABCanalysis::ABCanalysis(Data20162[,"rain"])$ABLimit

Trans2016[,"rain"]=Data20162[,"rain"]/ABCanalysis(Data20162[,"rain"])$ABLimit

Trans2016[Trans2016[,"rain"]>1,"rain"]=1.1

out=Distances::DistanceDistributions(Trans2016)
which(!is.finite(Trans))
out$ggobject

Distance=as.matrix(parallelDist::parDist(x = Trans2016,method = 'minkowski',p=0.5))
Distance=as.matrix(parallelDist::parDist(x = Trans2016,method = 'hellinger'))
MDplot(DatabionicSwarm::RobustNormalization(Distance[upper.tri(Distance,diag = F)]))

outFULL=DatabionicSwarm::Pswarm(Distance)
genUFULL=DatabionicSwarm::GeneratePswarmVisualization(Data = Trans2016,outFULL$ProjectedPoints,outFULL$LC,PlotIt = T)
Cls2016=ProjectionBasedClustering::interactiveClustering(genUFULL$Umatrix,genUFULL$Bestmatches)


ClusterCount(Cls2016)
DecisionTree=DecisionTree(Data2016,Cls = Cls2016,PlotIt=TRUE,CartPKG = F)

colnames(Data2016)
deT=DecisionTree(Data2016[,c(13,3,6)],Cls = Cls2016,PlotIt=TRUE,CartPKG = F)

Rules=DecisionTree2Rules(DecisionTree,digits = 1)

ClusterApply(Data2016[,"nnit13"],Cls2016,FUN = mean)

ClassMDplot(Data2016[,"nnit13"],Cls2016,MinimalAmoutOfData = 100,Ordering = "Average")
ClassMDplot(Data2016[,"con47"],Cls2016,MinimalAmoutOfData = 100,Ordering = "Average")


outModel=ProjectionBasedClustering::IPBC(Trans2016)
Cls2016=outModel$Cls
Cls2016[Cls2016==3]=2
ClassMDplot(Data2016[,"nnit13"],Cls2016,MinimalAmoutOfData = 100,Ordering = "Average")
ClassMDplot(Data2016[,"con47"],Cls2016,MinimalAmoutOfData = 100,Ordering = "Average")

DecisionTree=DecisionTree(Data2016[,c(13,3,6)],Cls = Cls2016,PlotIt=TRUE,CartPKG = T)
Heatmap(Trans2016,Cls2016)

ScriptPfad="/Hydrologie2019/08Analyse/20NeuesFramework.R"
Distance=as.matrix(parallelDist::parDist(x = Trans2016))

Distvec=Distance[upper.tri(Distance,diag = F)]
gmm=AdaptGauss(Distvec)
setwd(ReDi("ExplainableAI4TimeSeries2020/04DBS",Disk))

imx2016=ProjectionBasedClustering::interactiveGeneralizedUmatrixIsland(outModel$Umatrix,outModel$Bestmatches,Cls = Cls2016)

DistanceMethod="Euclidean"
BayesBoundaries2016=AdaptGauss::BayesDecisionBoundaries(gmm$Means,gmm$SDs,gmm$Weights)
FCPS::ClusterIntraDistances(Distance,Cls = Cls2016,PlotIt = T)
intradist2016=FCPS::ClusterIntraDistances(Distance,Cls = Cls2016,PlotIt = F)

#E:\Subversion\PRO\Research\ExplainableAI4TimeSeries2020\04PBC
save(file="Hydro2016.rda",Data2016,DF2016,Trans2016,DecisionTree,Rules,Distance,gmm,outModel,Cls2016,ScriptPfad,imx2016,DistanceMethod,BayesBoundaries2016,intradist2016)