#21GenerateFigures4NewData_AppendixKandJ.R
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
## 2015 ----
Disk="E"
setwd(ReDi("ExplainableAI4TimeSeries2020/04PBC",Disk))

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

deT2015=DecisionTree(Data[,c(13,3,6)],Cls = Cls,PlotIt=TRUE,CartPKG = F)
Rules2015=DecisionTree2Rules(deT2015,digits = 1)
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
setwd(ReDi("ExplainableAI4TimeSeries2020/04PBC",Disk))

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

deT2016=DecisionTree(Data2016[,c(13,3,6)],Cls = Cls2016,PlotIt=TRUE,CartPKG = F)

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

Rules2016=DecisionTree2Rules(deT2016,digits = 1)

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
