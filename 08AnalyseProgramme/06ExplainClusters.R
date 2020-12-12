#06ExplainClusters.R
library(DatabionicSwarm)
library(GeneralizedUmatrix)
library(ProjectionBasedClustering)
library(DataVisualizations)
library(rpart)
library(rpart.plot)
library(evtree)
library(FCPS)
library(ggplot2)
trainbestCART=function (Data, Names, Cls, MinClassSize = 10, PruningLevel = -1, 
          NrOfCrossvalidations = 10) 
{
  if (missing(Data)) {
    stop("Data not given")
  }
  if (missing(Cls)) {
    stop("Cls not given")
  }
  if (!is.matrix(Data)) {
    Data = as.matrix(Data)
    warning("Data is not a matrix")
  }
  if (missing(Names)) {
    warning("Names not given")
    Names = colnames(Data)
  }
  if (!is.vector(Cls)) {
    Cls = as.vector(Cls)
    warning("Cls is not a vector")
  }
  if (!is.vector(Names)) {
    Names = as.vector(Names)
    warning("Names is not a vector")
  }
  Data = as.matrix(Data)
  Cls = as.vector(Cls)
  Names = as.vector(Names)
  Data = cbind(Data, Cls)
  Names = append(Names, "Cls", after = length(Names))
  Data = data.frame(Data)
  names(Data) = make.names(Names)
  CartTree = rpart(formula = Cls ~ ., data = Data, method = "class", 
                   minsplit = MinClassSize, parms = list(split = "gini"), 
                   xval = NrOfCrossvalidations)
  if (PruningLevel == -1) {
    results = printcp(CartTree)
    results = data.frame(results)
    xerror = sort(na.last = NA, results$xerror, index.return = TRUE)
    newCP = results$CP[xerror$ix[1]]
    if (xerror$x[1] == xerror$x[2]) 
      newCP = results$CP[xerror$ix[2]]
    BestCartTree = prune(CartTree, cp = newCP)
  }
  else {
    return("PruningLevel has to be -1, otherwise its not implemented")
  }
  requireNamespace("rpart.plot")
  rpart.plot::prp(CartTree, main = "No. of incorrect classifications/No. of observations", 
                  type = 1, branch.lwd = 2, extra = 3, nn = TRUE, fallen.leaves = TRUE, 
                  branch = 0.5, faclen = 0, trace = 1, shadow.col = "gray", 
                  branch.lty = 3, split.cex = 1.2, split.prefix = "is ", 
                  split.suffix = "?", split.box.col = "lightgray", split.border.col = "darkgray", 
                  digits = 3)
  return(invisible(BestCartTree))
}
CART2Rules=function (Tree, digits) 
{
  requireNamespace("rpart.utils")
  if (!inherits(Tree, "rpart")) 
    stop("Not a legitimate \"rpart\" object")
  isleaf = Tree$frame$var == "<leaf>"
  idl = which(isleaf)
  output = vector(mode = "character", length = length(idl))
  NamesInRuleSet = vector(mode = "character")
  paths = rpart.utils::rpart.rules(Tree)
  rules = rpart.utils::rpart.lists(Tree)
  for (num in 1:length(idl)) {
    i = idl[num]
    leaf = as.numeric(row.names(Tree$frame[i, ]))
    class = Tree$frame$yval[i]
    output[num] = paste0("RuleNr ", class, ".", num, ":   Case is_a ", 
                         class, "$if  ")
    appendnext = ""
    path = paths[[leaf]]
    splitted = strsplit(x = (strsplit(path, ",")[[1]]), split = "")
    for (r in splitted) {
      output[num] = paste0(output[num], appendnext)
      number = as.numeric(r[2])
      direction = r[1]
      varname = attributes(rules[[direction]][number])$names
      if (!any(NamesInRuleSet == varname)) {
        NamesInRuleSet = append(NamesInRuleSet, varname)
      }
      varval = rules[[direction]][number][[1]]
      comparator = attributes(varval)$compare
      if (!missing(digits)) {
        varval = round(varval, digits)
      }
      output[num] = paste0(output[num], varname, " ", comparator, 
                           " ", varval)
      appendnext = "  and $  "
    }
    output[num] = paste0(output[num], " .$")
  }
  return(list(RuleSet = sort(output), NamesInRuleSet = NamesInRuleSet))
}

EvolutionTree <- function(Data, Cls){
  Cls <- factor(Cls)
  fullData = cbind(data.frame(Data) , Cls)
  fit <- evtree::evtree(Cls ~ ., fullData)
  plot(fit)
  return(fit)
}
#################################  ####################################################################################
# Knowledge Acquisition ----
########################################################################################################################
Disk="F"
setwd(ReDi("ExplainableAI4TimeSeries2020/04DBS",Disk))
load('HydrologieTaeglich_hellinger3Clusters.rda')

#ClstTrue=RenameDescendingClassSize(ClstTrue2)
plotTopographicMap(resUmatrix$Umatrix,resUmatrix$Bestmatches,ClstTrue,Imx = imx,BmSize = 1.1,NoLevels=10)
norm=NormalizeUmatrix(Trans4,resUmatrix$Umatrix,resUmatrix$Bestmatches)
plotTopographicMap(norm,resUmatrix$Bestmatches,ClstTrue,Imx = imx,BmSize = 1.1,NoLevels=10)

Header=colnames(Trans3)

#evolution tree
EvolutionTree(Trans3,ClstTrue)

#best cart
cart=trainbestCART(Trans3,colnames(Trans3),ClstTrue)

Clstsummed=ClstTrue
Clstsummed[Clstsummed>3]=4

ClassMDplot(Trans3[,1],Clstsummed,ClassNames = names,MinimalAmoutOfData = 20,PlotLegend=F,ColorSequence = GeneralizedUmatrix::DefaultColorSequence)$ggobject+theme_bw()+ggtitle('Class MDplot of Nitrate')+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+ylab('PDE')

ClassMDplot(Trans3[,7],Clstsummed,MinimalAmoutOfData = 20,PlotLegend=F,ColorSequence=GeneralizedUmatrix::DefaultColorSequence)$ggobject+theme_bw()+ggtitle('Class MDplot of Electric Conductivity')+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+ylab('PDE')


mback=Trans3
namesback=names(backtrafo$MinX)
for(i in 1:ncol(Trans3)){
  ind=which(namesback==colnames(Trans3)[i])
  if(length(ind)==1)
    mback[,i]=Trans3[,i]*backtrafo$Denom[ind]+backtrafo$MinX[ind]
}

range(mback[,1])
range(mback[,7])
mback[,'rain']=Trans[,"rain"]
Header=colnames(mback)
cart=trainbestCART(mback,Header,ClstTrue)
cc=ClstTrue
cc[cc>4]=4

plot(Trans$Time,cc,col=cc)


#cart
cart=trainbestCART(mback,Header,cc)
comments="06ExplainClusters.R;BackTransformed Data, of which Transormed data was used in DBS"

##For other XAIs
#procedure requires https://github.com/aultsch/DataIO
WriteLRN(FileName = "HydrologyToExplain.lrn",Data = mback,Header = colnames(mback),CommentOrDigits = comments)
#evolution tree
EvolutionTree(mback,cc)
ind=which(cc<4)
cart=trainbestCART(mback[ind,],Header,cc[ind])
#bloeder cart wechselt dann das features ohne sinn
colnames(mback)
cart=trainbestCART(mback[ind,-2],Header[-2],cc[ind])

##
require(RWeka) #java >8 needs to be installed
DF=data.frame(Cls=as.factor(as.character(ClstTrue)),mback)
DecisionRules=RWeka::JRip(Cls~.,data = DF)
#no relevant rules
summary(DecisionRules)

rules=CART2Rules(cart)
rules

names=c("Dry days with warm water of hgl",
        "Duality",
        "Dry days with cold water",
        "Outliers")
ClstTrue2=ClstTrue
ClstTrue2[ClstTrue2>3]=4
names=c('C1: DryDaysWarmWater','C2: Duality','C3: DryDaysColdWater','C4: Outliers')
ClassMDplot(mback[,1],ClstTrue2,MinimalAmoutOfData = 10,PlotLegend=F,ClassNames = names,ColorSequence = GeneralizedUmatrix::DefaultColorSequence)$ggobject+theme_bw()+ggtitle('Class MDplot of Nitrate')+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+ylab('PDE')

ClassMDplot(mback[,7],ClstTrue2,MinimalAmoutOfData = 10,PlotLegend=F,ClassNames = names,ColorSequence = GeneralizedUmatrix::DefaultColorSequence)$ggobject+theme_bw()+ggtitle('Class MDplot of Electric Conductivity')+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+ylab('PDE')