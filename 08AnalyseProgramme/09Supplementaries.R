#09Supplementariies
setwd(ReDi("ExplainableAI4TimeSeries2020/04DBS"))
load('HydrologieTaeglich_hellinger3Clusters.rda')
library(DataVisualizations)
library(FCPS)

## SI-A:Preprocessed Features ----
library(ggplot2)
MDplot(Trans3)+theme_bw()

#View(ContingencyTableSummary(RenameDescendingClassSize(ClstTrue[ind]),RenameDescendingClassSize(ClsVeri)))

## SI-C: Statistical Testing Taböe4 and Table 5 ----
ClassCount(ClstTrue)
ind1=which(ClstTrue==1)
ind2=which(ClstTrue==2)
ind3=which(ClstTrue==3)
ind4=which(ClstTrue==4)

Kl1=Trans3[ind1,7]
Kl2=Trans3[ind2,7]
Kl3=Trans3[ind3,7]
Kl4=Trans3[ind4,7]

ks.test(Kl1,Kl3)
ks.test(Kl1,Kl2)
ks.test(Kl2,Kl3)
ks.test(Kl2,Kl4)
ks.test(Kl3,Kl4)

length(Kl1)
length(Kl2)
length(Kl3)
length(Kl4)
length(Kl5)
InspectVariable(Kl4)
InspectVariable(Kl3)
InspectVariable(Kl2)
InspectVariable(Kl1)

Kl1=Trans3[ind1,1]
Kl2=Trans3[ind2,1]
Kl3=Trans3[ind3,1]
Kl4=Trans3[ind4,1]
length(Kl1)
length(Kl2)
length(Kl3)
length(Kl4)
length(Kl5)
InspectVariable(Kl4)
InspectVariable(Kl3)
InspectVariable(Kl2)
InspectVariable(Kl1)

ks.test(Kl1,Kl3)
ks.test(Kl1,Kl2)
ks.test(Kl2,Kl3)
ks.test(Kl2,Kl4)
ks.test(Kl3,Kl4)

kruskal.test(x = Trans3[,7],ClstTrue)

kruskal.test(list(Kl1,Kl2,Kl3,Kl4))

ks.test(alternative = "less",Kl1,Kl2,exact = T)

wilcox.test(alternative = "greater",Kl1,Kl2)

ks.test(alternative = "less",Kl1,Kl2,exact = T)


## SI-D: Sil
DataVisualizations::Silhouetteplot(Trans3,ClstTrue)

## SI-E: Distinction of Classes ----
setwd(ReDi('HydrologieSchwarmClustering2016/01Transformierte'))
V=TSAT::ReadDates('HydrologieAggregatedByMean2013bis2014.csv',ReDi('ExplainableAI4TimeSeries2020/09Originale'))
Time=as.Date(as.vector(as.matrix(V[,1])))
Trans=as.matrix(V[,2:15])
names=c('C1:DryDaysWarmWater','C2:Duality','C3:DryDaysColdWater')
cc=ClstTrue[ClstTrue<4]
ClassMDplot(Trans[ClstTrue==2,'rain'],cc,MinimalAmoutOfData = 50,ClassNames = names,PlotLegend=F)$ggobject+theme_bw()+ggtitle('Class MDplot of Rain')+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+ylab('PDE')

ClassMDplot(Trans[ClstTrue==2,'Wt18'],cc,MinimalAmoutOfData = 50,ClassNames = names,PlotLegend=F)$ggobject+theme_bw()+ggtitle('Class MDplot of Water Tempereture')+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+ylab('PDE')


out=ClassPDEplotMaxLikeli(Trans[ClstTrue<4,'Wt18'],cc,ClassNames = names[1:3],lwd=1.5)
out$ggobject+ggtitle('Class PDEplot of Water Temperature')+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+
  ylab('PDE')+xlab('Wt18 in ?C')+
  geom_vline(xintercept = 12.5, color = "red",size=1.5,linetype='dashed')

range(Trans[ClstTrue==3,'Wt18'])
ks.test(Trans[ClstTrue==1,'Wt18'],Trans[ClstTrue==2,'Wt18'],alternative = 'two.sided')


out=ClassPDEplotMaxLikeli(Trans[ClstTrue<3,'rain'],ClassNames = names[1:2],ClstTrue[ClstTrue<3],lwd=1.5)
out$ggobject+
  ggtitle('Class PDEplot of Rain')+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45,hjust = 1))+
  ylab('PDE')+xlab('Rain in mm/d')+xlim(c(0,10))+ylim(c(0,1))

ks.test(Trans[ClstTrue==1,'rain'],Trans[ClstTrue==2,'rain'],alternative = 'two.sided')

range(Trans[ClstTrue==3,'rain'])

range(Trans[ClstTrue==4,'St24'])
range(Trans[ClstTrue!=4,'St24'])

range(Trans[ClstTrue==4,'GWl3'])
range(Trans[ClstTrue!=4,'GWl3'])
range(Trans[ClstTrue==4,'Smoist24'])
range(Trans[ClstTrue!=4,'Smoist24'])

###
## Time Variations (unpublished) ----
setwd(ReDi("ExplainableAI4TimeSeries2020/04DBS"))
load('HydrologieTaeglich_hellinger3Clusters.rda')
cc=ClstTrue
cc[cc>4]=4
Trans=as.data.frame(Trans)
Trans$Time=Time
plot(Trans$Time,Trans$nnit13,col=cc)
requireNamespace("plotly")
TransNew=Trans
TransNew$Cls=as.factor(cc)
TransNew=TransNew[!is.na(TransNew$Time),]
TransNew=TransNew[order(TransNew$Time),]
p <- plotly::plot_ly(data=TransNew, x = ~Time, y = ~nnit13,color=~Cls,type = 'scatter',mode='lines+markers')
p


p <- plotly::plot_ly(data=TransNew)
p <- plotly::add_lines(p, x = ~Time, y = ~nnit13,color=~Cls, name = 'Nitrate in Siemens per m')
p <- plotly::add_lines(p, x = ~Time, y = ~Wt18,color=~as.factor(as.numeric(as.character(Cls))+5), name = 'Water temperature in ?C')
p <- plotly::layout(p, title = 'Temporal Variations', 
                    yaxis = list(side = "right", title = 'Nitrate in Siemens per m', showgrid = FALSE),
                    yaxis2 = list(overlaying = "y",side = "left", title = 'Water temperature in ?C', showgrid = FALSE),
                    xaxis = list(title = 'Time in days'))

p

TransNew=TransNew[as.numeric(as.character(TransNew$Cls))<4,]
X=TransNew$Time; Y1=TransNew$nnit13; Y2=TransNew$Wt13; xlab = "Time in days"; y1lab = 'Nitrate in Siemens per m'; y2lab = 'Water temperature in ?C'; 
main = "Temporal Variations of nitrate w.r.t temperature"; cols = c("black", "blue")

DualaxisClassplot(TransNew$Time,Y1=TransNew$nnit13,
                  Y2=TransNew$Wt13,Cls1 = as.numeric(TransNew$Cls),
                  Cls2 = as.numeric(TransNew$Cls),Colors =c('magenta','yellow','black','cyan','red','green'),
                  xlab = "Time in days", y1lab = 'Nitrate in Siemens per m', y2lab = 'Water temperature in ?C', 
                  main = "Temporal Variations of nitrate w.r.t temperature")

Classplot(TransNew$Time,Y=TransNew$con47,
          Cls = as.numeric(TransNew$Cls),
          xlab = "Time in days", ylab = 'C in Siemens per m', 
          main = "Temporal Variations of C")
#vielleicht ist besser die selben farben zu nutzen, aber statdessen die punkt form anzupassen?
DualaxisClassplot(TransNew$Time,Y1=TransNew$con47,
                  Y2=TransNew$Wt13,Cls1 = as.numeric(TransNew$Cls),
                  Cls2 = as.numeric(TransNew$Cls),Colors =c('magenta','yellow','black','cyan','red','green'),
                  xlab = "Time in days", y1lab = 'Electic Conductivity in Siemens per m', y2lab = 'Water temperature in ?C', 
                  main = "Temporal Variations of nitrate w.r.t temperature",Showgrid = F)

# p <- plotly::plot_ly(type='scatter',mode='markers',colors=c('magenta','yellow','black','cyan','red','green'))
# p <- plotly::add_trace(p, x = ~X, y = ~Y1,color=~as.factor(TransNew$Cls), name = y1lab)
# p <- plotly::add_trace(p, x = ~X, y = ~Y2,color=~as.factor(as.numeric(as.character(TransNew$Cls))+10), name = y2lab, 
#                        yaxis = "y2")
# p <- plotly::layout(p, title = main, yaxis2 = list(overlaying = "y", 
#                                                    side = "right", title = y2lab, showgrid = FALSE), xaxis = list(title = xlab, 
#                                                                                                                   showgrid = FALSE), yaxis = list(title = y1lab, showgrid = FALSE))
# p

p <- plotly::add_trace(p, x = ~Trans$Time, y = ~Trans$nnit13,color~as.factor(cc), name = 'asf')# line = list(color = cc))

p
p <- plotly::layout(p, title = main, xaxis = list(title = xlab, 
                                                  showgrid = FALSE), yaxis = list(title = ylab, showgrid = FALSE))
p

####




