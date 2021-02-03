##02AggregateData2Daily.R
#install IO procedures
devtools::install_github('mthrun/TSAT')
devtools::install_github('aultsch/DataIO')
#ReDi is an internal path function, please set all paths manually
library(dbt.DataIO)
library(TSAT)
#load raw data
setwd(ReDi('ExplainableAI4TimeSeries2020/90RawData'))
V=dbt.DataIO::ReadLRN('HydrologieTime')
Key1=V$Key
df1=as.data.frame(V$Data)
V=dbt.DataIO::ReadLRN('HydrologieVars.lrn')
Key2=V$Key#
all(Key1==Key2)
Data=V$Data
Time=strptime(paste(df1$YYYY,df1$MM,df1$DD,df1$hh,df1$mm),format = '%Y %m %d %H %M')
#aggregate
fun=mean
for(i in 1:ncol(Data)){
  if(i==1)
    Completedata=TSAT::aggregateTime2Days(as.POSIXlt(Time),as.numeric(Data[,i]),FUN = fun,na.rm=TRUE,Header = colnames(Data)[i])
  else
    Completedata=merge(Completedata,TSAT::aggregateTime2Days(as.POSIXlt(Time),as.numeric(Data[,i]),FUN = fun,na.rm=TRUE,Header = colnames(Data)[i]),
               by.x='Time',by.y='Time',all.x = T,all.y = T)
}
#write out
TSAT::WriteDates('HydrologieAggregatedByMean2013bis2014.csv',Completedata,OutDirectory = ReDi('ExplainableAI4TimeSeries2020/09Originale'),Comments = '02AggregateData2Daily.R; alles mean',CleanNames = T)
