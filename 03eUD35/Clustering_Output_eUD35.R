
Clustering_Output_eUD35 <- function(file1, file2)
{
#V=Clustering_Output_eUD35(file1,file2)
# Reads the text output of eUD35
#
#INPUT
# file1		the clustering text file
# file2 	the patterns text file
#
#OPTIONAL
# none
#
#OUTPUT
# LIST V of
# Cls				[1:n] numerical vector of the classification for the n datapoints used in eUD3.5
# Patterns			[1:p] character vector of strings, each ith strings, i<p, defines a pattern
#
# Author: Hamza Tayyab, 10/2020
# Reference: none
# Packages required: none

  # Reading Clusters and patterns from separate files
  Clusters <- read.delim(file1)
  pattern <- read.delim(file2)
  
  pat <- pattern$Patterns
  Cls = Clusters$Cluster
  Cls=gsub("cluster","",Cls)
  Cls=as.numeric(Cls)
  # Initializing an empty list
  #patterns_list <- list()
  
  #patterns_list[[1]] <- pat[1]
  #for(p in 2:length(pat)){
  #  patterns_list[[p]] <- pat[p]
  #}
  out=list(Cls=Cls,Patterns=pat)
  return(out)
}

