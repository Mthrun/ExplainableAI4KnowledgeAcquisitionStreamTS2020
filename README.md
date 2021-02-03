[![DOI](https://zenodo.org/badge/250265216.svg)](https://zenodo.org/badge/latestdoi/250265216)

# Installation Guide


This is the installation guide for the framework of "Explainable AI Framework for Multivariate Hydrochemical Time Series."

## Table of contents

1. [Description](#description)
2. [Installation](#installation)
3. [General Remarks](#Generalremarks)
4. [References](#references)

## Description

An explainable AI (XAI) framework based on swarm intelligence is introduced in the following scripts 
and applied on the use case of multivariate time series in order to explain states of water bodies. 
The scripts require several packages to be installed.

## Installation

#### Installation using CRAN
All packages used for the framework are available on CRAN. Installation can be performed automatically with all dependencies via

```R
#general plotting
install.packages("ggplot2",dependencies = T)
install.packages("rgl",dependencies = T)

#packages used in scripts
install.packages("GeneralizedUmatrix",dependencies = T)
install.packages("DatabionicSwarm",dependencies = T)
install.packages("DataVisualizations",dependencies = T)
install.packages("AdaptGauss",dependencies = T)
install.packages("ABCanalysis",dependencies = T)
install.packages("diptest",dependencies = T)
install.packages("ProjectionBasedClustering",dependencies = T)
install.packages("rpart",dependencies = T)
install.packages("rpart",dependencies = T)
install.packages("rpart.plot",dependencies = T)
install.packages("evtree",dependencies = T)

#optional package for linear model
install.packages("PPCI",dependencies = T)

#dependencies are optional because there are many other clustering
#algorithms included
install.packages("FCPS",dependencies = T)

```

#### Installations using Github
There are two packages required if ones would like to reproduce the data aggregation step in "02AggregateData2Daily.R".
Otherwise one can start with "03Transformations.R" and use a simple read.table instead of TSAT::ReadDates.
This step would depend on the source of data and, hence, is optional.
Please note, that dependecies have to be installed manually.

```R
remotes::install_github('mthrun/TSAT')
remotes::install_github('aultsch/DataIO')
```

#### Installation using R Studio
Please note, that dependecies have to be installed manually.

*Tools -> Install Packages -> Repository (CRAN) -> Name of Package*

# General Remarks

The scripts in 08AnalyseProgramme are ordered consecutively. Figure 1 in [Thrun et al., 2020] provides an overview of the steps to be taken.
Raw data is stored in 90RawData, aggregated data in 09Originale and preprocessed data in 01Transformierte. Other folders are explained in the scripts.

## Additional information

| Authors website  | http://www.deepbionics.org/           						   |
| ---------------- |--------------------------------------------------------------:|
| License          | GPL-3                                 						   |
| Dependencies     | R (>= 3.5.0)                          						   |
| Bug reports      | https://github.com/Mthrun/ExplainableAI4TimeSeries2020/issues |


## References

1. Thrun, M. C., Ultsch, A., & Breuer, L.: Explainable AI Framework for Multivariate Hydrochemical Time Series, Machine Learning and Knowledge Extraction (MAKE), DOI: accepted, MDPI, 2021. 
2. [Thrun/Stier, 2020]  Thrun, M. C., & Stier, Q.: Fundamental Clustering Algorithms Suite, SoftwareX, Vol. 13(C), pp. 100642, DOI: 10.1016/j.softx.2020.100642, Elsevier, in press, 2021. 
3. [Thrun/Ultsch, 2020b]  Thrun, M. C., & Ultsch, A.: Swarm Intelligence for Self-Organized Clustering, Artificial Intelligence, Vol. 290, pp. 103237, DOI: 10.1016/j.artint.2020.103237, Elsevier, 2020. 
4. [Thrun/Ultsch, 2020c]  Thrun, M. C., & Ultsch, A. : Using Projection based Clustering to Find Distance and Density based Clusters in High-Dimensional Data, Journal of Classification, pp. 1-33, DOI 10.1007/s00357-020-09373-2, Springer, 2020. 

