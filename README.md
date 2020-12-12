[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/DatabionicSwarm)](https://cran.r-project.org/package=DatabionicSwarm)
[![DOI](https://zenodo.org/badge/113846715.svg)](https://zenodo.org/badge/latestdoi/113846715)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/grand-total/DatabionicSwarm?color=blue)](https://r-pkg.org/pkg/DatabionicSwarm)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/last-month/DatabionicSwarm?color=green)](https://r-pkg.org/pkg/DatabionicSwarm)

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

*Tools -> Install Packages -> Repository (CRAN) -> FCPS*

# General Remarks

The scripts in 08AnalyseProgramme are ordered consecutively. Figure 1 in [Thrun et al., 2020] provides an overview of the steps to be taken.
Raw data is stored in 90RawData, aggregated data in 09Originale and preprocessed data in 01Transformierte. Other folders are explained in the scripts.

## Additional information

| Authors website  | http://www.deepbionics.org/           |
| ---------------- |--------------------------------------:|
| License          | GPL-3                                 |
| Dependencies     | R (>= 3.5.0)                          |
| Bug reports      | https://github.com/Mthrun/FCPS/issues |


## References

1. [Thrun et al., 2020]  Thrun, M. C., Ultsch, A., & Breuer, L.: Explainable AI Framework for Multivariate Hydrochemical Time Series, Machine Learning and Knowledge Extraction (MAKE), Vol. under major revision, 2020.
2. [Thrun/Stier, 2020]  Thrun, M. C., & Stier, Q.: Fundamental Clustering Algorithms Suite SoftwareX, accepted, pp., 2020.
3. [Thrun/Ultsch, 2020b]  Thrun, M. C., & Ultsch, A.: Swarm Intelligence for Self-Organized Clustering, Journal of Artificial Intelligence, in press, DOI: 10.1016/j.artint.2020.103237, 2020.
4. [Thrun/Ultsch, 2020c]  Thrun, M. C., & Ultsch, A. : Using Projection based Clustering to Find Distance and Density based Clusters in High-Dimensional Data, Journal of Classification, DOI 10.1007/s00357-020-09373-2, accepted, Springer, 2020.

