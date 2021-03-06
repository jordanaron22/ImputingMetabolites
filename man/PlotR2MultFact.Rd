% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ComparingImpMethods.R
\name{PlotR2MultFact}
\alias{PlotR2MultFact}
\title{Plots R2 of multiple factor imputation}
\usage{
PlotR2MultFact(r2.df)
}
\arguments{
\item{r2.df}{R2 dataframe, output from CalculateR2()}
}
\value{
Boxplot of R2 values by imputation type
}
\description{
Plots R2 in boxplot of different imputation methods.
Useful when underlying data is generated from multiple (2+) factors
}
