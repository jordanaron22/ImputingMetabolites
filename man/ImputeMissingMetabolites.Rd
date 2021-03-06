% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ComparingImpMethods.R
\name{ImputeMissingMetabolites}
\alias{ImputeMissingMetabolites}
\title{Imputes missing metabolites based on various correlation matrices}
\usage{
ImputeMissingMetabolites(
  sample_size = 1000,
  num_of_metabos = 20,
  num_to_remove = 5,
  factor_mat,
  metabolites_array_obs,
  metabolites_array_true,
  corr_avg,
  corr_array_true,
  corr_array_fa1,
  corr_array_fa2,
  corr_array_fa3,
  corr_array_fa4,
  corr_array_mice
)
}
\arguments{
\item{sample_size}{Number of individuals generated from correlation matrices}

\item{num_of_metabos}{Number of metabolites}

\item{num_to_remove}{Number of metabolites removed at a time from true data in CreateCorrelation()}

\item{factor_mat}{Matrix where 3 columns are three factors and rows are values of each factor}

\item{metabolites_array_obs}{Metabolite data with removed metabolites generated in CreateCorrelation()}

\item{metabolites_array_true}{Metabolite data generated in CreateCorrelation()}

\item{corr_avg}{Avgerage correlation. Ignores heterogeneity. Output from CalcCorrAvg()}

\item{corr_array_true}{True correlation matrix. An output from CreateCorrelation()}

\item{corr_array_fa1}{MICE and 1 factor factor analysis. An output from FactorAnalysis()}

\item{corr_array_fa2}{MICE and 2 factor factor analysis. An output from FactorAnalysis()}

\item{corr_array_fa3}{MICE and 3 factor factor analysis. An output from FactorAnalysis()}

\item{corr_array_fa4}{MICE and 4 factor factor analysis. An output from FactorAnalysis()}

\item{corr_array_mice}{Imputed with MICE. First output from ReorganizeAndMICE()}
}
\value{
data frame of correlation values by metabolite, factor, and imputation type
}
\description{
Imputes missing metabolites based on different correlation matrices calculated previously.
The different methods are: average correlation (ignoring heterogeneity), true correlation,
MICE imputation and 1-4 factor analysis, and MICE imputation. Outputs into dataframe. Takes time
(maybe a lot if large sample size) but useful for ggplot.
}
