# ImputingMetabolites
Compare R2 Values of Different Imputation Methods on Simulated Metabolite Data

## Sample R code

```{r}
#######
num_of_metabos <- 20
sample_size <- 1000
num_to_remove <- 4
num_of_predictors <- 20
m_0 <- CreateBlockMatrix(num_of_metabos,.6,.6,1,10)
m_1 <- CreateBlockMatrix(num_of_metabos,.3,.4,0,5)
m_1[c(6:10,16:20),c(6:10,16:20)] <- 0
m_2 <- CreateBlockMatrix(num_of_metabos,0,.3,0,5)
m_2[c(1:5,11:15),c(1:5,11:15)] <- 0
m_3 <- CreateBlockMatrix(num_of_metabos,.0,.0,0,4)
factor1 <- seq(-1,1,.25)
factor2 <- seq(-1,1,.25)
factor2 <- 0
factor3 <- 0
factor_mat <- expand.grid(factor1,factor2,factor3)
# #####

CorrelationAndData <- CreateCorrelation(m_0,m_1,m_2,m_3,
                                        num_of_metabos,
                                        sample_size,
                                        num_to_remove,
                                        factor_mat)

metabolites_array_true <- CorrelationAndData[[1]][[1]]
metabolites_array_obs <- CorrelationAndData[[1]][[2]]
corr_array_true <- CorrelationAndData[[2]][[1]]
corr_array_obs <- CorrelationAndData[[2]][[2]]

MICEArrays <- ReorganizeAndMICE(corr_array_obs,
                                factor_mat,
                                num_of_metabos,
                                num_of_predictors)

corr_array_mice <- MICEArrays[[1]]
reorganized_corr_mice_med <- MICEArrays[[2]]
factor_analysis_outputs <- FactorAnalysis(reorganized_corr_mice_med,num_of_metabos)
faRes1 <- factor_analysis_outputs[[1]][[1]]
faRes2 <- factor_analysis_outputs[[1]][[2]]
faRes3 <- factor_analysis_outputs[[1]][[3]]
faRes4 <- factor_analysis_outputs[[1]][[4]]
corr_array_fa1 <- factor_analysis_outputs[[2]][[1]]
corr_array_fa2 <- factor_analysis_outputs[[2]][[2]]
corr_array_fa3 <- factor_analysis_outputs[[2]][[3]]
corr_array_fa4 <- factor_analysis_outputs[[2]][[4]]

corr_avg <- CalcCorrAvg(corr_array_obs,num_of_metabos)

imputed.df <- ImputeMissingMetabolites(sample_size,num_of_metabos,num_to_remove,factor_mat,
                                       metabolites_array_obs,metabolites_array_true,corr_avg,corr_array_true,
                                       corr_array_fa1,corr_array_fa2,corr_array_fa3,corr_array_fa4,corr_array_mice)
r2.df <- CalculateR2(imputed.df,factor_mat)

if (length(factor2) == 1){
  PlotR2OneFact(r2.df)
}

PlotR2MultFact(r2.df)
```
