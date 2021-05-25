# ImputingMetabolites
Compare R2 Values of Different Imputation Methods on Simulated Metabolite Data

## Sample R code

```{r}
####### Initial Variables
#Assortment of user defined variables
#Can change around depending on analysis needed
#These are the *default* settings that have been tested to work
num_of_metabos <- 20
sample_size <- 1000
num_to_remove <- 5
num_of_predictors <- 20

#Correlation matrices used to generate data
m_0 <- CreateBlockMatrix(num_of_metabos,.6,.6,1,10)
m_1 <- CreateBlockMatrix(num_of_metabos,.3,.4,0,5)
m_1[c(6:10,16:20),c(6:10,16:20)] <- 0
m_2 <- CreateBlockMatrix(num_of_metabos,0,.3,0,5)
m_2[c(1:5,11:15),c(1:5,11:15)] <- 0
m_3 <- CreateBlockMatrix(num_of_metabos,.0,.0,0,4)
factor1 <- seq(-1,1,.25)
#If following line is run, data is generated according to 1 factor model
factor2 <- 0
# #Uncomment following line to generate data according to 2 factor model
# factor2 <- seq(-1,1,.25)
factor3 <- 0
#Similarly factor 3 can be changed to non-zero values to generate data according to a 3 factor model
factor_mat <- expand.grid(factor1,factor2,factor3)

######Data generation, imputation, and analysis pipeline
# Creates observed and true data and observed and true correlation matrices
CorrelationAndData <- CreateCorrelation(m_0,m_1,m_2,m_3,
                                        num_of_metabos,
                                        sample_size,
                                        num_to_remove,
                                        factor_mat)

#Sorting output into right variable names
metabolites_array_true <- CorrelationAndData[[1]][[1]]
metabolites_array_obs <- CorrelationAndData[[1]][[2]]
corr_array_true <- CorrelationAndData[[2]][[1]]
corr_array_obs <- CorrelationAndData[[2]][[2]]

#Calculates average correlation (ignores heterogeneity)
corr_avg <- CalcCorrAvg(corr_array_obs,num_of_metabos)

#Preps data for MICE, outputs MICE correlation data and input for factor analysis
MICEArrays <- ReorganizeAndMICE(corr_array_obs,
                                factor_mat,
                                num_of_metabos,
                                num_of_predictors)

#This is un-fischer transformed
corr_array_mice <- MICEArrays[[1]]
#This is still fischer tranaformed, and will be untransformed in the next step in the pipeline (in FactorAnalysis())
reorganized_corr_mice_med <- MICEArrays[[2]]

#Factor analysis step. Runs 1-4 factor factor analysis.
#Outputs intermediaries (faRes*) to check loadings, scores,...
#Also outputs factor analysis correlation matrices in 3d array described in function description
#For example corr_array_fa3[,,10] is the 3 factor factor analysis array for study 10
factor_analysis_outputs <- FactorAnalysis(reorganized_corr_mice_med,num_of_metabos)
faRes1 <- factor_analysis_outputs[[1]][[1]]
faRes2 <- factor_analysis_outputs[[1]][[2]]
faRes3 <- factor_analysis_outputs[[1]][[3]]
faRes4 <- factor_analysis_outputs[[1]][[4]]
corr_array_fa1 <- factor_analysis_outputs[[2]][[1]]
corr_array_fa2 <- factor_analysis_outputs[[2]][[2]]
corr_array_fa3 <- factor_analysis_outputs[[2]][[3]]
corr_array_fa4 <- factor_analysis_outputs[[2]][[4]]

#Imputes missing metabolites
#A lot of imputs, correlation matrices for all tested methods
#output is a data frame which has imputed metabolite values by factor and method
imputed.df <- ImputeMissingMetabolites(sample_size,num_of_metabos,num_to_remove,factor_mat,
                                       metabolites_array_obs,metabolites_array_true,corr_avg,corr_array_true,
                                       corr_array_fa1,corr_array_fa2,corr_array_fa3,corr_array_fa4,corr_array_mice)

#Calculates R2 values from imputed.df
#Output is another data frame, easy to work with in ggplot
r2.df <- CalculateR2(imputed.df,factor_mat)

#Plots R2 values in line plot and boxplot
#First PlotR2OneFact produces a lineplot and should only be used to visualize when data is generated under 1 factor
if (length(factor2) == 1){
  PlotR2OneFact(r2.df)
}

#PlotR2MultFact produces boxplots for each method
PlotR2MultFact(r2.df)
```
