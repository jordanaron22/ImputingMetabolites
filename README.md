# ImputingMetabolites
Compare R2 values of different imputation methods on simulated metabolite data. This repository mainly contains the r package described in the read me. Two additional folders, scripting and markdown, are also included. Scripting contains three files: an r script that will run this code on the hpc, a swarm file, and a folder named sample analysis. Inside sample analysis are two folders which contain data outputs from the script file (one 2factor, one 1factor) and a markdown that will produce plots from those files. The folder markdown contains past markdown documents that I've presented. Since creating them I've updated the general pipeline so those won't reproducible exactly, but can easily be reproduced using the script/package I've inlcuded. 

## R Package

### 1\. Installation 

```{r}
install.packages("devtools")
devtools::install_github("jordanaron22/ImputingMetabolites")
library(ImputingMetabolites)
```

### 2\. Example

#### 2.1 Initial User Defined Variables

```{r}
####### Initial Variables
#Assortment of user defined variables
#Can change around depending on analysis needed
#These are the *default* settings that have been tested to work
num_of_metabos <- 20
sample_size <- 1000
num_to_remove <- 4
#Set num_of_predictors 0 for most accurate MICE prediction
#Takes more time so I often set it to 20 for general debugging
#Determines amount of other values MICE considers during its imputation
#Set to 0 uses quickpred() function
num_of_predictors <- 20

```

##### 2.2 Define Number of Factors

Choose either 1 factor or 2 factors, not both.

###### 2.2.1 1 Factor

```{r}

###############################################
###Either run 1 Factor or 2 Factors, not both##
###############################################

#######1 Factor #######
low_val_1 <- 0.1
high_val_1 <- 0.4
low_val_2 <- 0
high_val_2 <- 0

#Correlation matrices used to generate data
m_0 <- CreateBlockMatrix(num_of_metabos,.6,.6,1,10)
m_1 <- CreateBlockMatrix(num_of_metabos,low_val_1,high_val_1,0,10)
m_2 <- CreateBlockMatrix(num_of_metabos,0,0,0,10)
m_3 <- CreateBlockMatrix(num_of_metabos,.0,.0,0,10)
#Setting factor2 and factor3 to 0 reduces the data generation to a 1 factor model
factor1 <- seq(-1,1,.25)
factor2 <- 0
factor3 <- 0
factor_mat <- expand.grid(factor1,factor2,factor3)
```

###### 2.2.1 2 Factors

```{r}
###############################################
###Either run 1 Factor or 2 Factors, not both##
###############################################

#######2 Factor #######
low_val_1 <- 0.1
high_val_1 <- 0.4
low_val_2 <- 0
high_val_2 <- 0.3

#Correlation matrices used to generate data
m_0 <- CreateBlockMatrix(num_of_metabos,.6,.6,1,10)
m_1 <- CreateBlockMatrix(num_of_metabos,low_val_1,high_val_1,0,5)
m_1[c(6:10,16:20),c(6:10,16:20)] <- 0
m_2 <- CreateBlockMatrix(num_of_metabos,low_val_2,high_val_2,0,5)
m_2[c(1:5,11:15),c(1:5,11:15)] <- 0
m_3 <- CreateBlockMatrix(num_of_metabos,.0,.0,0,4)
factor1 <- seq(-1,1,.25)
factor2 <- seq(-1,1,.25)
factor3 <- 0
#Similarly factor 3 can be changed to non-zero values to generate data according to a 3 factor model
factor_mat <- expand.grid(factor1,factor2,factor3)
#######
```

##### 2.3 Imputation and Analysis Pipeline

```{r}
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
```


##### 2.4 Visualization

```{r}
#Plots R2 values in line plot and boxplot
#First PlotR2OneFact produces a lineplot and should only be used to visualize when data is generated under 1 factor
if (length(factor2) == 1){
  PlotR2OneFact(r2.df)
}

#PlotR2MultFact produces boxplots for each method
PlotR2MultFact(r2.df)

```

### 3\. Scripting 

For scripting, I'm not sure how installing non CRAN packages work, so I would just create a script file with the functions. I read in the uniform distribution values for in the swarm file. To do this I include this block of code before everything else 

```{r}
args <- commandArgs(trailingOnly=T)
seed <- as.numeric(args[1])
set.seed(seed)
low_val_1 <- as.numeric(args[2])
high_val_1 <- as.numeric(args[3])
low_val_2 <- as.numeric(args[4])
high_val_2 <- as.numeric(args[5])
```

Then to save the files I include this block at the very end, which saves all of the r2.df outputs individually in a folder named after the inputs for the uniform distribution.

```{r}
block_vals <- paste0(low_val_1*10,high_val_1*10,low_val_2*10,high_val_2*10)
if (!dir.exists(block_vals)) {dir.create(block_vals)}
name <- paste0(block_vals,"/R24FactGrid",block_vals,"Seed",seed,".rda",sep="")
save(r2.df, file = name)
```

Finally, to generate the swarm file run 

```{r}
SwarmGenerator(500,"ComparingImpMethods.swarm","ComparingImpMethods.R",.3,.4,0,.2)
```

This will generate a swarm file named "ComparingImpMethods.swarm" that runs 500 instances of the R file named "ComparingImpMethods.R" with seeds 1-500. From the previous line, low_val_1 = .3, high_val_1=.4, low_val_2 = 0, and high_val_2 = .2. It will save r2.df from each instance into a folder named 3402.

Sample script can be found in the scripting folder.
