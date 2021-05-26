#' Creates Autoregressive Correlation Matrix
#'
#' Creates autoregressive correlation matrix of variable size with variable rho and diagonal values
#' @param num_of_metabos Number of metabolites (dimension of matrix)
#' @param rho Degree of autoregressiveness
#' @param diag_val Diagonal value. Set to 1 if m_0, 0 else
#' @return Autoregressive square matrix of dimension num_of_metabo
#' @export
CreateAutoRegressiveMat <- function(num_of_metabos, rho,diag_val) {
  exponent <- abs(matrix(1:num_of_metabos - 1, nrow = num_of_metabos, ncol = num_of_metabos, byrow = TRUE) - (1:num_of_metabos - 1))
  ar1 <- rho^exponent
  diag(ar1) <- diag_val
  return(ar1)
}

ImputeMissingMu <- function(metabolites,corr_mat,set_of_missing,set_of_obs,sample_size){
  missing_obs_cov <- corr_mat[set_of_missing,set_of_obs]
  inverted_obs_var <- solve(corr_mat[set_of_obs,set_of_obs],tol = 1e-40)
  for (ind in 1:sample_size){
    adj_mean <- matrix(c(metabolites[ind,set_of_obs]))
    mu_estimated <- ((missing_obs_cov %*% inverted_obs_var) %*% adj_mean)
    for (est_ind in 1:length(set_of_missing)){
      metabolites[ind,set_of_missing[est_ind]] <- mu_estimated[est_ind]
    }
  }
  return(metabolites)
}

#' Creates Block Diagonal Correlation Matrix
#'
#' Creates block diagonal correlation matrix of variable size with blocks and diagonal values.
#' Correlation values are drawn from a random uniform distribution between low_val and high_val
#' @param num_of_metabos Number of metabolites (dimension of matrix)
#' @param low_val lower bound for block correlation. Correlation is uniform between low_val and high_val
#' @param high_val upper bound for block correlation. Correlation is uniform between low_val and high_val
#' @param diag_val Diagonal value. Set to 1 if m_0, 0 else
#' @param size_of_block Block size. Default is 10. num_of_metabos must be divisible by size_of_block
#' @return Block diagonal matrix of dimension num_of_metabo
#' @import stats
#' @import MASS
#' @import ggplot2
#' @export
CreateBlockMatrix <- function(num_of_metabos,low_val,high_val,diag_val,size_of_block = 10){
  block <- matrix(runif(size_of_block^2,low_val,high_val),size_of_block,size_of_block)
  block[lower.tri(block)] <- 0
  block <- block + t(block)
  diag(block) <- diag_val
  return(diag(num_of_metabos/size_of_block) %x% block)
}

CreateMissingAndObsList <- function(fact_len,num_to_remove,num_of_metabos){
  missing_list <- list()
  obs_list <- list()
  starting_metabo <- floor((num_of_metabos - (num_to_remove*4))/2) + 1
  random_metabos <- c(starting_metabo:num_of_metabos)

  for (i in 1:fact_len){
    missing_list <- list.append(missing_list,random_metabos[1:num_to_remove])
    obs_list <- list.append(obs_list,c(1:num_of_metabos)[!(c(1:num_of_metabos) %in% random_metabos[1:num_to_remove])])
  }

  for (i in 1:fact_len){
    missing_list <- list.append(missing_list,random_metabos[(num_to_remove+1):(2*num_to_remove)])
    obs_list <- list.append(obs_list,c(1:num_of_metabos)[!(c(1:num_of_metabos) %in% random_metabos[(num_to_remove+1):(2*num_to_remove)])])
  }

  for (i in 1:fact_len){
    missing_list <- list.append(missing_list,random_metabos[(2*num_to_remove+1):(3*num_to_remove)])
    obs_list <- list.append(obs_list,c(1:num_of_metabos)[!(c(1:num_of_metabos) %in% random_metabos[(2*num_to_remove+1):(3*num_to_remove)])])
  }

  for (i in 1:fact_len){
    missing_list <- list.append(missing_list,random_metabos[(3*num_to_remove+1):(4*num_to_remove)])
    obs_list <- list.append(obs_list,c(1:num_of_metabos)[!(c(1:num_of_metabos) %in% random_metabos[(3*num_to_remove+1):(4*num_to_remove)])])
  }

  return(list(missing_list,obs_list))
}

Vec2Cor <- function(vector,num_of_metabos){
  corr_mat <- matrix(0,num_of_metabos,num_of_metabos)
  corr_mat[upper.tri(corr_mat, diag=FALSE)] <- vector
  corr_mat <- corr_mat + t(corr_mat)
  diag(corr_mat) <- 1
  return(corr_mat)
}

#' Plots R2 of multiple factor imputation
#'
#' Plots R2 in boxplot of different imputation methods.
#' Useful when underlying data is generated from multiple (2+) factors
#' @param r2.df R2 dataframe, output from CalculateR2()
#' @return Boxplot of R2 values by imputation type
#' @export
PlotR2MultFact <- function(r2.df){

  p <- ggplot(r2.df, aes(y = r2.df$R2, x = r2.df$Type)) +
    geom_boxplot() +
    scale_x_discrete(labels=c("True Correlation","One Factor","Two Factors","Three Factors","Four Factors","Mice","Ignoring Heterogeneity")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
    xlab("Imputation Type")

  return(p)
}

#' Plots R2 of one factor imputation
#'
#' Plots R2 in line chart of different imputation methods.
#' Useful when underlying data is generated from a single factor.
#' @param r2.df R2 dataframe, output from CalculateR2()
#' @return Line chart of R2 values by imputation type
#' @export
PlotR2OneFact <- function(r2.df){
  p <- ggplot(r2.df, aes(x=r2.df$Factor1, y=r2.df$R2, group=r2.df$Type)) +
    geom_line(aes(color = r2.df$Type))+
    geom_point(aes(color = r2.df$Type)) +
    scale_color_discrete(name = "", labels = c("True Imputed  \n vs. True Value",
                                               "Factor1 Imputed \n vs True Value",
                                               "Factor2 Imputed \n vs True Value",
                                               "Factor3 Imputed \n vs True Value",
                                               "Factor4 Imputed \n vs True Value",
                                               "Mice Imputed \n vs True Value",
                                               "Observed Imputed \n vs. True Value")) +
    theme_bw() +
    theme(legend.key.size = unit(1.5, "cm"),
          legend.key.width = unit(0.5,"cm"))
  return(p)
}

#' Creates true data, obs data, true correlation, and obs correlation
#'
#' First step of imputation pipeline. Correlation structure is m_0 + (m_1xF_1) + (m_2xF_2) + (m_3xF_3)
#' factor_mat is a matrix which contains the values for the 3 different factors. Each column corresponds to a factor and
#' each row is a value. By setting the third column (factor) to zero, we can generate data from a 2 factor model instead of 3.
#' Suppose that factor_mat has 30 rows (30 different combinations of factors), then 120 (30 x 4) data sets are generated in total.
#' num_to_remove metabolites are removed in four increments. Set num_to_remove=5 and that the first row of factor mat is a, b, c.
#' Under the same factors, in the first data set metabolites 1:5 are removed, in the second data set 6:10 are removed, in the third
#' 11:15, and in the fourth 16:20. This then repeats with the second row of factor mat.
#'
#' @param m_0 Residual correlation matrix
#' @param m_1 Correlation matrix for factor 1
#' @param m_2 Correlation matrix for factor 2
#' @param m_3 Correlation matrix for factor 3
#' @param num_of_metabos Number of metabolites
#' @param sample_size Number of individuals generated from correlation matrices
#' @param num_to_remove Number of metabolites to remove. Must be smaller than 1/4*num_of_metabolites
#' @param factor_mat Matrix where 3 columns are three factors and rows are values of each factor
#' @return true data, obs data, true correlation, and obs correlation
#' @import rlist
#' @export
CreateCorrelation <- function(m_0,m_1,m_2,m_3,
                              num_of_metabos=20,
                              sample_size=1000,
                              num_to_remove=5,
                              factor_mat){

  num_of_study <- dim(factor_mat)[1] * 4
  z_i <- matrix(rep(t(factor_mat),4),ncol=ncol(factor_mat),byrow=TRUE)

  corr_array_true <- array(NA, dim = c(num_of_metabos,num_of_metabos,num_of_study))
  for (i in 1:num_of_study){
    corr_array_true[,,i] <- as.matrix(nearPD(m_0 + (m_1*z_i[i,1]) + (m_2*z_i[i,2]) + (m_3*z_i[i,3]),corr=TRUE)$mat)
  }

  metabolites_array_true <- array(NA, dim = c(sample_size,num_of_metabos,num_of_study))
  for (i in 1:num_of_study){
    metabolites_array_true[,,i] <- mvrnorm(n = sample_size, rep(0,num_of_metabos), corr_array_true[,,i])
  }

  MissingObs <- CreateMissingAndObsList(dim(factor_mat)[1],num_to_remove, num_of_metabos)
  missing_list <- MissingObs[[1]]
  obs_list <- MissingObs[[2]]

  metabolites_array_obs <- metabolites_array_true

  for (i in 1:num_of_study){
    metabolites_array_obs[,missing_list[[i]],i] <- NA
  }

  corr_array_obs <- array(NA, dim = c(num_of_metabos,num_of_metabos,num_of_study))

  for (i in 1:num_of_study){
    corr_array_obs[,,i] <- cor(metabolites_array_obs[,,i],method = "spearman")
  }

  data_to_return <- list(metabolites_array_true,metabolites_array_obs)
  corr_to_return <- list(corr_array_true,corr_array_obs)

  return(list(data_to_return,corr_to_return))
}

#' Preps data and runs Mice
#'
#' From CreateCorrelation() correlation matrices are in 3d arrays with dimensions num_of_metabo x num_of_metabo x total number of studies
#' First we change the 3d array into a 2d matrix with dimensions (num_of_metabos x (num_of_metabos-1)/2 by total number of studies
#' Each column corresponds to an individual study and each row to a pairwise correlation of metabolites. Then we use the fischer
#' transform and run MICE on the transpose of the prepped 2d matrix. 5 MICE iputations are done and then the median is taken for each pairwise value
#' corr_array_mice is a 3d array with dimensions num_of_metabo x num_of_metabo x total number of studies that has been un-fischer transformed.
#' reorganized_corr_imp_med is a 2d matrix with dimensions (num_of_metabos x (num_of_metabos-1)/2 by total number of studies that is still fischer transformed.
#' reorganized_corr_imp_med will be the input for the next step of the imputation pipeline.
#'
#' @param corr_array_obs Observed correlation array. Output from CreateCorrelation()
#' @param factor_mat Matrix where 3 columns are three factors and rows are values of each factor
#' @param num_of_metabos Number of metabolites
#' @param num_of_predictors Number of other metabolites for MICE to use, set to 0 to use quickpred()
#'
#' @return un-fischer transformed 3d array corr_array_mice and fischer transformed 2d matrix reorganized_corr_imp_med
#' @import mice
#' @import Matrix
#' @import reshape2
#' @export
ReorganizeAndMICE <- function(corr_array_obs,
                              factor_mat,
                              num_of_metabos = 20,
                              num_of_predictors = 0){

  num_of_study <- dim(factor_mat)[1] * 4
  z_i <- matrix(rep(t(factor_mat),4),ncol=ncol(factor_mat),byrow=TRUE)

  corr_obs.df <- melt(corr_array_obs[,,])

  corr_obs.df <- corr_obs.df[corr_obs.df$Var1 < corr_obs.df$Var2,]
  corr_obs.df$value <- log((1 + corr_obs.df$value)/(1 - corr_obs.df$value)) / 2

  num_of_imputs <- 5
  num_of_pairs <- num_of_metabos * (num_of_metabos-1)/2
  reorganized_corr_obs <- matrix(-1,num_of_pairs,num_of_study)
  reorganized_corr_imp <- matrix(0,num_of_pairs,num_of_study)
  reorganized_corr_imp_med <- matrix(0,num_of_pairs,num_of_study)
  reorganized_corr_array_imp <- array(0,dim = c(num_of_pairs,num_of_study,num_of_imputs))

  for (i in 1:num_of_study){
    reorganized_corr_obs[,i] <- corr_obs.df$value[corr_obs.df$Var3 == i]
  }

  if (!num_of_predictors){
    pred <- quickpred(t(reorganized_corr_obs))
  } else {
    grid <- seq(0, 1, 0.05)
    result <- apply(sapply(grid, function(x) { rowSums(quickpred(t(reorganized_corr_obs), mincor=x)) } ), 1, function(x) {min(which(x<=num_of_predictors))})
    pred <- quickpred(t(reorganized_corr_obs), mincor=grid[result])
  }

  reorganized_corr_obs_mice  <- mice(t(reorganized_corr_obs),
                                     m=num_of_imputs,
                                     threshold=.999999,
                                     predictorMatrix = pred,
                                     method = c("pmm"),
                                     remove.collinear = FALSE)
  for (i in 1:num_of_imputs){
    reorganized_corr_array_imp[,,i] <- t(as.matrix(complete(reorganized_corr_obs_mice,i)))
  }

  for (i in 1:dim(reorganized_corr_imp)[1]){
    for (j in 1:dim(reorganized_corr_imp)[2]){
      reorganized_corr_imp_med[i,j] <- median(reorganized_corr_array_imp[i,j,])
    }
  }

  ###########Take median of all mice outputs and impute
  corr_array_mice <- array(0, dim = c(num_of_metabos,num_of_metabos,num_of_study))
  reorganized_corr_imp_med <- as.matrix(reorganized_corr_imp_med)

  reorganized_corr_imp_med_detran <- (exp(2 * reorganized_corr_imp_med) - 1) / (exp(2 * reorganized_corr_imp_med) + 1)
  # reorganized_corr_imp_med_detran <- reorganized_corr_imp_med

  for (i in 1:num_of_study){
    corr_array_mice[,,i] <- as.matrix(nearPD(Vec2Cor(reorganized_corr_imp_med_detran[,i],num_of_metabos),corr=TRUE)$mat)
    diag(corr_array_mice[,,i]) <- diag(corr_array_mice[,,i]) + .01
    corr_array_mice[abs(corr_array_mice)<.01] <- 0
  }



  return(list(corr_array_mice,reorganized_corr_imp_med))
}

#' Runs factor analysis
#'
#' Runs factor analysis with 1-4 factors. Converts input from 2d to 3d form (such as in the description of ReogranizeAndMICE()).
#' Outputs both factor analysis intermediates (loadings, scores, ...) and estimated correlation for each study.
#' Input is second output from ReorganizeAndMICE(), where the columns are studies and the rows are pairwise correlation.
#' Runs factor analysis on transpose of input.
#' @param reorganized_corr_mice_med Second output from ReorganizeAndMICE(). Fischer transformed 2d matrix
#' @param num_of_metabos Number of metabolites
#' @return list of lists. First list is factor analysis intermediates for 1-4 factors, second is 3d corr matrices for 1-4 factors
#' @import psych
#' @export
FactorAnalysis <- function(reorganized_corr_mice_med,num_of_metabos = 20){

  num_of_study <- dim(reorganized_corr_mice_med)[2]
  reorganized_corr_fa1 <- array(0, dim = c(dim(reorganized_corr_mice_med)[1],dim(reorganized_corr_mice_med)[2]))
  reorganized_corr_fa2 <- array(0, dim = c(dim(reorganized_corr_mice_med)[1],dim(reorganized_corr_mice_med)[2]))
  reorganized_corr_fa3 <- array(0, dim = c(dim(reorganized_corr_mice_med)[1],dim(reorganized_corr_mice_med)[2]))
  reorganized_corr_fa4 <- array(0, dim = c(dim(reorganized_corr_mice_med)[1],dim(reorganized_corr_mice_med)[2]))

  faRes1  <- fa(t(reorganized_corr_mice_med),
                nfactors = 1)
  faRes2  <- fa(t(reorganized_corr_mice_med),
                nfactors = 2)
  faRes3  <- fa(t(reorganized_corr_mice_med),
                nfactors = 3)
  faRes4  <- fa(t(reorganized_corr_mice_med),
                nfactors = 4)

  mult <- sqrt(apply(reorganized_corr_mice_med,1,var))
  for (i in 1:num_of_study){
    reorganized_corr_fa1[,i] <- rowMeans(reorganized_corr_mice_med) + mult*faRes1$loadings %*% faRes1$scores[i,]
    reorganized_corr_fa2[,i] <- rowMeans(reorganized_corr_mice_med) + mult*faRes2$loadings %*% faRes2$scores[i,]
    reorganized_corr_fa3[,i] <- rowMeans(reorganized_corr_mice_med) + mult*faRes3$loadings %*% faRes3$scores[i,]
    reorganized_corr_fa4[,i] <- rowMeans(reorganized_corr_mice_med) + mult*faRes4$loadings %*% faRes4$scores[i,]
  }

  reorganized_corr_fa1 <- (exp(2 * reorganized_corr_fa1) - 1) / (exp(2 * reorganized_corr_fa1) + 1)
  reorganized_corr_fa2 <- (exp(2 * reorganized_corr_fa2) - 1) / (exp(2 * reorganized_corr_fa2) + 1)
  reorganized_corr_fa3 <- (exp(2 * reorganized_corr_fa3) - 1) / (exp(2 * reorganized_corr_fa3) + 1)
  reorganized_corr_fa4 <- (exp(2 * reorganized_corr_fa4) - 1) / (exp(2 * reorganized_corr_fa4) + 1)

  corr_array_fa1 <- array(0, dim = c(num_of_metabos,num_of_metabos,num_of_study))
  corr_array_fa2 <- array(0, dim = c(num_of_metabos,num_of_metabos,num_of_study))
  corr_array_fa3 <- array(0, dim = c(num_of_metabos,num_of_metabos,num_of_study))
  corr_array_fa4 <- array(0, dim = c(num_of_metabos,num_of_metabos,num_of_study))
  for (i in 1:num_of_study){
    corr_array_fa1[,,i] <- as.matrix(nearPD(Vec2Cor(reorganized_corr_fa1[,i],num_of_metabos),corr=TRUE)$mat)
    corr_array_fa2[,,i] <- as.matrix(nearPD(Vec2Cor(reorganized_corr_fa2[,i],num_of_metabos),corr=TRUE)$mat)
    corr_array_fa3[,,i] <- as.matrix(nearPD(Vec2Cor(reorganized_corr_fa3[,i],num_of_metabos),corr=TRUE)$mat)
    corr_array_fa4[,,i] <- as.matrix(nearPD(Vec2Cor(reorganized_corr_fa4[,i],num_of_metabos),corr=TRUE)$mat)
  }

  for (i in 1:num_of_study){
    diag(corr_array_fa1[,,i]) <- diag(corr_array_fa1[,,i]) + .01
    diag(corr_array_fa2[,,i]) <- diag(corr_array_fa2[,,i]) + .01
    diag(corr_array_fa3[,,i]) <- diag(corr_array_fa3[,,i]) + .01
    diag(corr_array_fa4[,,i]) <- diag(corr_array_fa4[,,i]) + .01

    corr_array_fa1[abs(corr_array_fa1)<.01] <- 0
    corr_array_fa2[abs(corr_array_fa2)<.01] <- 0
    corr_array_fa3[abs(corr_array_fa3)<.01] <- 0
    corr_array_fa4[abs(corr_array_fa4)<.01] <- 0

  }
  return(list(list(faRes1,faRes2,faRes3,faRes4),list(corr_array_fa1,corr_array_fa2,corr_array_fa3,corr_array_fa4)))
}

#' Calculates Correlation Ignoring Heterogeneity
#'
#' @param corr_array_obs Output from CreateCorrelation()
#' @param num_of_metabos Number of metabolites
#' @return Correlation matrix with dimension num_of_metabo by num_of_metabo, ignores heterogeneity
#' @export
CalcCorrAvg <- function(corr_array_obs,num_of_metabos){
  corr_avg <- matrix(NA,num_of_metabos,num_of_metabos)

  for (i in 1:num_of_metabos){
    for (j in 1:num_of_metabos){
      corr_avg[i,j] <- mean(corr_array_obs[i,j,],na.rm = T)
    }
  }
  return(corr_avg)
}

#' Imputes missing metabolites based on various correlation matrices
#'
#' Imputes missing metabolites based on different correlation matrices calculated previously.
#' The different methods are: average correlation (ignoring heterogeneity), true correlation,
#' MICE imputation and 1-4 factor analysis, and MICE imputation. Outputs into dataframe. Takes time
#' (maybe a lot if large sample size) but useful for ggplot.
#'
#' @param num_of_metabos Number of metabolites
#' @param sample_size Number of individuals generated from correlation matrices
#' @param num_to_remove Number of metabolites removed at a time from true data in CreateCorrelation()
#' @param factor_mat Matrix where 3 columns are three factors and rows are values of each factor
#' @param metabolites_array_obs Metabolite data with removed metabolites generated in CreateCorrelation()
#' @param metabolites_array_true Metabolite data generated in CreateCorrelation()
#' @param corr_avg Avgerage correlation. Ignores heterogeneity. Output from CalcCorrAvg()
#' @param corr_array_true True correlation matrix. An output from CreateCorrelation()
#' @param corr_array_fa1 MICE and 1 factor factor analysis. An output from FactorAnalysis()
#' @param corr_array_fa2 MICE and 2 factor factor analysis. An output from FactorAnalysis()
#' @param corr_array_fa3 MICE and 3 factor factor analysis. An output from FactorAnalysis()
#' @param corr_array_fa4 MICE and 4 factor factor analysis. An output from FactorAnalysis()
#' @param corr_array_mice Imputed with MICE. First output from ReorganizeAndMICE()
#' @return data frame of correlation values by metabolite, factor, and imputation type
#' @export
ImputeMissingMetabolites <- function(sample_size = 1000,num_of_metabos = 20,num_to_remove = 5,
                                     factor_mat,
                                     metabolites_array_obs,
                                     metabolites_array_true,
                                     corr_avg,
                                     corr_array_true,
                                     corr_array_fa1,
                                     corr_array_fa2,
                                     corr_array_fa3,
                                     corr_array_fa4,
                                     corr_array_mice){

  num_of_study <- dim(factor_mat)[1] * 4
  z_i <- matrix(rep(t(factor_mat),4),ncol=ncol(factor_mat),byrow=TRUE)

  MissingObs <- CreateMissingAndObsList(dim(factor_mat)[1],num_to_remove, num_of_metabos)
  missing_list <- MissingObs[[1]]
  obs_list <- MissingObs[[2]]

  metabolites_array_obs_imput <- array(NA, dim = c(sample_size,num_of_metabos,num_of_study))
  metabolites_array_true_imput <- array(NA, dim = c(sample_size,num_of_metabos,num_of_study))
  metabolites_array_fact_imput1 <- array(NA, dim = c(sample_size,num_of_metabos,num_of_study))
  metabolites_array_fact_imput2 <- array(NA, dim = c(sample_size,num_of_metabos,num_of_study))
  metabolites_array_fact_imput3 <- array(NA, dim = c(sample_size,num_of_metabos,num_of_study))
  metabolites_array_fact_imput4 <- array(NA, dim = c(sample_size,num_of_metabos,num_of_study))
  metabolites_array_mice_imput <- array(NA, dim = c(sample_size,num_of_metabos,num_of_study))
  for (i in 1:num_of_study){
    metabolites_array_obs_imput[,,i]  <- ImputeMissingMu(metabolites_array_obs[,,i],corr_avg,missing_list[[i]],obs_list[[i]],sample_size)
    metabolites_array_true_imput[,,i] <- ImputeMissingMu(metabolites_array_obs[,,i],corr_array_true[,,i],missing_list[[i]],obs_list[[i]],sample_size)
    metabolites_array_fact_imput1[,,i] <- ImputeMissingMu(metabolites_array_obs[,,i],corr_array_fa1[,,i],missing_list[[i]],obs_list[[i]],sample_size)
    metabolites_array_fact_imput2[,,i] <- ImputeMissingMu(metabolites_array_obs[,,i],corr_array_fa2[,,i],missing_list[[i]],obs_list[[i]],sample_size)
    metabolites_array_fact_imput3[,,i] <- ImputeMissingMu(metabolites_array_obs[,,i],corr_array_fa3[,,i],missing_list[[i]],obs_list[[i]],sample_size)
    metabolites_array_fact_imput4[,,i] <- ImputeMissingMu(metabolites_array_obs[,,i],corr_array_fa4[,,i],missing_list[[i]],obs_list[[i]],sample_size)
    metabolites_array_mice_imput[,,i] <- ImputeMissingMu(metabolites_array_obs[,,i],corr_array_mice[,,i],missing_list[[i]],obs_list[[i]],sample_size)
  }

  imputed.df <- data.frame(ObsImp=double(),
                           TrueImp=double(),
                           FactImp1=double(),
                           FactImp2=double(),
                           FactImp3=double(),
                           FactImp4=double(),
                           MiceImp=double(),
                           TrueVal=double(),
                           MetNum=integer(),
                           Zind1=double(),
                           Zind2=double(),
                           Zind3=double(),
                           stringsAsFactors=FALSE)

  for(study_ind in 1:num_of_study){

    if (study_ind == floor(num_of_study * .1)){
      print("10% Done")
    }
    if (study_ind == floor(num_of_study * .5)){
      print("50% Done")
    }
    if (study_ind == floor(num_of_study * .9)){
      print("90% Done")
    }


    for (metabo_num in missing_list[[study_ind]]){
      obs_imp <- metabolites_array_obs_imput[,metabo_num,study_ind]
      true_imp <- metabolites_array_true_imput[,metabo_num,study_ind]
      fact_imp1 <- metabolites_array_fact_imput1[,metabo_num,study_ind]
      fact_imp2 <- metabolites_array_fact_imput2[,metabo_num,study_ind]
      fact_imp3 <- metabolites_array_fact_imput3[,metabo_num,study_ind]
      fact_imp4 <- metabolites_array_fact_imput4[,metabo_num,study_ind]
      mice_imp <- metabolites_array_mice_imput[,metabo_num,study_ind]
      true_val <- metabolites_array_true[,metabo_num,study_ind]
      to_add <- cbind(as.double(obs_imp),as.double(true_imp),as.double(fact_imp1),as.double(fact_imp2),as.double(fact_imp3),as.double(fact_imp4),as.double(mice_imp),as.double(true_val),as.integer(metabo_num),as.double(z_i[study_ind,1]),as.double(z_i[study_ind,2]),as.double(z_i[study_ind,3]))
      imputed.df <-rbind(imputed.df,to_add)
    }
  }

  colnames(imputed.df) <- c("ObsImp","TrueImp","FactImp1","FactImp2","FactImp3","FactImp4","MiceImp","TrueVal","Metabolite","Factor1","Factor2","Factor3")

  return(imputed.df)
}

#' Calculates R2 across factors for different imputation types
#'
#' Creates data frame of R2 values across factors and imputation type. May take time with large sample size.
#' @param imputed.df dataframe of imputed values. Output from ImputeMissingMetabolites()
#' @param factor_mat Matrix where 3 columns are three factors and rows are values of each factor
#' @return data frame of R2 values by imputation type
#' @export
CalculateR2 <- function(imputed.df,factor_mat){
  r2.df <- data.frame(R2=double(),
                      Zind1=double(),
                      Zind2=double(),
                      Zind3=double(),
                      group=character(),
                      stringsAsFactors=FALSE)

  num_of_study <- dim(factor_mat)[1] * 4
  z_i <- matrix(rep(t(factor_mat),4),ncol=ncol(factor_mat),byrow=TRUE)

  for(factor_ind in 1:dim(factor_mat)[1]){

    if (factor_ind == floor(dim(factor_mat)[1] * .1)){
      print("10% Done")
    }
    if (factor_ind == floor(dim(factor_mat)[1] * .5)){
      print("50% Done")
    }
    if (factor_ind == floor(dim(factor_mat)[1] * .9)){
      print("90% Done")
    }

    true_r <- cor(imputed.df[imputed.df$Factor1 == z_i[factor_ind,1] & imputed.df$Factor2 == z_i[factor_ind,2] & imputed.df$Factor3 == z_i[factor_ind,3],]$TrueImp,
                  imputed.df[imputed.df$Factor1 == z_i[factor_ind,1] & imputed.df$Factor2 == z_i[factor_ind,2] & imputed.df$Factor3 == z_i[factor_ind,3],]$TrueVal)^2
    to_add <- cbind(as.double(true_r),as.double(z_i[factor_ind,1]),as.double(z_i[factor_ind,2]),as.double(z_i[factor_ind,3]),"TrueIvsTrueV")
    r2.df <-rbind(r2.df,to_add)

    obs_r <- cor(imputed.df[imputed.df$Factor1 == z_i[factor_ind,1] & imputed.df$Factor2 == z_i[factor_ind,2] & imputed.df$Factor3 == z_i[factor_ind,3],]$ObsImp,
                 imputed.df[imputed.df$Factor1 == z_i[factor_ind,1] & imputed.df$Factor2 == z_i[factor_ind,2] & imputed.df$Factor3 == z_i[factor_ind,3],]$TrueVal)^2
    to_add <- cbind(as.double(obs_r),as.double(z_i[factor_ind,1]),as.double(z_i[factor_ind,2]),as.double(z_i[factor_ind,3]),"ObsIvsTrueV")
    r2.df <-rbind(r2.df,to_add)

    fact_r <- cor(imputed.df[imputed.df$Factor1 == z_i[factor_ind,1] & imputed.df$Factor2 == z_i[factor_ind,2] & imputed.df$Factor3 == z_i[factor_ind,3],]$FactImp1,
                  imputed.df[imputed.df$Factor1 == z_i[factor_ind,1] & imputed.df$Factor2 == z_i[factor_ind,2] & imputed.df$Factor3 == z_i[factor_ind,3],]$TrueVal)^2
    to_add <- cbind(as.double(fact_r),as.double(z_i[factor_ind,1]),as.double(z_i[factor_ind,2]),as.double(z_i[factor_ind,3]),"FactI1vsTrueV")
    r2.df <-rbind(r2.df,to_add)

    fact_r <- cor(imputed.df[imputed.df$Factor1 == z_i[factor_ind,1] & imputed.df$Factor2 == z_i[factor_ind,2] & imputed.df$Factor3 == z_i[factor_ind,3],]$FactImp2,
                  imputed.df[imputed.df$Factor1 == z_i[factor_ind,1] & imputed.df$Factor2 == z_i[factor_ind,2] & imputed.df$Factor3 == z_i[factor_ind,3],]$TrueVal)^2
    to_add <- cbind(as.double(fact_r),as.double(z_i[factor_ind,1]),as.double(z_i[factor_ind,2]),as.double(z_i[factor_ind,3]),"FactI2vsTrueV")
    r2.df <-rbind(r2.df,to_add)

    fact_r <- cor(imputed.df[imputed.df$Factor1 == z_i[factor_ind,1] & imputed.df$Factor2 == z_i[factor_ind,2] & imputed.df$Factor3 == z_i[factor_ind,3],]$FactImp3,
                  imputed.df[imputed.df$Factor1 == z_i[factor_ind,1] & imputed.df$Factor2 == z_i[factor_ind,2] & imputed.df$Factor3 == z_i[factor_ind,3],]$TrueVal)^2
    to_add <- cbind(as.double(fact_r),as.double(z_i[factor_ind,1]),as.double(z_i[factor_ind,2]),as.double(z_i[factor_ind,3]),"FactI3vsTrueV")
    r2.df <-rbind(r2.df,to_add)

    fact_r <- cor(imputed.df[imputed.df$Factor1 == z_i[factor_ind,1] & imputed.df$Factor2 == z_i[factor_ind,2] & imputed.df$Factor3 == z_i[factor_ind,3],]$FactImp4,
                  imputed.df[imputed.df$Factor1 == z_i[factor_ind,1] & imputed.df$Factor2 == z_i[factor_ind,2] & imputed.df$Factor3 == z_i[factor_ind,3],]$TrueVal)^2
    to_add <- cbind(as.double(fact_r),as.double(z_i[factor_ind,1]),as.double(z_i[factor_ind,2]),as.double(z_i[factor_ind,3]),"FactI4vsTrueV")
    r2.df <-rbind(r2.df,to_add)

    mice_r <- cor(imputed.df[imputed.df$Factor1 == z_i[factor_ind,1] & imputed.df$Factor2 == z_i[factor_ind,2] & imputed.df$Factor3 == z_i[factor_ind,3],]$MiceImp,
                  imputed.df[imputed.df$Factor1 == z_i[factor_ind,1] & imputed.df$Factor2 == z_i[factor_ind,2] & imputed.df$Factor3 == z_i[factor_ind,3],]$TrueVal)^2
    to_add <- cbind(as.double(mice_r),as.double(z_i[factor_ind,1]),as.double(z_i[factor_ind,2]),as.double(z_i[factor_ind,3]),"MiceIvsTrueV")
    r2.df <-rbind(r2.df,to_add)

  }

  colnames(r2.df) <- c("R2","Factor1","Factor2","Factor3","Type")
  r2.df$R2 <- as.double(r2.df$R2)
  r2.df$Factor1 <- as.double(r2.df$Factor1)
  r2.df$Factor2 <- as.double(r2.df$Factor2)
  r2.df$Factor3 <- as.double(r2.df$Factor3)
  r2.df$Type <- factor(r2.df$Type, levels = c("TrueIvsTrueV","FactI1vsTrueV","FactI2vsTrueV", "FactI3vsTrueV","FactI4vsTrueV","MiceIvsTrueV","ObsIvsTrueV"))

  return(r2.df)
}

#' Generates Swarm File
#' @param seed_max Max seed number. Runs this many simulations
#' @param swarm_name Name of swarm file output
#' @param r_name Name of r file to run
#' @param low_val_1 low value for m_1
#' @param high_val_1 high value for m_1
#' @param low_val_2 low value for m_2
#' @param high_val_2 high value for m_2
#' @return writes out swarm file
#' @export
SwarmGenerator <- function(seed_max, swarm_name,r_name, low_val_1,high_val_1,low_val_2,high_val_2){
  output_swarm <- ""
  for (seed in 1:seed_max){
    output_swarm <- paste(output_swarm,"module load R; Rscript ",r_name," ",
                          seed," ",
                          low_val_1," ",
                          high_val_1," ",
                          low_val_2," ",
                          high_val_2,";\n", sep = "")
  }

  writeLines(output_swarm, swarm_name)
}
