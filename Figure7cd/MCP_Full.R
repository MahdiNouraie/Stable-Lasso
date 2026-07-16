# This script generates Figures 4(c) and 4(d) from the paper
options(warn=-1) # Turn off warnings
if (!requireNamespace("MASS")) {install.packages("MASS")} # install package if not already installed
if (!requireNamespace("glmnet")) {install.packages("glmnet")} # install package if not already installed
if (!requireNamespace("ncvreg")) {install.packages("ncvreg")} # install package if not already installed
if (!requireNamespace("ggplot2")) {install.packages("ggplot2")} # install package if not already installed
if (!requireNamespace("reshape2")) {install.packages("reshape2")} # install package if not already installed
if (!requireNamespace("devtools")) {install.packages("devtools")} # install package if not already installed
if (!requireNamespace("parallel")) {install.packages("parallel")} # install package if not already installed
devtools::install_github("JLSteenwyk/ggpubfigs", force = TRUE)

# Load necessary library
library(MASS) # load the MASS package for mvrnorm
library(glmnet) # load the glmnet package for Lasso
library(ncvreg) # Load the ncvreg package for SCAD
library(ggplot2) # Load the ggplot2 package
library(reshape2) # Load the reshape2 package
library(parallel) # Load the parallel package for parallel processing
library(ggpubfigs) # Load the ggpubfigs package for color-blind friendly colors

#sessionInfo()
#R version 4.4.2 (2024-10-31)
#Platform: aarch64-apple-darwin20
#Running under: macOS Sequoia 15.4.1

#Matrix products: default
#BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
#LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

#locale:
#  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#time zone: Australia/Sydney
#tzcode source: internal

#attached base packages:
#  [1] parallel  stats     graphics  grDevices utils     datasets  methods  
#[8] base     

#other attached packages:
#  [1] ggpubfigs_0.0.1 reshape2_1.4.4  ggplot2_3.5.2   ncvreg_3.15.0  
#[5] glmnet_4.1-8    Matrix_1.7-1    MASS_7.3-64    

#loaded via a namespace (and not attached):
#  [1] generics_0.1.3     shape_1.4.6.1      stringi_1.8.7     
#[4] lattice_0.22-6     digest_0.6.37      magrittr_2.0.3    
#[7] grid_4.4.2         RColorBrewer_1.1-3 iterators_1.0.14  
#[10] pkgload_1.4.0      fastmap_1.2.0      foreach_1.5.2     
#[13] plyr_1.8.9         processx_3.8.4     sessioninfo_1.2.3 
#[16] pkgbuild_1.4.6     survival_3.8-3     ps_1.8.1          
#[19] urlchecker_1.0.1   promises_1.3.2     purrr_1.0.2       
#[22] scales_1.4.0       codetools_0.2-20   cli_3.6.5         
#[25] shiny_1.10.0       rlang_1.1.6        ellipsis_0.3.2    
#[28] splines_4.4.2      withr_3.0.2        remotes_2.5.0     
#[31] cachem_1.1.0       devtools_2.4.5     tools_4.4.2       
#[34] memoise_2.0.1      dplyr_1.1.4        httpuv_1.6.15     
#[37] curl_6.0.1         vctrs_0.6.5        R6_2.6.1          
#[40] mime_0.12          lifecycle_1.0.4    stringr_1.5.1     
#[43] fs_1.6.5           htmlwidgets_1.6.4  usethis_3.1.0     
#[46] miniUI_0.1.1.1     callr_3.7.6        desc_1.4.3        
#[49] pkgconfig_2.0.3    pillar_1.11.0      later_1.4.1       
#[52] gtable_0.3.6       glue_1.8.0         profvis_0.4.0     
#[55] Rcpp_1.1.0         tibble_3.3.0       tidyselect_1.2.1  
#[58] rstudioapi_0.17.1  farver_2.1.2       xtable_1.8-4      
#[61] htmltools_0.5.8.1  compiler_4.4.2   


# Stability measure (2018) from "https://github.com/nogueirs/JMLR2018/blob/master/R/getStability.R"
getStability <- function(X,alpha=0.05) {
  ## the input X is a binary matrix of size M*d where:
  ## M is the number of bootstrap replicates
  ## d is the total number of features
  ## alpha is the level of significance (e.g. if alpha=0.05, we will get 95% confidence intervals)
  ## it's an optional argument and is set to 5% by default
  ### first we compute the stability
  
  M<-nrow(X)
  d<-ncol(X)
  hatPF<-colMeans(X)
  kbar<-sum(hatPF)
  v_rand=(kbar/d)*(1-kbar/d)
  stability<-1-(M/(M-1))*mean(hatPF*(1-hatPF))/v_rand ## this is the stability estimate
  
  ## then we compute the variance of the estimate
  ki<-rowSums(X)
  phi_i<-rep(0,M)
  for(i in 1:M){ 
    phi_i[i]<-(1/v_rand)*((1/d)*sum(X[i,]*hatPF)-(ki[i]*kbar)/d^2-(stability/2)*((2*kbar*ki[i])/d^2-ki[i]/d-kbar/d+1))
  }
  phi_bar=mean(phi_i)
  var_stab=(4/M^2)*sum((phi_i-phi_bar)^2) ## this is the variance of the stability estimate
  
  ## then we calculate lower and upper limits of the confidence intervals
  z<-qnorm(1-alpha/2) # this is the standard normal cumulative inverse at a level 1-alpha/2
  upper<-stability+z*sqrt(var_stab) ## the upper bound of the (1-alpha) confidence interval
  lower<-stability-z*sqrt(var_stab) ## the lower bound of the (1-alpha) confidence interval
  
  return(list("stability"=stability,"variance"=var_stab,"lower"=lower,"upper"=upper))
  
}
# Newton's method for optimization
newton <- function(f, fp, x, tol = 1e-3, m = 100) {
  iter <- 0
  oldx <- x
  x <- oldx + 10 * tol
  
  while (abs(x - oldx) > tol) {
    iter <- iter + 1
    if (iter > m) {
      warning("Maximum iterations (m = ", m,
              ") reached: returning last estimate", call. = FALSE)
      break
    }
    oldx <- x
    x <- x - f(x) / fp(x)
  }
  
  return(x)
}
# AirHOLP (2025) from "https://github.com/Logic314/Air-HOLP"
AirHOLP <- function(X, y, Threshold = min(dim(X)[2]-1, dim(X)[1]/log(dim(X)[1])),
                    r0 = 10, adapt = TRUE, iter = 10, Lambda, U, XU) {
  # Arguments:-
  # X: matrix of features (Matrix)
  # y: response vector (Vector)
  # Threshold: screening threshold (Integer)
  # r0: initial penalties (Vector)
  # adapt: if >= 1 adaptive penalty will be used (Binary)
  # iter: maximum number of iteration for adaptive penalty selection (Integer)
  # Lambda: eigenvalues of XXT, if missing the function will compute it (Vector)
  # U: eigenvectors of XXT, if missing the function will compute it (Matrix)
  # XU: X transpose times U, if missing the function will compute it (Matrix)
  
  # Output:-
  # index_r: ranking of features by Air-HOLP (Matrix)
  # index_r0: ranking of features by Ridge-HOLP (Matrix)
  # Beta_r: regression coefficients of Air-HOLP (Matrix)
  # Beta_r0: regression coefficients of Ridge-HOLP (Matrix)
  # r: selected penalty parameters by Air-HOLP (Vector)
  # iter_last: number of iterations used in Air-HOLP (Vector)
  
  n <- dim(X)[1] # sample size
  p <- dim(X)[2] # number of features
  q <- length(r0) # number of penalty parameters
  iter_temp2 <- 0*(1:q) # used for calculating iter_last
  iter_temp1 <- iter_temp2 - 1 # used for calculating iter_last
  
  # Standardizing X and y:
  X <- X - matrix(rep(colMeans(X),each = n),n,p)
  X <- X/matrix(rep(sqrt(colMeans(X^2)),each = n),n,p)
  y <- (y - mean(y))/sd(y)
  
  if(adapt){
    # Main computations:
    if(missing(Lambda)|missing(U)){
      XXT <- tcrossprod(X)
      eXXT <- eigen(XXT)
      Lambda <- eXXT$values
      U <- eXXT$vectors
    }
    if(missing(XU)){
      XU <- crossprod(X,U)
    }
    Dn <- diag(Lambda)
    UTy <- crossprod(U,y)
    yUD2UTy <- UTy^2*(Lambda^2)
    
    # Penalty selection:
    r_max <- 1000*sqrt(n) # maximum penalty
    max.iter <- 30 # maximum number of iterations for Newtons method
    index_r <- matrix(1:(p*q), nrow = p, ncol = q)
    index_r0 <- index_r
    Beta_r <- index_r
    Beta_r0 <- index_r
    r <- r0
    r_temp <- r0
    for (j in 1:iter) {
      for (i in 1:q) {
        # Initial screening:
        Beta_temp <- XU%*%((Lambda+r[i])^(-1)*UTy)
        index_temp <- match(1:p,rank(-abs(Beta_temp), na.last = NA,
                                     ties.method = c("random")))
        Xs <- X[,index_temp[1:Threshold]] # Screened features
        if(j<2) {
          Beta_r0[,i] <- Beta_temp
          index_r0[,i] <- rank(-abs(Beta_temp), na.last = NA,
                               ties.method = c("random"))
        }
        
        # Estimating the expected response:
        ys <- Xs%*%(solve(crossprod(Xs) +
                            diag(Threshold)*10^-12)%*%crossprod(Xs,y))
        
        # MSE functions:
        ysUDUTy <- t(crossprod(ys,U)*Lambda)*UTy
        Z <- function(lam) { # The function we minimize
          t((Lambda+lam)^-2)%*%yUD2UTy - 2*t((Lambda+lam)^-1)%*%ysUDUTy
        }
        Z1 <- function(lam) { # First derivative
          -2*t((Lambda+lam)^-3)%*%yUD2UTy + 2*t((Lambda+lam)^-2)%*%ysUDUTy
        }
        Z2 <- function(lam) { # Second derivative
          6*t((Lambda+lam)^-4)%*%yUD2UTy - 4*t((Lambda+lam)^-3)%*%ysUDUTy
        }
        
        # MSE minimization:
        sol <- newton(Z1, Z2, 0.0001, tol = 0.001, m = max.iter)
        r[i] <- sol
        if(r[i] > r_max) {r[i] <- r_max}
        if(r[i] < 0.0001) {r[i] <- 0.0001}
        if(Z(r_max) < Z(r[i])) {r[i] <- r_max} # Checking boundaries
        if(Z(0.0001) < Z(r[i])) {r[i] <- 0.0001}
        
        # Feature screening:
        Beta_r[,i] <- XU%*%((Lambda+r[i])^(-1)*UTy)
        index_r[,i] <- rank(-abs(Beta_r[,i]), na.last = NA,
                            ties.method = c("random"))
        
        # Calculations for the number of iterations:
        if(abs(r[i] - r_temp[i]) < 0.01*r[i]){ # Checking relative error
          iter_temp1[i] <- j
          iter_temp2[i] <- iter_temp2[i] + 1
        }
      }
      if(sum(abs(r - r_temp) < 0.01*r) == q){ # Checking relative error
        break
      }
      r_temp <- r
    }
    iter_last <- iter_temp1 - iter_temp2 + 1 # Number of iterations
    AirHOLP <- list(index_r = index_r, index_r0 = index_r0, Beta_r = Beta_r,
                    Beta_r0 = Beta_r0, r = r, iter_last = iter_last)
  } else{
    if(q < 2) {
      # Feature screening:
      if(missing(Lambda)|missing(U)){
        Beta_r0 <- crossprod(X, solve(tcrossprod(X)+r0*diag(n),y))
      } else{
        UTy <- crossprod(U,y)
        Beta_r0 <- XU%*%((Lambda+r0)^(-1)*UTy)
      }
      index_r0 <- rank(-abs(Beta_r0), na.last = NA, ties.method = c("random"))
      AirHOLP <- list(index_r0 = index_r0, Beta_r0 = Beta_r0)
    } else{
      # Main computations:
      if(missing(Lambda)|missing(U)){
        XXT <- tcrossprod(X)
        eXXT <- eigen(XXT)
        Lambda <- eXXT$values
        U <- eXXT$vectors
      }
      if(missing(XU)){
        XU <- crossprod(X,U)
      }
      Dn <- diag(Lambda)
      UTy <- crossprod(U,y)
      
      # Feature screening:
      index_r <- matrix(1:(p*q), nrow = p, ncol = q)
      for (i in 1:q) {
        Beta_r0[,i] <- XU%*%((Lambda+r0[i])^(-1)*UTy)
        index_r0[,i] <- rank(-abs(Beta_r0[,i]), na.last = NA,
                             ties.method = c("random"))
      }
      AirHOLP <- list(index_r0 = index_r0, Beta_r0 = Beta_r0)
    }
  }
}
# Function to generate data
generate_data <- function(n, p) {
  if (p %% 5 != 0) {
    stop("p must be divisible by 5")
  }
  
  # Number of variables per group
  group_size <- p / 5
  
  # Correlation values for each group
  rho_values <- c(0.8, 0.85, 0.9, 0.95, 0.99)
  
  # Initialize covariance matrix
  Sigma <- matrix(0, nrow = p, ncol = p)
  
  for (i in 1:5) {
    start_idx <- (i - 1) * group_size + 1
    end_idx <- i * group_size
    
    # Set correlation within the group
    Sigma[start_idx:end_idx, start_idx:end_idx] <- rho_values[i]
    diag(Sigma)[start_idx:end_idx] <- 1  # Set diagonal to 1
  }
  
  # Generate correlated predictors
  X <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  
  # Generate noise
  epsilon <- rnorm(n, mean = 0, sd = 1)
  
  # Define beta coefficients
  beta <- rep(0, p)  # Initialize all coefficients as zero
  active_values <- c(3, 2.5, 2, 1.5, 1)  # Active variable values
  active_indices <- seq(group_size, p, by = group_size)  # Last element of each group
  beta[active_indices] <- active_values  # Assign active variable values
  
  # Generate response variable
  response <- X %*% beta + epsilon
  
  # Combine predictors and response into a data frame
  data <- cbind(Y = response, X)
  
  return(data)
}
n <- 100  # Number of samples 
p <- 1000  # Number of predictors

# Stability Selection with MCP
time_taken <- system.time({results <- mclapply(1:100, function(i) {
  set.seed(26 + i) # Set seed for reproducibility
  
  Stability_values <- c() # Initialize an empty vector to store stability values
  TP6 <- c(); FP6 <- c(); FN6 <- c(); TN6 <- c() # Initialize vectors for TP, FP, FN, TN under thereshold 0.6
  TP7 <- c(); FP7 <- c(); FN7 <- c(); TN7 <- c() # Initialize vectors for TP, FP, FN, TN under thereshold 0.7
  TP8 <- c(); FP8 <- c(); FN8 <- c(); TN8 <- c() # Initialize vectors for TP, FP, FN, TN under thereshold 0.8
  TP9 <- c(); FP9 <- c(); FN9 <- c(); TN9 <- c() # Initialize vectors for TP, FP, FN, TN under thereshold 0.9
  
  # Generate dataset
  data <- generate_data(n, p)
  
  x <- as.matrix(data[, 2:ncol(data)]) # Predictors
  y <- data[,1] # Response
  
  x <- scale(x)  # Standardize the predictors
  y <- scale(y, scale = FALSE) # Center the response
  colnames(x) <- paste0("V", 1:ncol(x)) # Set column names for predictors
  
  cv_fit <- cv.ncvreg(x, y, penalty = "MCP") # Cross-validation
  candidate_set <- cv_fit$lambda  # Candidate lambda values
  
  S_list <- vector("list", length(candidate_set)) # Initialize a list to store selection matrix for each lambda
  names(S_list) <- paste0("lambda_", seq_along(candidate_set)) # Name the list entries
  
  B <- 100 # Number of subsamples
  
  # Stability Selection for each lambda in candidate_set
  for (lambda_idx in seq_along(candidate_set)) {
    
    lambda <- candidate_set[lambda_idx]  # Current lambda value
    S <- matrix(0, nrow = B, ncol = p)  # Initialize selection matrix for the current lambda
    colnames(S) <- colnames(x) # Set column names of S to predictor names
    
    for (j in 1:B) { 
      set.seed(26 + i + j) # Set seed for reproducibility
      # sub-sampling the data (half of the original data without replacement)
      sample_index <- sample(1:nrow(data), nrow(data) / 2, replace = FALSE) 
      
      # Prepare the response and predictors
      x_sub <- x[sample_index, ] # Predictors
      y_sub <- y[sample_index] # Response
      
      # Fit MCP-penalized regression
      MCP_fit <- ncvreg(x_sub, y_sub, penalty = "MCP", lambda = lambda)
      
      # Extract significant predictors (ignoring the intercept, hence [-1])
      significant_predictors <- ifelse(coef(MCP_fit) != 0, 1, 0)[-1]
      
      # Store the significant predictors in matrix S
      S[j, ] <- significant_predictors
    }
    
    # Store the matrix S for the current lambda in the corresponding list entry
    S_list[[lambda_idx]] <- S
  }
  
  stability_results <- lapply(S_list, function(S) { # Loop through regularization set and compute stability values
    getStability(S) # Compute stability values
  })
  
  stab_values <- c() # Initialize an empty vector to store stability values of each lambda
  for (k in 1:length(candidate_set)) {
    temp <- stability_results[[paste0("lambda_", k)]]$stability # Extract stability values
    stab_values <- c(stab_values, temp) # Append stability values to the vector
  }
  
  # Finding lambda stable_1sd
  max_stability <- max(stab_values, na.rm = T) # Find the maximum stability value
  stability_1sd_threshold <- max_stability - sd(stab_values, na.rm = T) # Define the stability threshold as max stability - 1SD
  index_of_stable_1sd <- max(which(stab_values >= stability_1sd_threshold)) # since candidate values are sorted decreasingly, we take the last index to get the minimum value
  Stability_values <- c(Stability_values, stab_values[index_of_stable_1sd]) # store stability value
  
  important_vars <- paste0("V", seq(p / 5, p, by = p / 5))  # Last element of each group
  non_important_vars <- setdiff(paste0("V", 1:p), important_vars) # Non-important variables
  
  S_stable_1sd <- as.data.frame(S_list[[index_of_stable_1sd]])  # Extract the selection matrix
  col_means <- colMeans(S_stable_1sd) # Calculate column means of the selection matrix
  # Filter predictors with selection frequency > 0.6
  S_stable_1sd_filtered <- S_stable_1sd[, col_means > 0.6, drop = FALSE]
  # Get selected variable names
  selected_vars <- names(round(colMeans(S_stable_1sd_filtered), 3))
  # True Positives (TP): Important variables selected
  tp_count <- sum(selected_vars %in% important_vars)
  TP6 <- c(TP6, tp_count)
  # False Positives (FP): Non-important variables selected
  fp_count <- sum(selected_vars %in% non_important_vars)
  FP6 <- c(FP6, fp_count)
  # False Negatives (FN): Important variables not selected
  fn_count <- sum(!(important_vars %in% selected_vars))
  FN6 <- c(FN6, fn_count)
  # True Negatives (TN): Non-important variables not selected
  tn_count <- sum(!(non_important_vars %in% selected_vars))
  TN6 <- c(TN6, tn_count)
  # Filter predictors with selection frequency > 0.7
  S_stable_1sd_filtered <- S_stable_1sd[, col_means > 0.7, drop = FALSE]
  # Get selected variable names
  selected_vars <- names(round(colMeans(S_stable_1sd_filtered), 3))
  # True Positives (TP): Important variables selected
  tp_count <- sum(selected_vars %in% important_vars)
  TP7 <- c(TP7, tp_count)
  # False Positives (FP): Non-important variables selected
  fp_count <- sum(selected_vars %in% non_important_vars)
  FP7 <- c(FP7, fp_count)
  # False Negatives (FN): Important variables not selected
  fn_count <- sum(!(important_vars %in% selected_vars))
  FN7 <- c(FN7, fn_count)
  # True Negatives (TN): Non-important variables not selected
  tn_count <- sum(!(non_important_vars %in% selected_vars))
  TN7 <- c(TN7, tn_count)
  # Filter predictors with selection frequency > 0.8
  S_stable_1sd_filtered <- S_stable_1sd[, col_means > 0.8, drop = FALSE]
  # Get selected variable names
  selected_vars <- names(round(colMeans(S_stable_1sd_filtered), 3))
  # True Positives (TP): Important variables selected
  tp_count <- sum(selected_vars %in% important_vars)
  TP8 <- c(TP8, tp_count)
  # False Positives (FP): Non-important variables selected
  fp_count <- sum(selected_vars %in% non_important_vars)
  FP8 <- c(FP8, fp_count)
  # False Negatives (FN): Important variables not selected
  fn_count <- sum(!(important_vars %in% selected_vars))
  FN8 <- c(FN8, fn_count)
  # True Negatives (TN): Non-important variables not selected
  tn_count <- sum(!(non_important_vars %in% selected_vars))
  TN8 <- c(TN8, tn_count)
  # Filter predictors with selection frequency > 0.9
  S_stable_1sd_filtered <- S_stable_1sd[, col_means > 0.9, drop = FALSE]
  # Get selected variable names
  selected_vars <- names(round(colMeans(S_stable_1sd_filtered), 3))
  # True Positives (TP): Important variables selected
  tp_count <- sum(selected_vars %in% important_vars)
  TP9 <- c(TP9, tp_count)
  # False Positives (FP): Non-important variables selected
  fp_count <- sum(selected_vars %in% non_important_vars)
  FP9 <- c(FP9, fp_count)
  # False Negatives (FN): Important variables not selected
  fn_count <- sum(!(important_vars %in% selected_vars))
  FN9 <- c(FN9, fn_count)
  # True Negatives (TN): Non-important variables not selected
  tn_count <- sum(!(non_important_vars %in% selected_vars))
  TN9 <- c(TN9, tn_count)
  
  list(TP6 = TP6, FP6 = FP6, FN6 = FN6, TN6 = TN6,
       TP7 = TP7, FP7 = FP7, FN7 = FN7, TN7 = TN7,
       TP8 = TP8, FP8 = FP8, FN8 = FN8, TN8 = TN8,
       TP9 = TP9, FP9 = FP9, FN9 = FN9, TN9 = TN9,
       Stability_values = Stability_values)}
  , mc.cores = detectCores() - 1)})
print(time_taken[3]/3600) # Print the time taken for the computation

TP6 <- as.vector(do.call(rbind, lapply(results, `[[`, "TP6")))
TP7 <- as.vector(do.call(rbind, lapply(results, `[[`, "TP7")))
TP8 <- as.vector(do.call(rbind, lapply(results, `[[`, "TP8")))
TP9 <- as.vector(do.call(rbind, lapply(results, `[[`, "TP9")))
FP6 <- as.vector(do.call(rbind, lapply(results, `[[`, "FP6")))
FP7 <- as.vector(do.call(rbind, lapply(results, `[[`, "FP7")))
FP8 <- as.vector(do.call(rbind, lapply(results, `[[`, "FP8")))
FP9 <- as.vector(do.call(rbind, lapply(results, `[[`, "FP9")))
FN6 <- as.vector(do.call(rbind, lapply(results, `[[`, "FN6")))
FN7 <- as.vector(do.call(rbind, lapply(results, `[[`, "FN7")))
FN8 <- as.vector(do.call(rbind, lapply(results, `[[`, "FN8")))
FN9 <- as.vector(do.call(rbind, lapply(results, `[[`, "FN9")))
TN6 <- as.vector(do.call(rbind, lapply(results, `[[`, "TN6")))
TN7 <- as.vector(do.call(rbind, lapply(results, `[[`, "TN7")))
TN8 <- as.vector(do.call(rbind, lapply(results, `[[`, "TN8")))
TN9 <- as.vector(do.call(rbind, lapply(results, `[[`, "TN9")))
Stability_values <- as.vector(do.call(rbind, lapply(results, `[[`, "Stability_values")))
# Remove unnecessary objects and run the garbage collector
rm(results, time_taken); gc() 

# Stability Selection with Stable MCP
time_taken <- system.time({results <- mclapply(1:100, function(i) {
  set.seed(26 + i) # Set seed for reproducibility
  
  Stability_values2 <- c() # Initialize an empty vector to store stability values
  TP62 <- c(); FP62 <- c(); FN62 <- c(); TN62 <- c() # Initialize vectors for TP, FP, FN, TN under thereshold 0.6
  TP72 <- c(); FP72 <- c(); FN72 <- c(); TN72 <- c() # Initialize vectors for TP, FP, FN, TN under thereshold 0.7
  TP82 <- c(); FP82 <- c(); FN82 <- c(); TN82 <- c() # Initialize vectors for TP, FP, FN, TN under thereshold 0.8
  TP92 <- c(); FP92 <- c(); FN92 <- c(); TN92 <- c() # Initialize vectors for TP, FP, FN, TN under thereshold 0.9
  
  # Generate dataset
  data <- generate_data(n, p)
  
  x <- as.matrix(data[, 2:ncol(data)]) # Predictors
  y <- data[,1] # Response
  
  x <- scale(x)  # Standardize the predictors
  y <- scale(y, scale = FALSE) # Center the response
  colnames(x) <- paste0("V", 1:ncol(x)) # Set column names for predictors
  
  Threshold <- nrow(x)/log(nrow(x))  # Screening threshold
  AHOLP <- AirHOLP(x, y, Threshold = Threshold, r0 = 10, adapt = TRUE, iter = 10)
  ranked_features <- AHOLP$index_r  # Ranking of features
  w <- 1 - (1/ranked_features)
  
  cv_fit <- cv.ncvreg(x, y, penalty = "MCP") # Cross-validation
  candidate_set <- cv_fit$lambda  # Candidate lambda values
  
  S_list <- vector("list", length(candidate_set)) # Initialize a list to store selection matrix for each lambda
  names(S_list) <- paste0("lambda_", seq_along(candidate_set)) # Name the list entries
  
  B <- 100 # Number of subsamples
  
  # Stability Selection for each lambda in candidate_set
  for (lambda_idx in seq_along(candidate_set)) {
    
    lambda <- candidate_set[lambda_idx]  # Current lambda value
    S <- matrix(0, nrow = B, ncol = p)  # Initialize selection matrix for the current lambda
    colnames(S) <- colnames(x) # Set column names of S to predictor names
    
    for (j in 1:B) { 
      set.seed(26 + i + j) # Set seed for reproducibility
      # sub-sampling the data (half of the original data without replacement)
      sample_index <- sample(1:nrow(data), nrow(data) / 2, replace = FALSE) 
      
      # Prepare the response and predictors
      x_sub <- x[sample_index, ] # Predictors
      y_sub <- y[sample_index] # Response
      
      # Fit MCP-penalized regression
      MCP_fit <- ncvreg(x_sub, y_sub, penalty = "MCP", lambda = lambda, penalty.factor = w)
      
      # Extract significant predictors (ignoring the intercept, hence [-1])
      significant_predictors <- ifelse(coef(MCP_fit) != 0, 1, 0)[-1]
      
      # Store the significant predictors in matrix S
      S[j, ] <- significant_predictors
    }
    
    # Store the matrix S for the current lambda in the corresponding list entry
    S_list[[lambda_idx]] <- S
  }
  
  stability_results <- lapply(S_list, function(S) { # Loop through regularization set and compute stability values
    getStability(S) # Compute stability values
  })
  
  stab_values <- c() # Initialize an empty vector to store stability values of each lambda
  for (k in 1:length(candidate_set)) {
    temp <- stability_results[[paste0("lambda_", k)]]$stability # Extract stability values
    stab_values <- c(stab_values, temp) # Append stability values to the vector
  }
  
  # Finding lambda stable_1sd
  max_stability <- max(stab_values, na.rm = T) # Find the maximum stability value
  stability_1sd_threshold <- max_stability - sd(stab_values, na.rm = T) # Define the stability threshold as max stability - 1SD
  index_of_stable_1sd <- max(which(stab_values >= stability_1sd_threshold)) # since candidate values are sorted decreasingly, we take the last index to get the minimum value
  Stability_values2 <- c(Stability_values2, stab_values[index_of_stable_1sd]) # store stability value
  
  important_vars <- paste0("V", seq(p / 5, p, by = p / 5))  # Last element of each group
  non_important_vars <- setdiff(paste0("V", 1:p), important_vars) # Non-important variables
  S_stable_1sd <- as.data.frame(S_list[[index_of_stable_1sd]])  # Extract the selection matrix
  col_means <- colMeans(S_stable_1sd) # Calculate column means of the selection matrix
  # Filter predictors with selection frequency > 0.6
  S_stable_1sd_filtered <- S_stable_1sd[, col_means > 0.6, drop = FALSE]
  # Get selected variable names
  selected_vars <- names(round(colMeans(S_stable_1sd_filtered), 3))
  # True Positives (TP): Important variables selected
  tp_count <- sum(selected_vars %in% important_vars)
  TP62 <- c(TP62, tp_count)
  # False Positives (FP): Non-important variables selected
  fp_count <- sum(selected_vars %in% non_important_vars)
  FP62 <- c(FP62, fp_count)
  # False Negatives (FN): Important variables not selected
  fn_count <- sum(!(important_vars %in% selected_vars))
  FN62 <- c(FN62, fn_count)
  # True Negatives (TN): Non-important variables not selected
  tn_count <- sum(!(non_important_vars %in% selected_vars))
  TN62 <- c(TN62, tn_count)
  # Filter predictors with selection frequency > 0.7
  S_stable_1sd_filtered <- S_stable_1sd[, col_means > 0.7, drop = FALSE]
  # Get selected variable names
  selected_vars <- names(round(colMeans(S_stable_1sd_filtered), 3))
  # True Positives (TP): Important variables selected
  tp_count <- sum(selected_vars %in% important_vars)
  TP72 <- c(TP72, tp_count)
  # False Positives (FP): Non-important variables selected
  fp_count <- sum(selected_vars %in% non_important_vars)
  FP72 <- c(FP72, fp_count)
  # False Negatives (FN): Important variables not selected
  fn_count <- sum(!(important_vars %in% selected_vars))
  FN72 <- c(FN72, fn_count)
  # True Negatives (TN): Non-important variables not selected
  tn_count <- sum(!(non_important_vars %in% selected_vars))
  TN72 <- c(TN72, tn_count)
  # Filter predictors with selection frequency > 0.8
  S_stable_1sd_filtered <- S_stable_1sd[, col_means > 0.8, drop = FALSE]
  # Get selected variable names
  selected_vars <- names(round(colMeans(S_stable_1sd_filtered), 3))
  # True Positives (TP): Important variables selected
  tp_count <- sum(selected_vars %in% important_vars)
  TP82 <- c(TP82, tp_count)
  # False Positives (FP): Non-important variables selected
  fp_count <- sum(selected_vars %in% non_important_vars)
  FP82 <- c(FP82, fp_count)
  # False Negatives (FN): Important variables not selected
  fn_count <- sum(!(important_vars %in% selected_vars))
  FN82 <- c(FN82, fn_count)
  # True Negatives (TN): Non-important variables not selected
  tn_count <- sum(!(non_important_vars %in% selected_vars))
  TN82 <- c(TN82, tn_count)
  # Filter predictors with selection frequency > 0.9
  S_stable_1sd_filtered <- S_stable_1sd[, col_means > 0.9, drop = FALSE]
  # Get selected variable names
  selected_vars <- names(round(colMeans(S_stable_1sd_filtered), 3))
  # True Positives (TP): Important variables selected
  tp_count <- sum(selected_vars %in% important_vars)
  TP92 <- c(TP92, tp_count)
  # False Positives (FP): Non-important variables selected
  fp_count <- sum(selected_vars %in% non_important_vars)
  FP92 <- c(FP92, fp_count)
  # False Negatives (FN): Important variables not selected
  fn_count <- sum(!(important_vars %in% selected_vars))
  FN92 <- c(FN92, fn_count)
  # True Negatives (TN): Non-important variables not selected
  tn_count <- sum(!(non_important_vars %in% selected_vars))
  TN92 <- c(TN92, tn_count)
  
  list(TP62 = TP62, FP62 = FP62, FN62 = FN62, TN62 = TN62,
       TP72 = TP72, FP72 = FP72, FN72 = FN72, TN72 = TN72,
       TP82 = TP82, FP82 = FP82, FN82 = FN82, TN82 = TN82,
       TP92 = TP92, FP92 = FP92, FN92 = FN92, TN92 = TN92,
       Stability_values2 = Stability_values2)}
  , mc.cores = detectCores() - 1)})
print(time_taken[3]/3600) # Print the time taken for the computation

TP62 <- as.vector(do.call(rbind, lapply(results, `[[`, "TP62")))
TP72 <- as.vector(do.call(rbind, lapply(results, `[[`, "TP72")))
TP82 <- as.vector(do.call(rbind, lapply(results, `[[`, "TP82")))
TP92 <- as.vector(do.call(rbind, lapply(results, `[[`, "TP92")))
FP62 <- as.vector(do.call(rbind, lapply(results, `[[`, "FP62")))
FP72 <- as.vector(do.call(rbind, lapply(results, `[[`, "FP72")))
FP82 <- as.vector(do.call(rbind, lapply(results, `[[`, "FP82")))
FP92 <- as.vector(do.call(rbind, lapply(results, `[[`, "FP92")))
FN62 <- as.vector(do.call(rbind, lapply(results, `[[`, "FN62")))
FN72 <- as.vector(do.call(rbind, lapply(results, `[[`, "FN72")))
FN82 <- as.vector(do.call(rbind, lapply(results, `[[`, "FN82")))
FN92 <- as.vector(do.call(rbind, lapply(results, `[[`, "FN92")))
TN62 <- as.vector(do.call(rbind, lapply(results, `[[`, "TN62")))
TN72 <- as.vector(do.call(rbind, lapply(results, `[[`, "TN72")))
TN82 <- as.vector(do.call(rbind, lapply(results, `[[`, "TN82")))
TN92 <- as.vector(do.call(rbind, lapply(results, `[[`, "TN92")))
Stability_values2 <- as.vector(do.call(rbind, lapply(results, `[[`, "Stability_values2")))
# Remove unnecessary objects and run the garbage collector
rm(results, time_taken); gc() 




# Combine results into a data frame
results <- data.frame(
  'FP0.6' = FP6, 'TP0.6' = TP6, 'FN0.6' = FN6, 'TN0.6' = TN6, 
  'FP0.7' = FP7, 'TP0.7' = TP7, 'FN0.7' = FN7, 'TN0.7' = TN7, 
  'FP0.8' = FP8, 'TP0.8' = TP8, 'FN0.8' = FN8, 'TN0.8' = TN8, 
  'FP0.9' = FP9, 'TP0.9' = TP9, 'FN0.9' = FN9, 'TN0.9' = TN9, 
  'FP0.62' = FP62, 'TP0.62' = TP62, 'FN0.62' = FN62, 'TN0.62' = TN62, 
  'FP0.72' = FP72, 'TP0.72' = TP72, 'FN0.72' = FN72, 'TN0.72' = TN72, 
  'FP0.82' = FP82, 'TP0.82' = TP82, 'FN0.82' = FN82, 'TN0.82' = TN82, 
  'FP0.92' = FP92, 'TP0.92' = TP92, 'FN0.92' = FN92, 'TN0.92' = TN92,
  'S1' = Stability_values, 'S2' = Stability_values2
)

# Function to calculate F1-score
calculate_F1 <- function(results, thresholds) {
  F1 <- data.frame(Threshold = character(), F1 = numeric(), stringsAsFactors = FALSE)
  
  for (t in thresholds) {
    fp <- results[[paste0("FP0.", t)]]
    tp <- results[[paste0("TP0.", t)]]
    fn <- results[[paste0("FN0.", t)]]
    tn <- results[[paste0("TN0.", t)]]
    
    precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
    recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
    f1 <- ifelse(precision + recall == 0, 0, 2 * (precision * recall) / (precision + recall))
    
    F1 <- rbind(F1, data.frame(Threshold = t, F1 = f1))
  }
  
  return(F1)
}


tmp_df <- results[, c("S1", "S2")] # Extract stability values
tmp_df <- melt(tmp_df) # Reshape data for ggplot
colnames(tmp_df) <- c("Method", "Stability") # Rename columns
tmp_df$Method <- as.character(tmp_df$Method) # Convert to character
tmp_df$Method[tmp_df$Method == "S1"] <- "MCP" # Rename methods
tmp_df$Method[tmp_df$Method == "S2"] <- "Stable MCP" # Rename methods


# plot stability values across methods
ggplot(tmp_df, aes(x = Method, y = Stability, fill = Method)) +
  geom_boxplot() +
  scale_fill_manual(values = friendly_pal("ito_seven")) + # Use the color-blind palette
  labs(title = "Stability Comparison", x = "Method", y = "Stability") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),  # Center title and increase size
    axis.title = element_text(size = 16),  # Increase axis title size
    axis.text = element_text(size = 14, face = "bold"),   # Increase axis label size
    legend.title = element_text(size = 16),  # Increase legend title size
    legend.text = element_text(size = 14),    # Increase legend text size
    axis.text.x = element_text(angle = 90, hjust = 1)
  )
rm(tmp_df) # Remove temporary data frame


metrics_results1 <- calculate_F1(results, c(6, 62)) # Calculate F1-score for thresholds 0.6
metrics_results2 <- calculate_F1(results, c(7, 72)) # Calculate F1-score for thresholds 0.7
metrics_results3 <- calculate_F1(results, c(8, 82)) # Calculate F1-score for thresholds 0.8
metrics_results4 <- calculate_F1(results, c(9, 92)) # Calculate F1-score for thresholds 0.9
metrics_results <- rbind(metrics_results1, metrics_results2, metrics_results3, metrics_results4)
rm(metrics_results1, metrics_results2, metrics_results3, metrics_results4) # Remove temporary data frames


metrics_results[, 3] <- rep(c(1, 2), each = 100, times = 4) # Add method column
colnames(metrics_results) <- c("Threshold", "F1-score", "Method") # Rename columns

# Rename thresholds
metrics_results$Threshold[metrics_results$Threshold == 6] <- 0.6
metrics_results$Threshold[metrics_results$Threshold == 62] <- 0.6

metrics_results$Threshold[metrics_results$Threshold == 7] <- 0.7
metrics_results$Threshold[metrics_results$Threshold == 72] <- 0.7

metrics_results$Threshold[metrics_results$Threshold == 8] <- 0.8
metrics_results$Threshold[metrics_results$Threshold == 82] <- 0.8

metrics_results$Threshold[metrics_results$Threshold == 9] <- 0.9
metrics_results$Threshold[metrics_results$Threshold == 92] <- 0.9


metrics_results$Method <- as.character(metrics_results$Method) # Convert method column to character
metrics_results$Method[metrics_results$Method == 1] <- 'MCP' # Rename methods
metrics_results$Method[metrics_results$Method == 2] <- 'Stable MCP' # Rename methods


metrics_results$Method <- as.factor(metrics_results$Method) # Convert method column to factor
metrics_results$Threshold <- as.factor(metrics_results$Threshold) # Convert threshold column to factor

# Plot F1-score comparison
ggplot(metrics_results, aes(x = Threshold, y = `F1-score`, colour = Method, group = Method)) +
  stat_summary(fun = "mean", geom = "line", size = 1.2) +   # Line through mean points
  stat_summary(fun = "mean", geom = "point", size = 2.5) +  # Points at mean
  scale_colour_manual(values = friendly_pal("ito_seven")) + 
  labs(title = "F1-score Comparison", x = "Threshold", y = "F1-score Mean") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

write.csv(results, "MCP_Full.csv", row.names = FALSE) # Save results to CSV file

rm(list = ls()); gc() # Clear all objects from the environment and run garbage collector
