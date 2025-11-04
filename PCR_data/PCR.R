# This script reproduces the analysis on the PCR data
options(warn=-1) # Turn off warnings
if (!requireNamespace("glmnet")) {install.packages("glmnet")} # install package if not already installed
library(glmnet)
set.seed(26) # Set seed for reproducibility

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
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] glmnet_4.1-8 Matrix_1.7-1

#loaded via a namespace (and not attached):
#  [1] compiler_4.4.2    tools_4.4.2       survival_3.8-3    rstudioapi_0.17.1
#[5] Rcpp_1.1.0        splines_4.4.2     codetools_0.2-20  grid_4.4.2       
#[9] iterators_1.0.14  foreach_1.5.2     shape_1.4.6.1     lattice_0.22-6 


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


x <- read.table('/Users/48099783/Desktop/Thesis/My PhD/Stable_Lasso/PCR/Xgene.txt') # Load the predictor matrix
x <- t(as.matrix(x)) # Transpose the matrix to have genes as columns
x <- scale(x)  # Standardize the predictors
y <- unlist(read.table('/Users/48099783/Desktop/Thesis/My PhD/Stable_Lasso/PCR/Y3.txt')) # Load the response vector
y <- scale(y, scale = FALSE) # Center the response

gene_names <- unlist(read.table('/Users/48099783/Desktop/Thesis/My PhD/Stable_Lasso/PCR/gene_id.txt')) # Load gene names
colnames(x) <- gene_names # Assign gene names as column names of the predictor matrix
rm(gene_names) # Remove gene names to free up memory

n <- nrow(x) # Number of observations
p <- ncol(x) # Number of predictors


# Stability Selection with Lasso
cv_lasso <- cv.glmnet(x, y, nfolds = 10, alpha = 1) # Fit Lasso model with 10-fold CV
  
candidate_set <- cv_lasso$lambda # Candidate set of lambda values
  
S_list <- vector("list", length(candidate_set)) # Initialize a list to store selection matrix for each lambda
names(S_list) <- paste0("lambda_", seq_along(candidate_set)) # Name the list entries
  
B <- 100 # Number of subsamples
 
# Stability Selection for each lambda in candidate_set
for (lambda_idx in seq_along(candidate_set)) {
  
  lambda <- candidate_set[lambda_idx]  # Current lambda value
  S <- matrix(0, nrow = B, ncol = p)  # Initialize selection matrix for the current lambda
  colnames(S) <- colnames(x) # Set column names of S to predictor names
  
  for (j in 1:B) { 
    set.seed(26 + j) # Set seed for reproducibility
    # sub-sampling the data (half of the original data without replacement)
    sample_index <- sample(1:nrow(x), nrow(x) / 2, replace = FALSE) 
    
    # Prepare the response and predictors
    x_sub <- x[sample_index, ] # Predictors
    y_sub <- y[sample_index] # Response
    
    lasso_model <- glmnet(x_sub, y_sub, alpha = 1, lambda = lambda)
    # Extract significant predictors (ignoring the intercept, hence [-1])
    significant_predictors <- ifelse(coef(lasso_model) != 0, 1, 0)[-1]
    
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
lambda_stable_1sd <- candidate_set[index_of_stable_1sd] # lambda stable_1sd

S_stable_1sd <- as.data.frame(S_list[[index_of_stable_1sd]])  # Extract the selection matrix
col_means <- colMeans(S_stable_1sd) # Calculate column means of the selection matrix
# Filter predictors with selection frequency > 0.2
S_stable_1sd_filtered <- S_stable_1sd[, col_means > 0.2, drop = FALSE]

Lasso_results <- round(colMeans(S_stable_1sd_filtered), 3) # Print the selection frequencies of the filtered predictors

rm(lasso_model, S, S_list, S_stable_1sd, S_stable_1sd_filtered, stability_results,
   x_sub, col_means, index_of_stable_1sd, j, k, lambda, lambda_idx,
   max_stability, sample_index, significant_predictors, stab_values,
   stability_1sd_threshold, temp, y_sub, cv_lasso, lambda_stable_1sd); gc() # Remove unnecessary objects and run the garbage collector




# Stability Selection with Stable Lasso 
S_list <- vector("list", length(candidate_set)) # Initialize a list to store selection matrix for each lambda
names(S_list) <- paste0("lambda_", seq_along(candidate_set)) # Name the list entries

Threshold <- nrow(x)/log(nrow(x))  # Screening threshold
AHOLP <- AirHOLP(x, y, Threshold = Threshold, r0 = 10, adapt = TRUE, iter = 10)
ranked_features <- AHOLP$index_r  # Ranking of features
w <- 1 - (1/ranked_features) # define weights based on the ranking of features

# Stability Selection for each lambda in candidate_set
for (lambda_idx in seq_along(candidate_set)) {
    
  lambda <- candidate_set[lambda_idx]  # Current lambda value
  S <- matrix(0, nrow = B, ncol = p)  # Initialize selection matrix for the current lambda
  colnames(S) <- colnames(x) # Set column names of S to predictor names
    
  for (j in 1:B) { 
    set.seed(26 + j) # Set seed for reproducibility
    # sub-sampling the data (half of the original data without replacement)
    sample_index <- sample(1:nrow(x), nrow(x) / 2, replace = FALSE) 
      
    # Prepare the response and predictors
    x_sub <- x[sample_index, ] # Predictors
    y_sub <- y[sample_index] # Response
      
    lasso_model <- glmnet(x_sub, y_sub, alpha = 1, lambda = lambda, penalty.factor = w)
    # Extract significant predictors (ignoring the intercept, hence [-1])
    significant_predictors <- ifelse(coef(lasso_model) != 0, 1, 0)[-1]
      
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


# Finding lambda stable
max_stability <- max(stab_values, na.rm = T) # Find the maximum stability value
stable_values <- which(stab_values > 0.75) # Index of stable lambda values
lambda_stable <- min(candidate_set[stable_values]) # Minimum stable lambda value
index_of_lambda_stable <- which(candidate_set == lambda_stable) # Index of lambda_stable

S_stable <- S_list[[index_of_lambda_stable]] # Extract the selection matrix of lambda_stable
col_means <- colMeans(S_stable) # Calculate column means of the selection matrix
# Filter predictors with selection frequency > 0.2
S_stable_filtered <- S_stable[, col_means > 0.2, drop = FALSE]

Stable_Lasso_results <- round(colMeans(S_stable_filtered), 3) 



# Remove unnecessary objects and run the garbage collector
rm(lasso_model, S, S_list, S_stable_1sd, S_stable_1sd_filtered, stability_results,
   x_sub, col_means, j, k, lambda, lambda_idx, S_stable_filtered,
   max_stability, sample_index, significant_predictors, stab_values,
   temp, y_sub, Threshold, screening_vector, AHOLP, w, S_stable, B, candidate_set,
   index_of_lambda_stable, lambda_stable, lambda_stable_1sd, AirHOLP,
   newton, getStability, ranked_features, cv_lasso, stable_values, x, y, n, p); gc() # Remove unnecessary objects and run the garbage collector

print("Lasso Results:")
print(Lasso_results) 
print("Stable Lasso Results:")
print(Stable_Lasso_results) 
