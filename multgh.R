multgh <- function(N=30,g,h,Sigma,warnings = T){ 
  # Parameters:
  ############ N :                      the number of datapoints to generate
  ############ h:                a vector of h values for the simulated data
  ############ kurtosis:                a vector of kurtosis values for the simulated data
  ############ Sigma:                   the specified correlation matrix 
  ############ warnings:                set to F to ignore warnings
  
  library('rootSolve') ## For implementation of Newton's method for roots approximation.
  library('lattice')   ## For making a grid of points to look for convenient guesses to supply in the Newton's method
  library('pracma')    ## For numerical integration
  library('Matrix')    ## For matrix operations
  
  
  model <- function(x, parms = c(0,0)){
    ### The formula to calculate the skewness
    c(F1 = ((3*exp(x[1]^2 / (2-6*x[2])) + exp(9*x[1]^2 / (2-6*x[2])) - 
            3*exp(2*x[1]^2 /(1 - 3*x[2])) - 1) / (1-3*x[2])^0.5 - 3*(1 - 2*exp(x[1]^2 /(2-4*x[2]))+
            exp(2*x[1]^2/(1-2*x[2])))*(exp(x[1]^2 / (2-2*x[2])) -1)/((1-2*x[2])^0.5 * (1-x[2])^0.5) +
            2*(exp(x[1]^2/(2-2*x[2])) -1)^3 / (1-x[2])^1.5) / (x[1]^3 * (((1-2*exp(x[1]^2 /(2-4*x[2])) +
            exp(2*x[1]^2 / (1-2*x[2]))) / (1-2*x[2])^0.5 + (exp(x[1]^2/(2-2*x[2])) -1)^2 / (x[2] -1))/x[1]^2)^1.5) -
            parms[1],
    
    ### The formula to calculate the kurtosis
    F2 = (exp(8*x[1]^2/(1-4*x[2]))*(1 + 6*exp(6*x[1]^2/(4*x[2]-1)) + exp(8*x[1]^2/(4*x[2]-1)) -
          4*exp(7*x[1]^2/(8*x[2]-2)) - 4*exp(15*x[1]^2/(8*x[2]-2)))/(1-4*x[2])^0.5 -
            4*(3*exp(x[1]^2/(2-6*x[2])) + exp(9*x[1]^2/(2-6*x[2])) - 3*exp(2*x[1]^2/(1-3*x[2]))-1)*
            (exp(x[1]^2/(2-2*x[2]))-1)/((1-3*x[2])^0.5 * (1-x[2])^0.5) -
            6*(exp(x[1]^2/(2-2*x[2]))-1)^4/(x[2]-1)^2 -12*(1-2*exp(x[1]^2/(2-4*x[2])) +
              exp(2*x[1]^2/(1-2*x[2])))*(exp(x[1]^2/(2-2*x[2]))-1)^2/((1-2*x[2])^0.5*(x[2]-1))+
            3*(1-2*exp(x[1]^2/(2-4*x[2])) + exp(2*x[1]^2/(1-2*x[2])))^2/(2*x[2]-1))/
      ((1-2*exp(x[1]^2/(2-4*x[2])) + exp(2*x[1]^2/(1-2*x[2])))/(1-2*x[2])^0.5 +
         (exp(x[1]^2/(2-2*x[2]))-1)^2/(x[2]-1))^2 - 
      parms[2])}
  
  # The function to calculate the mean given g and h: 
  
  mu_gh <- function(g,h) (exp(g^2 *(2-2*h)^-1)-1)/(g*(1-h)^0.5)
  
  # The function to calculate the variance given g and h. 
  
  variance_gh <- function(g,h) ((1-2*exp(g^2 *(2 - 4*h)^-1) + exp(2*g^2 * (
    1-2*h)^-1))/(1-2*h)^0.5 + (exp(g^2 *(2 -2*h)^-1)-1)^2 / (h-1))/g^2
  
  # the function to generate g-h distribution.
  q_z <- function(g,h,z) (exp(g*z) - 1)*exp(0.5*h*z^2)/g 
  
  mu <- mu_gh(g,h)
  sd_s <- sqrt(variance_gh(g,h))
  
  # Now we'll compute the pairwise  intermediate correlations
  
  Sigma_intermediate <- matrix(nrow=length(h),ncol=length(h))
  for (i in 1:length(g)){
    for (j in (1:length(g))){
      if (i > j) { 
        ## we'll compare the pairwise correlations and construct the entire intermediate matrix
        q_i <- function(z) (q_z(g[i],h[i],z) - mu[i]) / sd_s[i] 
        q_j <- function(z) (q_z(g[j],h[j],z) - mu[j]) / sd_s[j]
        f_ij <- function(z_i,z_j,pij) (1/(2*pi*sqrt(1-pij^2))) * exp(-(2*(1-pij^2))^-1 *(z_i^2-2*pij*z_i*z_j +z_j^2))
        
        ## The integrand of the specified correlation formula. Note that this function's output is another function. 
        ## This is needed to calculate the integral using the integral2 function
        
        p_ij <- function(z_i,z_j,pij) function(z_i, z_j) q_i(z_i)*q_j(z_j)*f_ij(z_i,z_j,pij)

        ## The solver of this library implements a different algorithm (Gauss-Kronrod) than the paper's one (Adaptative Monte Carlo) 
    
        
        
        
        ### We need to calculate the intermediate correlation (whereas in Table 1 it is given) 
        ### by finding the root of the following function
        ### (where x is the intermediate correlation and parms is the specified correlation)
        
        #Z <- seq(-8,8,length.out = 10000)
        #D_i <-   sapply(Z,q_i) 
        #D_j <-   sapply(Z, q_j) 

        #xmin <- quantile(D_i,0.0)
        #xmax <- quantile(D_i,1)
        #ymin <- quantile(D_j,0.0)
        #ymax <- quantile(D_j,1)
        
        xmin <- -8
        xmax <- 8
        ymin <- -8
        ymax <- 8
        
        F_ij <- function(x,parms) (integral2(fun = p_ij(z_i,z_j,pij=x[1]),
                                    xmin = xmin,xmax = xmax, ymin= ymin, ymax = ymax)$Q - parms[1]) 
        
        
        ## Once again, we will use the same methods above to solve for a root of F_ij, and hence, find the intermediate correlation
        
        
        #x0 <- optimize(f = function(x) F_ij(x,Sigma[i,j]),lower = -0.99,upper = 0.99)$minimum
        
        #intermediate_corr <- uniroot(f = function(x) F_ij(x,Sigma[i,j]),lower = -0.99,upper = 0.99)
        #Sigma_intermediate[i,j] <- intermediate_corr$root
        
        intermediate_corr <- optimize(f = function(x) F_ij(x,Sigma[i,j])^2,lower = -0.99,upper = 0.99)
        Sigma_intermediate[i,j] <- intermediate_corr$minimum
        
      
      }
    }
  }
  
  diag(Sigma_intermediate) <- 1 # set the diagonal to ones
  Sigma_intermediate <- forceSymmetric(Sigma_intermediate,uplo="L") # as the correlation matrix is a sym. matrix.
  
  if (mean(eigen(x = Sigma_intermediate)$values>0) < 1){ # A way to evaluate if at least one eigenvalue is non-positive
    Sigma_intermediate <- nearPD(Sigma_intermediate)
    if (isTRUE(warning)) print('Warning!!! Non-positive definite matrix found. A correction was made to make it Positive Definite')
  }
  
  chol_decomposition <- chol(Sigma_intermediate) # Cholesky Decomposition
  V <- matrix(rnorm(length(h)*N),nrow=N,ncol=length(h)) # A matrix of standard normal values
  Z <- V %*% t(chol_decomposition)          # To apply formula (10) 
  output <- matrix(nrow = nrow(Z),ncol = ncol(Z))
  
  for (i in 1:ncol(Z)) {
    output[,i] <- q_z(g[i], h[i],Z[,i])
  }
  
  return (list(data = output, intermediate_matrix =Sigma_intermediate, chol = chol_decomposition, g=g,h=h))

  }