gh_roots <- function(skew, kurt, m = 10, tries = 100, no.digits = 2, method="Newton", positive = FALSE){
  ### A function that implements Newton's Homotopy Continuation to find the roots of a system of nonlinear equations.
  ### Parameters:
      ### skewness :          the skewness of the gh distribution
      ### kurtosis :          the kurtosis of the gh distribution
      ### m        :          number of divisions for the interval [0,1]
      ### tries    ;          number of points to use as x0 (initial guess)
  ### Output:
      ### A matrix containing all the roots found along with their associated values of skewness and kurtosis (rounded to 3 decimals)
  
  library('rootSolve')
  
  model <- function(x, parms = c(0,0)){
    ### The formula to calculate the skewness
    
    return(c(F1 = ((3*exp(x[1]^2 / (2-6*x[2])) + exp(9*x[1]^2 / (2-6*x[2])) - 
                      3*exp(2*x[1]^2 /(1 - 3*x[2])) - 1) / (1-3*x[2])^0.5 - 3*(1 - 2*exp(x[1]^2 /(2-4*x[2]))+
                                                                                 exp(2*x[1]^2/(1-2*x[2])))*(exp(x[1]^2 / (2-2*x[2])) -1)/((1-2*x[2])^0.5 * (1-x[2])^0.5) +
                     2*(exp(x[1]^2/(2-2*x[2])) -1)^3 / (1-x[2])^1.5) / (x[1]^3 * (((1-2*exp(x[1]^2 /(2-4*x[2])) +
                                                                                      exp(2*x[1]^2 / (1-2*x[2]))) / (1-2*x[2])^0.5 + (exp(x[1]^2/(2-2*x[2])) -1)^2 / (x[2] -1))/x[1]^2)^1.5) -
               parms[1],
             
             ### The formula to calculate the kurtosis
             F2 = (exp(8*x[1]^2/(1-4*x[2]))*(1 + 6*exp(6*x[1]^2/(4*x[2]-1)) + exp(8*x[1]^2/(4*x[2]-1)) -
                                               4*exp(7*x[1]^2/(8*x[2]-2)) - 4*exp(15*x[1]^2/(8*x[2]-2)))/(1-4*x[2])^0.5 -
                     4*(3*exp(x[1]^2/(2-6*x[2])) + exp(9*x[1]^2/(2-6*x[2])) - 3*exp(2*x[1]^2/(1-3*x[2]))-1)*
                     (exp(x[1]^2/(2-2*x[2]))-1)/((1-3*x[2])^0.5 * (1-x[2])^04) -
                     6*(exp(x[1]^2/(2-2*x[2]))-1)^4/(x[2]-1)^2 -12*(1-2*exp(x[1]^2/(2-4*x[2])) +
                                                                      exp(2*x[1]^2/(1-2*x[2])))*(exp(x[1]^2/(2-2*x[2]))-1)^2/((1-2*x[2])^0.5*(x[2]-1))+
                     3*(1-2*exp(x[1]^2/(2-4*x[2])) + exp(2*x[1]^2/(1-2*x[2])))^2/(2*x[2]-1))/
               ((1-2*exp(x[1]^2/(2-4*x[2])) + exp(2*x[1]^2/(1-2*x[2])))/(1-2*x[2])^0.5 +
                  (exp(x[1]^2/(2-2*x[2]))-1)^2/(x[2]-1))^2 - 
               parms[2]))
  }
  
  mu_gh <- function(g,h) (exp(g^2 *(2-2*h)^-1)-1)/(g*(1-h)^0.5)
  variance_gh <- function(g,h) ((1-2*exp(g^2 *(2 - 4*h)^-1) + exp(2*g^2 * (
    1-2*h)^-1))/(1-2*h)^0.5 + (exp(g^2 *(2 -2*h)^-1)-1)^2 / (h-1))/g^2
  omega <- function(x, parms) sum(model(x,parms)^2)
  
 
  if (method=="Newton") {
    # Newton's Homotopy Function
    homotopy = function(x,x0,t,parms) model(x,parms) - (1-t)*model(x0,parms)
  
    } else {
  # Fixed point homotopy
    homotopy = function(x,x0,t,parms) t*model(x,parms) + (1-t)*(x-x0)
  
    }
  
  homotopy.root = function(fn = homotopy, x0, m = 30, parms){
    delta = 1/m
    t=0
    for (k in 1:m){
      t = t + delta
      x0 = multiroot(function(x) homotopy(x,x0,t,parms),
                     start = x0,
                     positive = positive)$root
      #print(paste(t,"---",x0[1],",",x0[2],"---",model(x0)[1],",",model(x0)[2]))
      if (is.na(model(x0)[1])){
        break
      }
      
    }
    return(x0)
  }
  # Main
  hs = gs = solutions = list()
  g <- seq(0,10,length.out = 1000)
  h <- seq(0+1e-8, 0.25 - 1e-8,length.out = 200)
  X <- as.matrix(expand.grid(x=g,y=h))
  
  for (i in 1:length(skew)){
  skewness = skew[i]
  kurtosis = kurt[i]
  parms = c(skewness,kurtosis)
  om = apply(X,1,function(x) omega(x,parms))
  ind <- which(om %in% sort(unique(om))[1:tries])
  X_guess = X[ind,]
  sols = matrix(nrow = tries,ncol = 6)
  for (j in 1:tries){
    x0 = X_guess[j,]
    x0 = homotopy.root(homotopy,x0,m,parms)
    sols[j,] = c(x0,model(x0),mu_gh(x0[1],x0[2]),variance_gh(x0[1],x0[2]))
  }
  #sols = (unique(na.omit(round(sols,no.digits))))
  sols = na.omit(sols)
  dup = duplicated(round(sols,digits = no.digits))
  sols = matrix(sols[!dup,],ncol=6)
  
  if (length(sols) ==0) {
    print(paste("Failed to find solutions for a skewness of",skewness,"and a kurtosis of ",kurtosis))
  } else {
    print(paste(length(sols[,1]), "solution(s) found for a skewness of ",skewness,"and a kurtosis of ",kurtosis))
  }
  gs[[i]] = sols[,1]
  hs[[i]] = sols[,2]
  solutions[[i]] = sols
  }
  
  gs = t(expand.grid(gs))
  hs = t(expand.grid(hs))
  
  out = list(g = gs, h = hs, solutions = solutions )
  return (out)
}

