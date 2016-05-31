# function to test for normality and plot histograms.
# outputs variables names and the result of the normality test
# TRUE = anderson darling test significant (p < 0.05)
norm.tests <- function (dat, varlist, plot = T){
  output <- array()
  for (i in varlist){
    #x <- shapiro.test(dat[,i]) ## Shapiro test of normality
    x <- ad.test(dat[,i])       ## anderson darling test of normality
    output <- c(output, x$p.value)
    
    if(plot){
      qqnorm(dat[,i], main = colnames(dat)[i])
      qqline(dat[,i])
      hist(dat[,i], plot = T, main = colnames(dat)[i])
    }  
  }
  output <- output[-1]
  return(cbind(colnames(dat)[varlist],output<0.05))
}


# Function for Graphical Assessment of Multivariate Normality
graph.mvn <- function(dat, vars, title = NA){
  x <- as.matrix(dat[,vars]) # n x p numeric matrix
  center <- colMeans(x) # centroid
  n <- nrow(x); p <- ncol(x); cov <- cov(x); 
  d <- mahalanobis(x,center,cov) # distances 
  qqplot(qchisq(ppoints(n),df=p),d,
         main="QQ Plot Assessing Multivariate Normality", sub = title,
         ylab="Mahalanobis D2")
  abline(a=0,b=1) 
}  


# function for Graphical Assessment of Multivariate Normality by Stock
graph.mvn.Stock <- function(dat, vars){
  for (i in unique(dat$Stock)){
    graph.mvn(dat[dat$Stock ==i,], vars, title = i)
  }
}
