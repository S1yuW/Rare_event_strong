###################################################
# Let (\lambda_1, \lambda_2, \cdots, \lambda_n\right) be the beta-Jacobi ensemble J_n(p_1, p_2).
# Let \lambda_{(1)}< \cdots <\lambda_{(n)} be the order statistics of \lambda_1, \cdots, \lambda_n.
# In this R code, we construct strong efficient  importance Sampling estimator for P(p lambda_{(n)}/p_1 > x).
# Under the condition :n^5/p_1^3 \to 0 and n p_1^2 /p_2 \to 0
###################################################



########################

# Generate the tridiagonal matrix J. (The eigenvalues of J) ~ J_n(p_1, p_2).
FunJ <- function(beta, n, p1, p2){
  J <- array(0, c(n,n))
  c <- rep(0,n)
  s <- rep(0,n)
  for (i in 1:n){
    c[i] <- rbeta(1,beta*(p1-i+1)/2 , beta*(p2-i+1)/2 )
    s[i] <- rbeta(1,beta*(n-i)/2 , beta*(p1+p2 -n -i+1)/2 )
  }
  
  for (i in 2:n){
    J[i,i] <-s[i-1]*(1-c[i-1]) + c[i]*(1-s[i-1])
  }
  J[1,1] <- c[1]
  for (i in 3:n){
    J[i-1,i] <- sqrt(c[i-1]*(1- c[i-1])*s[i-1]*(1-s[i-2]))
    J[i,i-1] <- J[i-1,i] 
  }
  J[1,2] <- sqrt(c[1]*(1- c[1])*s[1]) 
  J[2,1] <- J[1,2] 
  return(J)
}


#Calculate Fn
Ffun <- function(beta, n, p1 ,p2, x){
  J <- FunJ(beta, n-1,p1-1,p2-1)
  eigenvalues <- eigen(J)$values
  ln1 <-eigenvalues[1]
  
  r <-  beta*(p1+p2)*(x-1)/2/x
  
  while (TRUE) {
     ln <- rexp(1,rate = r) +  max(ln1,x*p1/(p1+p2))
    if (ln < 1) {
      break
    }
  }
  
  r1n <- beta*(p1-n+1)/2
  r2n <- beta*(p2-n+1)/2
  logAn <- lgamma(1+beta/2)+ lgamma(beta*(p1+p2)/2)+ lgamma(beta*(p1+p2-1)/2) -lgamma(1+beta*n/2)- lgamma(beta*p1/2) - lgamma(beta*p2/2) - lgamma(beta*(p1+p2-n)/2)
  logunx <- (r1n-1)*log(ln) + (r2n-1)*log(1-ln)
  loghnx <- log(r) -log(1-exp(-r*(1- max(ln1,x*p1/(p1+p2))  ))) - r*ln + r*max(ln1,x*p1/(p1+p2))
  logprob <- log(n) + logAn + sum(beta*log(ln-eigenvalues)) + logunx - loghnx

  return(exp(logprob))

}


# Simulate iter_num i.i.d. copies using Ffun -> LL[iter]
CopiesIS <- function(beta, n, p1, p2, x, iter_num) {
  LL <- numeric(iter_num)
  
  for (iter in 1:iter_num) {
    if (is.na(x)) {
      # Handle case where x is NA
      break
    }
    
    LL[iter] <- Ffun(beta, n, p1, p2, x)
  }
  
  return(LL)
}


####### ggplot
#Given beta,n,p1,p2 and different x, the transformation curve of C.O.V. with respect to N
#  N <- c(seq(N1,N2,gap))
COVNSim <- function(beta, n, p1, p2, x, N1, N2, gap) {
  Nseq <- seq(N1, N2, gap)
  CN <- matrix(0, nrow = length(Nseq), ncol = 2)
  
  for (i in 1:length(Nseq)) {
    iter_num <- Nseq[i]
    CN[i, 1] <- iter_num
    LL <- CopiesIS(beta, n, p1, p2, x, iter_num)
    CN[i, 2] <- sqrt(mean(LL^2) - mean(LL)^2) / sqrt(iter_num) / mean(LL)
  }
  
  return(CN)
}


library(ggplot2)
library(patchwork)

############################################
# 这组参数设置满足条件H2，这个时候估计器是渐近有效的。
#H2 beta <- 100 p2 <- 10^9
#H1 beta <- 10 p2 <- 10^9
#H0 beta <- 1 p2 <- 10^9

#H2 
beta <- 1 
p2 <- 10^6

#H1 beta <- 1 p2 <- 10^7

#H0 beta <- 1 p2 <- 10^9

n <- 10
p1 <- 10^3


x <- seq(1.15,1.45,0.001)

y <- matrix(0, nrow = length(x), ncol = 3)

N <- 10000

for (i in 1: length(x) ){
  y[i,1] <- x[i]
  LL <- CopiesIS(beta,n,p1, p2,x[i],N)
  y[i,2]  <- mean(LL)
  y[i,3] <- sqrt(mean(LL^2) - mean(LL)^2) / sqrt(N) / mean(LL)
}
data <- data.frame(y)
colnames(data) <- c('x','logP','C.O.V')




# 第一幅图
plot1 <- ggplot(data, aes(x, C.O.V)) +
  geom_point(size = 1) +
  theme_bw() +
  scale_color_manual(values = cols) +
  scale_y_log10()

# 第二幅图
plot2 <- ggplot(data, aes(x, logP)) +
  geom_point(size = 1) +
  theme_bw() +
  scale_color_manual(values = cols) +
  scale_y_log10()

# 将两个图形左右拼在一起
combined_plot <- plot1 + plot2 +
  plot_layout(ncol = 2)

# 显示拼接后的图形
print(combined_plot)



ggsave('log1PH2.png', plot = combined_plot, width = 15, height = 5, dpi = 600 )
##############################



