library(ggplot2)

library(ggpubr)
#Given beta,n,p1,p2 and different x, the transformation curve of C.O.V. with respect to N
#  N <- c(seq(N1,N2,gap))

beta <- 1
n <- 10
p1 <- 1000
p2 <- 1e+09
N1 <- 100
N2 <- 2e+03
gap <- 100

x <- seq(1.2,1.4,0.05)
y <- array(0, c(length(x)*length(c(seq(N1,N2,gap))), 3))


for (i in 1:length(x)) {
  i1 <- (i-1)*length(c(seq(N1,N2,gap)))+1
  i2 <- i*length(c(seq(N1,N2,gap)))
  y[i1:i2,1] <- c(seq(N1,N2,gap))
  CNSim <- COVNSim(beta,n,p1, p2,x[i],N1, N2,gap)
  y[i1:i2,2]  <- CNSim[,2]
  y[i1:i2,3] <- x[i]
}

data <- data.frame(y)
colnames(data) <- c('N','C.O.V.','x')

data$x <- as.character(data$x)

#cols <- c(rgb(251,132,2, maxColorValue = 255),rgb(255,195,0, maxColorValue = 255),rgb(1,53,101, maxColorValue = 255),rgb(1,36,76, maxColorValue = 255),rgb(1,7,19, maxColorValue = 255))           

#cols <- c(rgb(251,132,2, maxColorValue = 255),rgb(255,195,0, maxColorValue = 255),rgb(40,114,113, maxColorValue = 255),rgb(38,70,83, maxColorValue = 255),rgb(1,36,76, maxColorValue = 255))           

cols <- c(rgb(251,132,2, maxColorValue = 255),rgb(255,195,0, maxColorValue = 255),rgb(1,86,153, maxColorValue = 255),rgb(95,198,201, maxColorValue = 255),rgb(79,89,100, maxColorValue = 255))           

p <- ggplot(data, aes(N,C.O.V.))  +
  geom_point(aes(color = x, shape= x),size = 2) + 
  theme_bw()  +
  scale_color_manual(values = cols)
p

ggsave('COVN.png', plot = p, width = 15, height = 5, dpi = 600 )
##############################

#################
beta <- 2
n <- 10
p1 <- 20
p2 <- 40
x <- c(seq(0.75,0.85,0.001))
iter_num <- 5000
y <- array(0, c(length(x),2 ) )
for (i in 1: length(x) ){
  
  y[i,1] <- x[i]
  a <- ((p1 + p2) * x[i] - p1) / sqrt(n * p1)
  
  LL <- CopiesIS(beta,n,p1, p2,a,iter_num)
  y[i,2]  <- mean(LL)
}
data <- data.frame(y)
colnames(data) <- c('x','logP')


p<-ggplot(data, aes(x, logP)) +
  geom_point(size = 0.6) +
  theme_bw() +
  scale_color_manual(values = cols) +
  scale_y_log10()


ggsave('logP.png', plot = p, width = 5, height = 5, dpi = 600 )
##############################

