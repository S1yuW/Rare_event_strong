library(ggplot2)
library(patchwork)

#################
beta <- 1
n <- 10
p1 <- 10^3
p2 <- 10^9

N1 <- 500
N2 <- 5000
gap <- 500

x <- c(1.00001,1.00002,,1.00003,,1.00004,,1.00005)

x <- c(seq(1.00001,1.01,0.005))
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

