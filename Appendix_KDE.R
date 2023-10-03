# ///////////////////////////////// #
#     Kernel Density Estimation     #  ---- 
#               Unibo               #
#     Academic year 2022-2023       #
#         Ruggero Zustovi           #
# ///////////////////////////////// #

rm(list=ls())
graphics.off()

## Data ----
data(faithful)
help(faithful)

hist(faithful$eruptions)
hist(faithful$waiting)

data1<-faithful$waiting
summary(data1)
hist(data1,xlab="Waiting time (in mins)",freq=F)
plot(density(data1))

Colors <- c("#FF0000","#000000","#242424","#494949", "#6D6D6D", "#929292","#B6B6B6","#DBDBDB")

## Plots (using different bandwidths) ----
# KDE using Gaussian kernel
hvect<-c(1,2,5,10)
kde1<-density(data1,bw=hvect[1])
plot(kde1,xlab="Waiting time (in mins)",
     main="KDE with different bandwidth values", 
     col=Colors[1], lwd = 2.25)


# plot(density(data1), col='red')
for(i in 2:7){
  kde1<-density(data1,bw=hvect[i])
  lines(kde1,col=Colors[i], lwd = 1.5)
}
legend("topleft",paste("h=",hvect),col=Colors,lty=rep(1,7))


## KDE with Bumps ----
# Generate sample data
set.seed(123)
data <- rnorm(10, mean = 0)
n <- length(data)

data.grid <- seq(from = min(data) - 1, to = max(data) + 1, by = 0.001)

# Compute kernel density estimation
density_data <- density(data, bw = 0.3)
data.h <- density_data$bw

gauss <- function(x) 1 / sqrt(2 * pi) * exp(-(x^2) / 2)
data.bumps <- sapply(data, function(a) gauss((data.grid - a) / data.h) / (n * data.h))

# Plot the KDE with bumps
plot(density_data, main = "Kernel density estimation with Bumps", 
     xlab = "X", ylab = "Density", lwd = 1.5)
rug(data, side = 1, col = "red", lwd = 0.7, ticksize = 0.03)

# Add kernel shape on the rug
for (i in 1:ncol(data.bumps)) {
  lines(data.grid, data.bumps[, i], col = "green", lwd = 1)
}




