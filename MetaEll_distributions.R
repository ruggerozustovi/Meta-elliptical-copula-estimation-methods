# ///////////////////////////////// #
#        Elliptical Copulas         #  ---- 
#               Unibo               #
#     Academic year 2022-2023       #
#         Ruggero Zustovi           #
# ///////////////////////////////// #

rm(list=ls())
graphics.off()

# Libraries ----

library(ElliptCopulas)
library(copula)
library(mvtnorm)  
library(MASS)
library(ggplot2)
library(ggExtra)
library(ggpubr)
library(gridExtra)


# Simulation of an elliptical copula ----
# set.seed(1234) # (reproducibility)

## 1. Copula distributions ----
gaussian = normalCopula(param=0.7, dim = 2)
studT = tCopula(param=0.7,dim=2,df=2)
copula.vec = c(gaussian, studT)
names(copula.vec) = list("Gaussian", "Student-t")

X = list()
X.df = list()
gplot_list = list()
## 2. Observations ----
for (i in 1:length(copula.vec) ){
  myMvdX <- mvdc(copula=copula.vec[[i]], margins=c("exp", "exp"),
                 paramMargins=list(list(rate=4),
                                   list(rate=2)))
  X[[i]] <- rMvdc(1000,myMvdX)
  X.df[[i]] <- data.frame(X[[i]])
  p <- ggplot(X.df[[i]], aes(X1, X2)) +
    geom_point() +
    ggtitle(paste("copula = ", names(copula.vec)[i]))
  gplot_list[[i]] <- p
}
ggMarginal(gplot_list[[1]], type = "histogram")
ggMarginal(gplot_list[[2]], type = "histogram")


#mtext("Observations", side =3,line = -2, outer = TRUE,col='black')
ggMarginal(p, type = "histogram")




## 3. Samples ----
#### samples from copulas ----
par(mfrow=c(1,2))

for (i in 1:length(copula.vec) ){
  u <- rCopula(1000,copula.vec[[i]])
  plot(u,xlab="u1",ylab="u2")
  contour(copula.vec[[i]],pCopula,col='darkblue',
          xlab="u1",ylab="u2")
  mtext(paste("copula = ", names(copula.vec)[i]), side =3,line = -2, outer = TRUE,col='black')
}

## 4. Contour plots ----

# cdf of copula
par(mfrow=c(1,2))

#### 3d plot ----
for (i in 1:length(copula.vec) ){
  persp(copula.vec[[i]],pCopula,col='lightgreen',main=paste("copula = ", names(copula.vec)[i]))
}
mtext("cdf of copulas", side =3,line = -2, outer = TRUE,col='black')

#### contour plot ----
for (i in 1:length(copula.vec) ){
  contour(copula.vec[[i]],pCopula,col='lightgreen',main=paste("copula = ", names(copula.vec)[i]))
}
mtext("cdf of copulas", side =3,line = -2, outer = TRUE,col='black')

#### density ----
persp(gaussian,dCopula,col='lightgreen',main="Gaussian")
persp(studT,dCopula,col='lightgreen',main="Student-t")
mtext("density of copulas", side =3,line = -2, outer = TRUE,col='black')





## ANOTHER PLOT ----
par(mfrow=c(1,1))

gaussian2 = normalCopula(param=0.7, dim = 2)
myMvdY <- mvdc(copula= gaussian2, margins=c("norm", "norm"),
               paramMargins=list(list(mean = 0, sd = 4),
                                 list(mean = 1, sd = 2)))
Y <- rMvdc(1000,myMvdY)
Y.df <- data.frame(Y)
meta_p1 <- ggplot(Y.df, aes(X1, X2)) +
  geom_point(color = "darkblue", fill = "cornflowerblue", shape = 21) +
  # ggtitle("Observations")+
  theme(plot.title = element_text(hjust = 0.5))

meta_p1 = ggMarginal(meta_p1, type = "histogram",
                     col = "darkblue", 
                     fill = "cornflowerblue")

UY=cbind(pnorm(Y[,1],mean=0,sd=4), pnorm(Y[,2],mean=1,sd=2))
plot(UY)

UY.df = data.frame(UY, pCopula(UY, gaussian2))
colnames(UY.df) = c("u1", "u2", "z")

meta_rank1 = ggplot(Y.df, aes(X1, X2)) +
  geom_density_2d(aes(colour = after_stat(level)), linewidth = 1)+
  labs(x = "X1",
       y = "X2") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14))



den1 = ggplot(UY.df, aes(x = u1, y = u2)) +
  geom_density_2d(aes(colour = after_stat(level)), linewidth = 1)+
  labs(x = "u1",
       y = "u2") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14))


myMvdY2 <- mvdc(copula= gaussian2, margins=c("beta", "gamma"),
               paramMargins=list(list(shape1 = 3, shape2 = 4),
                                 list(shape = 2, rate = 4)))
Y2 <- rMvdc(1000,myMvdY2)
Y2.df <- data.frame(Y2)
meta_p2 <- ggplot(Y2.df, aes(X1, X2)) +
  geom_point(color = "darkblue", fill = "cornflowerblue", shape = 21) +
  # ggtitle("Observations")+
  xlab("Y1") +
  ylab("Y2") +
  theme(plot.title = element_text(hjust = 0.5))

meta_p2 = ggMarginal(meta_p2, type = "histogram",
                     col = "darkblue", 
                     fill = "cornflowerblue")

UY2=cbind(pbeta(Y2[,1],shape1 = 3, shape2 = 4), 
          pgamma(Y2[,2],shape = 2, rate = 4))
plot(UY2)

UY2.df = data.frame(UY2, dCopula(UY2, gaussian2))
colnames(UY2.df) = c("u1", "u2", "z")

meta_rank2 <- ggplot(Y2.df, aes(X1, X2)) +
  geom_density_2d(aes(colour = after_stat(level)), linewidth = 1)+
  labs(x = "Y1",
       y = "Y2") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14))

den2 = ggplot(UY2.df, aes(x = u1, y = u2)) +
  geom_density_2d(aes(colour = after_stat(level)), linewidth = 1)+
  labs(x = "u1",
       y = "u2") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14))


fig1 = ggarrange(meta_p1, meta_p2)

annotate_figure(fig1,
                top = text_grob( "Observations",
                                  size = 18))


fig2 = ggarrange(den1, den2, 
          common.legend = TRUE,
          legend = "bottom")

annotate_figure(fig2, top = text_grob("Contour", size = 16))

par(mfrow = c(1,2))
den3d <- kde2d(Y.df$X1, Y.df$X2)
persp(den3d, theta=-60, phi=20, box = T, main = "Density plot",
      xlab="X1", ylab="X2", zlab="Density")
contour(den3d, xlab="X1", ylab="X2", main = "Contour plot",
        xlim = c(-10, 10), ylim = c(-4, 6), labcex = 0.6)

 
den3d2 <- kde2d(Y2.df$X1, Y2.df$X2)
persp(den3d2, theta=-60, phi=20, xlab="Y1", ylab="Y2", zlab="Density",
      main = "Density plot")
contour(den3d2, xlab="Y1", ylab="Y2", main = "Contour plot",
        xlim = c(0,0.8), ylim = c(-0.2, 1.5), labcex = 0.6)




