# ///////////////////////////////// #
#      Case studies & Examples      #  ---- 
#               Unibo               #
#     Academic year 2022-2023       #
#         Ruggero Zustovi           #
# ///////////////////////////////// #

rm(list=ls())
graphics.off()

# Libraries ----
library(evir)
library(spearmanCI)
library(copula)
library(RNOmni)
library(plot3D)
library(RColorBrewer)
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(ggpubr)

# Simulation ----
set.seed(1234)
Cop <- normalCopula(param=0.7, dim = 2)

## Samples ----
X_cop <- mvdc(copula=Cop, margins=c("exp", "exp"),
               paramMargins=list(list(rate=4),
                                 list(rate=2)))
  
X <- rMvdc(1000,X_cop)

Y_cop <- mvdc(copula=Cop, margins=c("norm", "norm"),
               paramMargins=list(list(mean=0, sd=2),
                                 list(mean=0, sd=1)))
Y <- rMvdc(1000,Y_cop)
par(mfrow=c(1,2))
plot(X,xlab="X1",ylab="X2")
plot(Y,xlab="Y1",ylab="Y2")
mtext("Observations", cex=2, side =3,line = -2, outer = TRUE,col='black')


#par(mfrow=c(1,2))
X.df = data.frame(X)
p1 <- ggplot(data = X.df, aes(X1, X2)) +
  geom_point()
p1 <- ggMarginal(p1, type = "densigram", fill = "grey")   

Y.df = data.frame(Y)
colnames(Y.df) = c("Y1", "Y2")
p2 <- ggplot(data = Y.df, aes(Y1, Y2)) +
  geom_point()
p2 <- ggMarginal(p2, type = "densigram", fill = "grey")   

figure <- ggarrange(p1,p2, ncol=2)
annotate_figure(figure, top = text_grob("Observations", 
                                      color = "black", 
                                      size = 16)
                )

## Dependence structure ----
UX=cbind(pexp(X[,1],rate=4), pexp(X[,2],rate=2))
UX.df = data.frame(UX)
colnames(UX.df) = c("UX1", "UX2")


UY=cbind(pnorm(Y[,1],mean=0,sd=2), pnorm(Y[,2],mean=0,sd=1))
UY.df = data.frame(UY)
colnames(UY.df) = c("UY1", "UY2")

plot(UX,xlab="UX1",ylab="UX2")
plot(UY,xlab="UY1",ylab="UY2")
mtext("Ranks", cex=2, side =3,line = -2, outer = TRUE,col='black')

r1 <- ggplot(data = UX.df, aes(UX1, UX2))+
  geom_point()
r1 <- ggMarginal(r1, type = "densigram", fill = "grey")   

r2 <- ggplot(data = UY.df, aes(UY1, UY2))+
  geom_point()
r2 <- ggMarginal(r2, type = "densigram", fill = "grey")   

figure2 <- ggarrange(r1,r2, ncol=2)
annotate_figure(figure2, top = text_grob("Ranks", 
                                        color = "black", 
                                        size = 16)
)



# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# Danish  ----
data(danish)

## Data ----
data(fire)


## Smaller dataset ----
# dataset containing only those claims where profits is different from zero
fire.red <- subset(fire, profits != 0)
attach(fire.red)

## Plots ----
par(mfrow=c(1,1))
plot(building,contents,pch=1,xlim=c(0,6),ylim=c(0,6), 
     xlab="Loss of Building",ylab="Loss of Contents")

plot(building,profits,pch=1,xlim=c(0,6),ylim=c(0,6), 
     xlab="Loss of Building",ylab="Loss of Profits")

plot(contents,profits,pch=1,xlim=c(0,6),ylim=c(0,6), 
     xlab="Loss of Contents",ylab="Loss of Profits")
#mtext("Danish Fire Insurance Claims", cex = 3,side =3,line = -3.5, outer = TRUE,col='black')

ggplot(data = fire.red, aes(building, contents))+
  geom_point()+
  xlab("Loss of Building")+
  ylab("Loss of Contents")+
  xlim(0,6)+
  ylim(0,6)+
  theme(axis.title = element_text(size = 18))

ggplot(data = fire.red, aes(building, profits))+
  geom_point()+
  xlab("Loss of Building")+
  ylab("Loss of Profits")+
  xlim(0,6)+
  ylim(0,6)+
  theme(axis.title = element_text(size = 18))

ggplot(data = fire.red, aes(contents, profits))+
  geom_point()+
  xlab("Loss of Contents")+
  ylab("Loss of Profits")+
  xlim(0,6)+
  ylim(0,6)+
  theme(axis.title = element_text(size = 18))


fire.build <- fire.red[,2]
fire.cont <- fire.red[,3]

r.build <- rank(fire.build)
r.build <- (r.build - min(r.build)) / (max(r.build) + min(r.build))

r.cont <- rank(fire.cont)
r.cont <- (r.cont - min(r.cont)) / (max(r.cont) + min(r.cont))

plot(r.build, r.cont, col = "blue",
     xlab = "loss of buildings",
     ylab = "loss of contents",
     main = "Normalized ranks of danish fire data")

ggplot(data = data.frame(r.build, r.cont),
       aes(r.build, r.cont))+
  geom_point()+
  xlab("Loss of Buildings")+
  ylab("Loss of Contents")+
  ggtitle("Normalized ranks of Danish fire data")+
  theme(axis.title = element_text(size = 16),
        plot.title = element_text(size = 18, hjust = 0.5))

x_c <- cut(r.build, 25)
y_c <- cut(r.cont, 25)

##  Calculate joint counts at cut levels:
z <- table(x_c, y_c)

zv <- c(r.build, r.cont)
zt <- matrix(c(r.build, r.cont), ncol = 2)

##  Plot as a 3D histogram:
hist3D(z=z, border="black", 
       theta=-30, phi=45, axes=TRUE,
       label=TRUE, nticks=5, 
       ticktype="detailed", space=0.1,
       main = "Histogram", xlab = "loss of building",
       ylab = "loss of content", cex.axis = 0.8,
       col = ramp.col(col = c("lightblue","blue", "darkblue")))

##  Plot as a 2D heatmap:
image2D(z=z, border="black",
        col = ramp.col(col = c("lightblue", "blue", "darkblue")))

