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

# (U1=FX1(X1),FX2(X2)) (probability integral transform)
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


## NON-PARAMETRIC ----
par(mfrow=c(1,1))
grid = seq(0, 10, by = 0.005)  # grid possible values
cov1 = rbind(c(1, 0.7),
             c(0.7, 1)) 
result = EllCopEst(dataU = UX.df, grid = grid, Sigma_m1 = solve(cov1),
                   h = 0.05, a = 1, startPoint = "identity")
plot(grid, result$g_d_norm, xlim=c(0,2), type = "l")

g1 = exp(- grid/2) / (2*pi)
g_gauss = DensityGenerator.normalize(grid, grid_g = g1, d = 2) 
lines(grid, g_gauss, col = "red")


g2 = (pi^(-1)*(gamma(3/2)) / (gamma(1/2))) * (1 + grid)^(-3/2)
g_cau = DensityGenerator.normalize(grid, grid_g = g2, d = 2)


resultY = EllCopEst(dataU = UY.df, grid = grid, Sigma_m1 = solve(cov1),
                   h = 0.05, a = 1, startPoint = "identity")
plot(grid, resultY$g_d_norm, xlim=c(0,2), type = "l")
lines(grid, g_gauss, col = "red")


colorss = c("est_g" = "red", "g_gauss" = "blue", "g_cau" = "#009900")

plot_df = cbind(grid, result$g_d_norm, resultY$g_d_norm, g_gauss, g_cau)
plot_df = data.frame(plot_df)
p1 <- ggplot(data = plot_df,aes(x = grid))+
  geom_line(aes(y = result$g_d_norm, color = "est_g"), size = 0.8) +
  geom_line(aes(y = g_gauss, color = "g_gauss"), size = 0.8) +
  geom_line(aes(y = g_cau, color = "g_cau"), size = 0.8) +
  xlim(0,3)+
  labs(x = "x", y = expression(hat(g)~(x))) +
  scale_color_manual(values = colorss) +
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )


p2 <- ggplot(data = plot_df,aes(x = grid))+
  geom_line(aes(y = resultY$g_d_norm, color = "est_g"), size = 0.8) +
  geom_line(aes(y = g_gauss, color = "g_gauss"), size = 0.8) +
  geom_line(aes(y = g_cau, color = "g_cau"), size = 0.8) +
  xlim(0,3)+
  labs(x = "y", y = expression(hat(g)~(y))) +
  scale_color_manual(values = colorss) +
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )

figure <- ggarrange(p1,p2, ncol=2)
annotate_figure(figure, top = text_grob("Density generators normal copula",
                                        size = 20))

