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


df <- fire.red[,2:4]
pseudo <- pobs(df)

plot(pseudo[,1], pseudo[,2], col = "blue",
     xlab = "loss of buildings",
     ylab = "loss of contents",
     main = "Normalized ranks of danish fire data"
)

x_c <- cut(pseudo[,1], 25)
y_c <- cut(pseudo[,2], 25)

##  Calculate joint counts at cut levels:
z <- table(x_c, y_c)

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



## NON PARAMETRIC ----

grid = seq(0, 10, by = 0.05)  # grid possible values

te_distr = TEllDistrEst(df, h = 0.05, grid = grid)

tau = KTMatrixEst(pseudo)
tau.x = KTMatrixEst(df)
# gen_est = EllCopEst(pseudo, grid = grid, Sigma_m1 = solve(tau), d = 2)

plot(grid, te_distr$estEllCopGen,  type = "l", xlim =c(0,2))
g1 = exp(- grid/2) / (2*pi)
g_norm = DensityGenerator.normalize(grid, grid_g = g1, d = 2) 
lines(grid, g_norm, col = "red")

t6.stud = (pi * 6)^(-(3/2))*(gamma(9/2) / gamma(3)) * (1 + (grid/6))^(-(9/2))
g_studt6 = DensityGenerator.normalize(grid, grid_g = t6.stud, d = 2)

t9.stud = (pi * 9)^(-(3/2))*(gamma(6) / gamma(9/2)) * (1 + (grid/9))^(-6)
g_studt9 = DensityGenerator.normalize(grid, grid_g = t9.stud, d = 2)



result = EllCopEst(dataU = pseudo[,1:2], Sigma_m1 = solve(tau[1:2, 1:2]),
                   h = 0.005, grid = grid)

plot(pseudo[,1], pseudo[,2])

est = result$g_d_norm
est_g = te_distr$estEllCopGen

df2 = cbind(grid, est, g_norm, g_studt6, g_studt9)
df2 = data.frame(df2)

colorss = c("est" = "black", 
            "g_norm" = "green", 
            "g_studt6" = "blue", 
            "g_studt9" = "red"
            )

ggplot(data = df2,aes(x = grid))+
  geom_line(aes(y = est_g, color = "est"), size = 1) +
  geom_line(aes(y = g_norm, color = "g_norm"), size = 0.8) +
  geom_line(aes(y = g_studt6, color = "g_studt6"), size = 0.8) +
  geom_line(aes(y = g_studt9, color = "g_studt9"), size = 0.8) +
  xlim(0,3)+
  labs(x = "x", y = expression(hat(g)~(x)), color = "generators") +
  scale_color_manual(values = colorss) +
  ggtitle("Danish fire data meta-elliptical copula generator function") +
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    plot.title = element_text(hjust = 0.5, size = 20)
  )



cop2 = EllCopSim(1000, d = 2, grid = grid, g_d = est, A = chol(tau[1:2, 1:2]))

df.cop2 = data.frame(cop2)
ggplot(data = df.cop2, aes(X1, X2))+
  geom_point()








