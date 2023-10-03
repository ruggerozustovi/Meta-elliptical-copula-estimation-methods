# CASE STUDY

rm(list=ls())
graphics.off()
cat("\f")  

# Libraries
library(quantmod)
library(ggplot2)
library(ggpubr)
library(copula)
library(ElliptCopulas)
library(VineCopula)
library(gridExtra)
library(ggExtra)
library(patchwork)
library(stats)
library(MASS)
library(plot3D)
library(goftest)
library(gofCopula)
library(tseries)
library(fitdistrplus)
library(forecast)


# source("garchAUTO.R")
Sys.setlocale("LC_TIME", "en_US.UTF-8")

# Define stock symbols and date range
symbols <- c("AAPL", "GOOGL", "AMZN", "MSFT", "META")
start_date <- "2021-01-01"
end_date <- "2023-01-01"

# Download data
getSymbols(symbols, from = start_date, to = end_date)

# APPLE VS GOOGLE ----
## Data ----
aapl_data <- AAPL
googl_data <- GOOGL


ag.data = cbind(aapl_data$AAPL.Adjusted, googl_data$GOOGL.Adjusted)
ag.data = data.frame(ag.data)
ag.data$date <- as.Date(row.names(ag.data))
colnames(ag.data) = c("aapl", "googl", "date")

p_apple1 <- ggplot(ag.data, aes(x = date, y = aapl)) +
  geom_line() +
  labs(x = "Date", y = "Apple stocks price")+
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
  ggtitle("Time series Apple stocks") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 18))
p_apple1

p_google1 <- ggplot(ag.data, aes(x = date, y = googl)) +
  geom_line() +
  labs(x = "Date", y = "Google stocks price")+
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
  ggtitle("Time series Google stocks") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 18))
p_google1

ggarrange(p_apple1, p_google1, ncol = 2)


# Log Returns ----
aapl_returns <- diff(log(aapl_data$AAPL.Adjusted))
aapl_returns = aapl_returns[-1]

googl_returns <- diff(log(googl_data$GOOGL.Adjusted))
googl_returns = googl_returns[-1]


ag.mat = cbind(as.numeric(aapl_returns),as.numeric( googl_returns))
ag.df = data.frame(ag.mat)
colnames(ag.df) = c("aapl", "googl")
ag.stock_names = c("aapl", "googl")

date.ag.mat = cbind(aapl_returns, googl_returns)
ag.df2 = data.frame(date.ag.mat)
colnames(ag.df2) = c("aapl", "googl")
ag.df2$date <- as.Date(row.names(ag.df2))


### TIME SERIES ----
p_apple2 <- ggplot(ag.df2, aes(x = date, y = aapl)) +
  geom_line() +
  labs(x = "Date", y = "Apple returns")+
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
  ggtitle("Time series daily returns Apple stocks") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 18))
p_apple2

p_google2 <- ggplot(ag.df2, aes(x = date, y = googl)) +
  geom_line() +
  labs(x = "Date", y = "Google returns")+
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
  ggtitle("Time series daily returns Google stocks") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 18))
p_google2

ggarrange(p_apple2, p_google2, ncol = 1)


#### Returns ----

# Not my function
ic_alpha= function(alpha, acf_res){
  return(qnorm((1 + (1 - alpha))/2)/sqrt(acf_res$n.used))
}

ggplot_acf_pacf= function(res_, lag, label, alpha= 0.05){
  df_= with(res_, data.frame(lag, acf))
  
  # IC alpha
  lim1= ic_alpha(alpha, res_)
  lim0= -lim1
  
  
  ggplot(data = df_, mapping = aes(x = lag, y = acf)) +
    geom_bar(stat = "identity", position = "identity") +
    labs(y= label) +
    geom_hline(aes(yintercept = lim1), linetype = 2, color = 'blue') +
    geom_hline(aes(yintercept = lim0), linetype = 2, color = 'blue')
}



retA.acf_ts= ggplot_acf_pacf(res_= acf(aapl_returns, plot= F)
                          , 20
                          , label= "ACF") +
  ggtitle("Autocorrelation plot Apple returns ") +
  theme(plot.title = element_text(hjust = 0.5))

retA2.acf_ts= ggplot_acf_pacf(res_= acf(aapl_returns^2, plot= F)
                             , 20
                             , label= "ACF") +
  ggtitle("Autocorrelation plot Apple squared-returns ") +
  theme(plot.title = element_text(hjust = 0.5))

retG.acf_ts= ggplot_acf_pacf(res_= acf(googl_returns, plot= F)
                          , 20
                          , label= "ACF") +
  ggtitle("Autocorrelation plot Google returns ") +
  theme(plot.title = element_text(hjust = 0.5))

retG2.acf_ts= ggplot_acf_pacf(res_= acf(googl_returns^2, plot= F)
                             , 20
                             , label= "ACF") +
  ggtitle("Autocorrelation plot Google squared-returns ") +
  theme(plot.title = element_text(hjust = 0.5))

ggarrange(retA2.acf_ts, retG2.acf_ts)


plot(googl_returns^2)
plot(aapl_returns^2)


adf.test(aapl_returns)
adf.test(googl_returns)

### Residuals ----
library(rugarch)
meanModel1 <- list(armaOrder = c(0,0))
varModel1 <- list(model = "sGARCH", garchOrder = c(1,1))

uspec1 <- ugarchspec(variance.model = varModel1, mean.model = meanModel1, distribution.model = "norm")
uspec2 <- ugarchspec(variance.model = varModel1, mean.model = meanModel1, distribution.model = "std")


fit.A <- ugarchfit(uspec2, data = ag.df[,1])
fit.G <- ugarchfit(uspec2, data = ag.df[,2])

fit <- c(fit.A, fit.G)


z1 = fit.A@fit$z
z2 = fit.G@fit$z
Z <- cbind(z1,z2)
U <- pobs(Z)
U = data.frame(U)
colnames(U) <- c("u1", "u2")
par(pty = "s")
plot(U, xlab = expression(hat(U)[1]), ylab = expression(hat(U)[2]))


#### indipendenza ----
A.acf_ts= ggplot_acf_pacf(res_= acf(Z[,1], plot= F)
                          , 20
                          , label= "ACF") +
    ggtitle("ACF plot Apple's arma+garch model residuals ") +
  theme(plot.title = element_text(hjust = 0.5))

G.acf_ts= ggplot_acf_pacf(res_= acf(Z[,2], plot= F)
                          , 20
                          , label= "ACF") +
  ggtitle("ACF plot Apple's arma+garch model squared-residuals ") +
  theme(plot.title = element_text(hjust = 0.5))

A2.acf_ts= ggplot_acf_pacf(res_= acf(Z[,1]^2, plot= F)
                          , 20
                          , label= "ACF") +
  ggtitle("ACF plot Google's arma+garch model residuals ") +
  theme(plot.title = element_text(hjust = 0.5))

G2.acf_ts= ggplot_acf_pacf(res_= acf(Z[,2]^2, plot= F)
                          , 20
                          , label= "ACF") +
  ggtitle("ACF plot Google's arma+garch model squared-residuals ") +
  theme(plot.title = element_text(hjust = 0.5))

ggarrange(A.acf_ts, A2.acf_ts, G.acf_ts, G2.acf_ts)

#### stazionarietà ----
Z_time = Z
Z_time = data.frame(cbind(ag.df2$date, z1, z2))
colnames(Z_time) = c("date", "z1", "z2")
Z_time$date <- as.Date(Z_time$date)

p_apple3 <- ggplot(Z_time, aes(x = date, y = z1)) +
  geom_line() +
  labs(x = "Date", y = "residuals")+
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
  ggtitle("Time series Apple model residuals") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 18))
p_apple3

p_google3 <- ggplot(Z_time, aes(x = date, y = z2)) +
  geom_line() +
  labs(x = "Date", y = "residuals")+
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
  ggtitle("Time series Google model residuals") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 18))
p_google3

ggarrange(p_apple3, p_google3, ncol = 1)

adf.test(z1)
adf.test(z2)

## Univariate Analysis ----

#### Histograms ----
gp1 <- ggplot(data.frame(x = z1), aes(x = x)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "grey") +
  stat_function(fun = dnorm, args = list(mean = mean(z1), sd = sd(z1)))+
  xlab("Returns") +
  labs(title = "Fitted Distribution for Apple returns")+
  theme(plot.title = element_text(hjust = 0.5))

gp2 <- ggplot(data.frame(x =z2), aes(x = x)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "grey") +
  stat_function(fun = dnorm, args = list(mean = mean(z2), sd = sd(z2)))+
  xlab("Returns") +
  labs(title = "Fitted Distribution for Google returns")+
  theme(plot.title = element_text(hjust = 0.5))


ggarrange(gp1,gp2)

#### MARGINALs ESTIMATION ----
library(metRology)
library(qqplotr)

#### apple residuals ----
A.par.est = fitdistr(z1, "normal")
fit.t.a <- fitdistrplus::fitdist(z1, "t.scaled", 
                               start = list(mean = mean(z1),
                                            sd = sd(z1),
                                            df = 30))

ggplot(mapping = aes(sample = z1))+
  stat_qq_band(distribution = "t.scaled", 
               dparams = fit.t.a$estimate)+
  stat_qq_point(distribution = "t.scaled", 
                dparams = fit.t.a$estimate,
                size = 2)+
  stat_qq_line(distribution = "t.scaled", 
               dparams = fit.t.a$estimate)+
  xlab("Theoretical quantile")+
  ylab("Sample quantile")+
  ggtitle("Student-t Q-Q plot of residuals")+
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title = element_text(size = 16))


##### google residuals ----
G.par.est = fitdistr(z2, "t")
G.par.est.n = fitdistr(z2, "normal")

fit.n.g <- fitdistrplus::fitdist(z2, "norm", 
                               start = list(mean = mean(z2),
                                            sd = sd(z2)))

fit.t <- fitdistrplus::fitdist(z2, "t.scaled", 
                      start = list(mean = mean(z2),
                                   sd = sd(z2),
                                   df = 30))

ggplot(mapping = aes(sample = z2))+
  stat_qq_band(distribution = "t.scaled", 
               dparams = fit.t$estimate)+
  stat_qq_point(distribution = "t.scaled", 
                dparams = fit.t$estimate,
                size = 2)+
  stat_qq_line(distribution = "t.scaled", 
               dparams = fit.t$estimate)+
  xlab("Theoretical quantile")+
  ylab("Sample quantile")+
  ggtitle("Student-t Q-Q plot of residuals")+
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title = element_text(size = 16))


#### Tests ----
shapiro.test(z1)
ks.test(z1, "pnorm")
qqnorm(z1)
qqline(z1, col = "red", lwd = 2)


shapiro.test(z2)
ks.test(z2, "pnorm")
qqnorm(z2)
qqline(z2, col = "red", lwd = 2)


## Dependence ----
p = ggplot(data = data.frame(Z), aes(z1, z2)) +
  geom_point() +
  ggtitle("Observations") +
  xlab("Apple residuals") +
  ylab("Google residuals") +
  theme( plot.title = element_text(hjust = 0.5, size = 14))
p = ggMarginal(p, type = "densigram", fill = "grey")


p2 = ggplot(data = U, aes(u1, u2))+
  geom_point() +
  ggtitle("Pseudo-observations") +
  theme( plot.title = element_text(hjust = 0.5, size = 14))
ggMarginal(p2, type = "histogram", fill = "grey")

fig1 = ggarrange(p, p2, ncol = 2)
fig1 + plot_annotation(
  title = "Dependence between residuals",
  theme = theme(
    plot.title = element_text(hjust = 0.5, size = 18)
  )
)

### Kendall's tau ----

ag.tau = KTMatrixEst(U)
ag.tau
stats::cor.test(U[,1], U[,2],
                alternative = "two.sided",
                method = "kendall",
                exact = NULL, conf.level = 0.95, continuity = FALSE)

sigma = sin((pi*ag.tau[1,2])/2)

cov1 = rbind(c(1, sigma),
             c(sigma, 1)) 


## PARAMETRIC ESTIMATION ----
# selectedCopula <- BiCopSelect(U[,1], U[,2],
#                               familyset = NA)
# selectedCopula

# selected.n.Copula <- BiCopEst(U[,1], U[,2],
#                               family = 1)
# selected.n.Copula
# selected.t.Copula <- BiCopEst(U[,1], U[,2],
#                               family = 2)
# selected.t.Copula

fitT <- fitCopula(ellipCopula("t", dim = 2), data = U, method = "ml")
summary(fitT)
fitT@estimate
dof = fitT@estimate[2]

fitnorm <- fitCopula(ellipCopula("norm", dim = 2), data = U, method = "ml")
summary(fitnorm)
fitnorm@estimate

sim_norm = rCopula(nrow(ag.df), copula = fitnorm@copula)
sim_t = rCopula(nrow(ag.df), copula = fitT@copula)

grid = seq(0, 10, by = 0.01)
sim2_norm = EllCopSim(n = 1000, grid=grid,d = 2, 
                      g_d = exp(- grid/2) / (2*pi),
                      A = chol(cov1))

t5.stud = (((pi * dof)^(-1)*(gamma((2+dof)/2)) / gamma(dof/2))) * (1 + (grid/dof))^(-((2+dof)/2))
g_studt = DensityGenerator.normalize(grid, grid_g = t5.stud, d = 2)


sim2_t5 = EllCopSim(n = 1000, grid=grid,d = 2, 
                    g_d = t5.stud,
                    A = chol(cov1))



n1 = ggplot(data = U, aes(u1, u2))+
  geom_point()+
  geom_point(data = as.data.frame(sim2_norm), 
             aes(V1, V2),
             color = "red")+
  ggtitle("pseudo-observations and normal copula")+
  theme(plot.title = element_text(hjust = 0.5))

t1 = ggplot(data = U, aes(u1, u2))+
  geom_point()+
  geom_point(data = as.data.frame(sim2_t5),
             aes(V1, V2),
             color = "blue")+
  ggtitle("pseudo-observations and t-copula")+
  theme(plot.title = element_text(hjust = 0.5))

ggarrange(n1, t1, ncol = 2)

#### 1) Goodness-of-Fit tests ----
set.seed(1234)
gof2_norm = gofKS("normal", Z, M = 50)
gof2_norm

gof2_stud5 = gofKS("t",  Z, df = dof, M = 50, df.est = F)
gof2_stud5


# NON-PARAMETRIC ----
# grid = seq(0, 10, by = 0.01)


te_dis = TEllDistrEst(Z, h = 0.05, grid = grid)

pseudo_U = pobs(Z)
h = bw.ucv(pseudo_U)
cop_dis = EllCopEst(U, 
                    grid = grid, 
                    Sigma_m1 = solve(cov1), 
                    h = h, niter = 10,
                    startPoint = "gaussian")
cop_est = cop_dis$g_d_norm

g1 = exp(- grid/2) / (2*pi)
g_norm = DensityGenerator.normalize(grid, grid_g = g1, d = 2)

t6.stud = (pi * 6)^(-(3/2))*(gamma(9/2) / gamma(3)) * (1 + (grid/6))^(-(9/2))
g_studt6 = DensityGenerator.normalize(grid, grid_g = t6.stud, d = 2)


est = te_dis$estEllCopGen
g_df = cbind(grid, cop_est, g_norm, g_studt, g_studt6)
g_df = data.frame(g_df)

colorss = c("est" = "black",
            "g_norm" = "red",
            "g_studt" = "blue"
)

p3 = ggplot(data = g_df, aes(x=grid))+
  geom_line(aes(y = cop_est, color = "est")) +
  geom_line(aes(y = g_norm, color = "g_norm")) +
  geom_line(aes(y = g_studt, color = "g_studt")) +
  xlim(0,2)+
  ggtitle("Estimation of the generator function") +
  theme(plot.title = element_text(hjust = 0.5, size = 18)) +
  labs(x = "x", y = expression(hat(g)~(x)), color = "generators") +
  scale_color_manual(values = colorss)

p3

g.na = cop_est[1:200]

sim_U = EllCopSim(n = 1000, d = 2,  grid[1:length(g.na)],
              A = chol(cov1),
              g_d = g.na)

sim_X = matrix(nrow = 1000, ncol = 2)
sim_X[,1] = (metRology::qt.scaled(sim_U[,1],
                                  df = fit.t.a$estimate[3],
                                  mean = fit.t.a$estimate[1],
                                  sd = fit.t.a$estimate[2]))
sim_X[,2] = (metRology::qt.scaled(sim_U[,2],
                                  df = fit.t$estimate[3],
                                  mean = fit.t$estimate[1],
                                  sd = fit.t$estimate[2]))


gdenT = ggplot() +
  geom_density_2d(data = data.frame(U), 
                  aes(u1, u2, colour = after_stat(level)))+
  labs(x = "u1",
       y = "u2") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.95, 0.05),
        legend.justification = c("bottom"),
        legend.box.just = "bottom",
        legend.margin = margin(6, 6, 6, 6))

gscatT = ggplot()+
  geom_point(data = data.frame(U), aes(u1, u2))+
  labs(x = "u1",
       y = "u2") +
  theme(plot.title = element_text(hjust = 0.5))


gdenE = ggplot() +
  geom_density_2d(data = data.frame(sim_U), 
                  aes(X1, X2, colour = after_stat(level)))+
  labs(x = "u1",
       y = "u2") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.95, 0.05),
        legend.justification = c("bottom"),
        legend.box.just = "bottom",
        legend.margin = margin(6, 6, 6, 6))

gscatE = ggplot()+
  geom_point(data = data.frame(sim_U), aes(X1, X2))+
  labs(x = "u1",
       y = "u2") +
  theme(plot.title = element_text(hjust = 0.5))



sim1 = rCopula(1000,fitT@copula)
colnames(sim1) = c("u1", "u2")

gden1 = ggplot() +
  geom_density_2d(data = data.frame(sim1), 
                  aes(u1, u2, colour = after_stat(level)))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.95, 0.05),
        legend.justification = c("bottom"),
        legend.box.just = "bottom",
        legend.margin = margin(6, 6, 6, 6))

gscat1 = ggplot() +
  geom_point(data = data.frame(sim1), aes(u1, u2))+
  theme(plot.title = element_text(hjust = 0.5))


sim2 = rCopula(1000,fitnorm@copula)
colnames(sim2) = c("u1", "u2")
gden2 = ggplot() +
  geom_density_2d(data = data.frame(sim2), 
                  aes(u1, u2, colour = after_stat(level)))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.95, 0.05),
        legend.justification = c("bottom"),
        legend.box.just = "bottom",
        legend.margin = margin(6, 6, 6, 6))


gscat2 = ggplot() +
  geom_point(data = data.frame(sim2), aes(u1, u2))+
  theme(plot.title = element_text(hjust = 0.5))

fig4 = ggarrange(gdenT, gscatT,
          gdenE, gscatE, 
          ncol = 2, nrow = 2,
          common.legend = TRUE,
          legend = "bottom")
annotate_figure(fig4,
                left = text_grob( "Non-parametric simulation               Observations",
                                  rot = 90,
                                  size = 16))

fig5 = ggarrange(gden1, gscat1,
          gden2, gscat2,
          ncol = 2, nrow = 2,
          common.legend = TRUE,
          legend = "bottom")

annotate_figure(fig5,
                left = text_grob( "Gaussian                            Student-t",
                                  rot = 90,
                                  size = 16))



















