# ///////////////////////////////// #
#         Simulation 1 multi        #  ---- 
#               Unibo               #
#     Academic year 2022-2023       #
#         Ruggero Zustovi           #
# ///////////////////////////////// #

rm(list=ls())
graphics.off()
cat("\f")   # it delete the console

# Libraries ----
library(ElliptCopulas)
library(mvtnorm)
library(copula)
library(stats)
library(gridExtra)
library(ggplot2)
library(ggpubr)

## Initialization ----
set.seed(1234)
d = 2 # dimension
n = 1000 # sample size
grid = seq(0, 10, by = 0.005)  # grid possible values
a = 1 # parameter a in the Liebscher's procedure
h = 0.05 # bandwidth for KDE

gplot_list = list()
te.gplot_list = list()

#### correlation matrix ----
cov1 = rbind(c(1, 0.3),
             c(0.3, 1))

#### Generators ----
g_list = list(g2 = sqrt(grid)/(1+grid^3),
              g4 = exp(-grid)*sin(2*grid)^2)


g_list.names = c("g2(x) = sqrt(x)/(1+x^3)",
                 "g4(x) = exp(-x)*sin(2*x)^2")
U = list()
X = list()
#g_norm = list()
te_distr = list()
df = list()
MISE_list = list()
MISE_h = list()
points = list()

h1 = c()
h2 = c()
h3 = c()
h4 = c()
bh = list()


## Data generation ----
for(i in 1:length(g_list)){
  set.seed(12345)
  ## Data generation ----
  g_norm = DensityGenerator.normalize(grid, grid_g = g_list[[i]], d = d) 
  
  U[[i]] = EllCopSim(n = n, d = d, grid = grid, g_d = g_norm, A = chol(cov1))
  
  X[[i]] = matrix(nrow = n, ncol = d)
  X[[i]][,1] = stats::qunif(U[[i]][,1], min = 0, max = 1)
  X[[i]][,2] = stats::qunif(U[[i]][,2], min = 0, max = 1)
  
  ## Estimation g meta-elliptical ----
  est = list()
  MISE = c()
  
  b = c(U[[i]][,1],U[[i]][,2])
  h1[i] = bw.nrd(b)
  h2[i] = bw.nrd0(b)
  h3[i] = bw.ucv(b)
  h4[i] = bw.bcv(b)
  
  bh[[i]] = c(h1[i], h2[i], h3[i], h4[i])
  
  color_vector <- c("blue", "green", "orange", "red")
  bandwidth.names = c("Scott","Silverman", "unbiased cv", "biased cv")
  colorss = c("Scott" = "blue",
              "Silverman" = "green",
              "unbiased cv" = "orange",
              "biased cv" = "red")
  
  for (k in 1:4){
    result = EllCopEst(dataU = U[[i]], grid = grid, Sigma_m1 = solve(cov1),
                       h = bh[[i]][k], a = a, 
                       startPoint = "identity", Kernel = "gaussian", niter = 10)
    est[[k]] <- result
    
    ### MISE ----
    MISE[k] <- round(mean((est[[k]]$g_d_norm - g_norm)^2), digits = 6)
  }
  
  MISE_h[[i]] = MISE
  
  n_lines = length(est)
  df[[i]] <- data.frame(grid, g_norm)
  p = ggplot(df[[i]], aes(grid, g_norm)) + geom_line(linewidth = 1.3) + xlim(0, 1.5)
  line_data = list()
  
  for (j in 1:n_lines) {
    line_data[[j]] <- data.frame(grid, g_d_norm = est[[j]]$g_d_norm)
  }
  
  p = p + 
    geom_line(data = line_data[[1]], aes(grid, g_d_norm, color = "Scott"), linewidth = 0.8) +
    geom_line(data = line_data[[2]], aes(grid, g_d_norm, color = "Silverman"), linewidth = 0.8) +
    geom_line(data = line_data[[3]], aes(grid, g_d_norm, color = "unbiased cv"), linewidth = 0.8) +
    geom_line(data = line_data[[4]], aes(grid, g_d_norm, color = "biased cv"), linewidth = 0.8) +
    ggtitle(g_list.names[i]) +
    scale_color_manual(values = colorss, name='Bandwidth') +
    labs(x = "x", y = expression(hat(g)~(x)), color = "Bandwidth") +
    theme(
      legend.position = c(.95, .95),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(6, 6, 6, 6),
      plot.title = element_text(hjust = 0.5, size = 20)
    )
  
  print(p)
  gplot_list[[i]] <- p
}

fig1 <- ggarrange(gplot_list[[1]],gplot_list[[2]])
fig1

for(j in 1:length(g_list)){
  diff_h = cbind(bh[[1]], bh[[2]])
  diff_mise = cbind(MISE_h[[1]], MISE_h[[2]])
  points[[j]] = data.frame(diff_h[,j], diff_mise[,j])
  colnames(points[[j]]) = c("diff_h", "diff_mise")
}



for(i in 1:length(g_list)){
  set.seed(12345)
  ## Data generation ----
  g_norm = DensityGenerator.normalize(grid, grid_g = g_list[[i]], d = d) 
  
  U[[i]] = EllCopSim(n = n, d = d, grid = grid, g_d = g_norm, A = chol(cov1))
  
  X[[i]] = matrix(nrow = n, ncol = d)
  X[[i]][,1] = stats::qunif(U[[i]][,1], min = 0, max = 1)
  X[[i]][,2] = stats::qunif(U[[i]][,2], min = 0, max = 1)
  
  ## Estimation g meta-elliptical ----
  est = list()
  MISE = c()
  h_range  = seq(0, 0.1, by = 0.001)
  h = sort(sample(h_range, size = 50, replace = TRUE))
  color_vector <- colors()[1:50]
  
  for (k in 1:length(h)){
    result = EllCopEst(dataU = U[[i]], grid = grid, Sigma_m1 = solve(cov1),
                       h = h[k], a = a, 
                       startPoint = "identity", Kernel = "gaussian",niter = 10)
    est[[k]] <- result
    
    ### MISE ----
    MISE[k] <- round(mean((est[[k]]$g_d_norm - g_norm)^2), digits = 6)
  }
  
  MISE_list[[i]] = MISE
  
  n_lines = length(est)
  df[[i]] <- data.frame(grid, g_norm)
  p = ggplot(df[[i]], aes(grid, g_norm), size = 1) + geom_line() + xlim(0, 2)
  
  for (j in 1:n_lines) {
    line_data <- data.frame(grid, g_d_norm = est[[j]]$g_d_norm)
    p <- p + geom_line(data = line_data, aes(grid, g_d_norm), 
                       color = color_vector[j]) +
      ggtitle(g_list.names[i]) +
      labs(x = "x", y = expression(hat(g)~(x))) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  gplot_list[[i]] <- p
  
}

grid_arrange <- do.call(grid.arrange, c(gplot_list, ncol = 2))

selected_bandwidth = c("Scott", "Silverman", "UCV", "BCV")
points.df1 <- data.frame(selected_bandwidth, points[[1]])
points.df2 <- data.frame(selected_bandwidth, points[[2]])


error.df1 = cbind(h, MISE_list[[1]])
error.df1 = data.frame(error.df1)
plot_h1 <- ggplot() +
  geom_line(data = error.df1, aes(h, MISE_list[[1]])) +
  geom_point(data = points.df1, aes(diff_h, diff_mise), color = "red") +
  geom_text(data = points.df1, aes(diff_h, diff_mise),
            label = selected_bandwidth,
            size = 3.5,
            nudge_y = 0.00005,
            # nudge_x = -0.0005,
            color = "black") +
  ggtitle(g_list.names[1]) +
  labs(x = "h", y = "MISE") +
  theme(plot.title = element_text(hjust = 0.5))
plot_h1

error.df2 = cbind(h, MISE_list[[2]])
error.df2 = data.frame(error.df2)
plot_h2 <- ggplot() +
  geom_line(data = error.df2, aes(h, MISE_list[[2]])) +
  geom_point(data = points.df2, aes(diff_h, diff_mise), color = "red") +
  geom_text(data = points.df2, aes(diff_h, diff_mise),
            label = selected_bandwidth,
            size = 3.5,
            nudge_y = 0.0005,
            nudge_x = 0.0008,
            color = "black") +          
  ggtitle(g_list.names[2]) +
  labs(x = "h", y = "MISE") +
  theme(plot.title = element_text(hjust = 0.5))


figure_h <- ggarrange(plot_h1, plot_h2, ncol = 2)
figure_h







