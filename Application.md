# Mixed Latent Position Cluster Model Simulation Studies and Real Data Application

This tutorial provides all the R code to fit the Mixed Latent Position Cluster Model to network data.
The detailed code for the simulation studies and real data application of the [Mixed Latent Position Cluster Model] paper is also provided.

We start from loading all the functions and packages we need:

``` r
rm(list=ls())
gc()
source("Functions.R")

library("igraph")
library("RColorBrewer")
library("pheatmap")
library("ggplot2")
library("latex2exp")
library("IMIFA")
library("grid")
library("gridExtra")

# Define colors
# display.brewer.all(type="seq")
# brewer.pal.info[brewer.pal.info$category == "seq",]
# display.brewer.pal(9,"RdBu")
My_colors <- c(brewer.pal(10,"RdBu")[c(4,7)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[3],
               brewer.pal(10,"RdBu")[c(2,9)],brewer.pal(10,"PRGn")[c(9,2)],brewer.pal(9,"YlOrBr")[6],
               brewer.pal(9,"RdPu")[5],brewer.pal(9,"GnBu")[5],brewer.pal(9,"Reds")[c(6,9)],brewer.pal(9,"Greys")[c(3,6,9)])
library("hues")
swatch(My_colors)
```

## 1. Two Motivating MLPCM Network Examples

The 1st 3-node motivating network example generated from the MLPCM shown in Section 2 of the paper can be reproduced by the following codes:

``` r
set.seed(1)
# Plot the motivation example of the new LPCM with 5 nodes
MixedLPCM_MotivationPlot_GraphMatrix_3NodesExample <- matrix(0,9,9)
MixedLPCM_MotivationPlot_GraphMatrix_3NodesExample[1,c(6,8)] <- rep(1,2)
MixedLPCM_MotivationPlot_GraphMatrix_3NodesExample[2,c(4,9)] <- rep(1,2)
MixedLPCM_MotivationPlot_GraphMatrix_3NodesExample[3,c(5,7)] <- rep(1,2)
library(mvtnorm)
MixedLPCM_MotivationPlot_positions_3NodesExample <- rbind(c(0,-1),c(1*sinpi(2/3),-1*cospi(2/3)),c(1*sinpi(4/3),-1*cospi(4/3)),
                                                          mvtnorm::rmvnorm(2,c(0,-1),(1/10)*diag(2)), # generate overt positions with precision 7 for gorup 1
                                                          mvtnorm::rmvnorm(2,c(1*sinpi(2/3),-1*cospi(2/3)),(1/8)*diag(2)),# with precision 4 for gorup 2
                                                          mvtnorm::rmvnorm(2,c(1*sinpi(4/3),-1*cospi(4/3)),(1/6)*diag(2)))# with precision 1 for gorup 3

library("igraph")
g_obs_3NodesExample <- graph_from_adjacency_matrix(MixedLPCM_MotivationPlot_GraphMatrix_3NodesExample,mode = "directed")
# Node sizes are proportional to their social impact, \tilde{\gamma}_i
V(g_obs_3NodesExample)$size <- c(rep(15,3),rep(10,6))
# Additional graphical settings
V(g_obs_3NodesExample)$frame.color <- "black"
library("latex2exp")
V(g_obs_3NodesExample)$label <- c(TeX(r"($u_1$)"),TeX(r"($u_2$)"),TeX(r"($u_3$)"),
                                  TeX(r"($v_{1 \leftarrow 2}$)"),TeX(r"($v_{1 \leftarrow 3}$)"),
                                  TeX(r"($v_{2 \leftarrow 1}$)"),TeX(r"($v_{2 \leftarrow 3}$)"),
                                  TeX(r"($v_{3 \leftarrow 1}$)"),TeX(r"($v_{3 \leftarrow 2}$)"))
V(g_obs_3NodesExample)$label.cex <- c(rep(1,3),rep(0.65,6))


# Node colors indicate the inferred clustering
Customized_colors <- My_colors[1:10]
V(g_obs_3NodesExample)$color <- adjustcolor(c(Customized_colors[6:8],rep(Customized_colors[1],2),rep(Customized_colors[2],2),rep(Customized_colors[3],2)), alpha.f = .7)

# Assign interactions (0 or 1) to the network
E(g_obs_3NodesExample)$weight <- c(0); beta_ex <- 1.7
exp_beta_ex_minus_dist2_list <- c()
for (i in 1:nrow(get.edgelist(g_obs_3NodesExample))){
  exp_beta_ex_minus_dist2 <- exp(beta_ex-dist(MixedLPCM_MotivationPlot_positions_3NodesExample[get.edgelist(g_obs_3NodesExample)[i,],])^2)
  exp_beta_ex_minus_dist2_list <- c(exp_beta_ex_minus_dist2_list,exp_beta_ex_minus_dist2)
  E(g_obs_3NodesExample)$weight[i] <- round(exp_beta_ex_minus_dist2) # specify the edge weights by round(PoisMean)
}
exp_beta_ex_minus_dist2_list
E(g_obs_3NodesExample)$weight
E(g_obs_3NodesExample)$color <- colorRampPalette(brewer.pal(9,"Greys")[c(3,6)])(max(E(g_obs_3NodesExample)$weight)+1)[E(g_obs_3NodesExample)$weight+1]

# Make the position plot
par(mfrow=c(1,1),mai = c(0.05, 0.05, 0.05, 0.05),mgp = c(1,0.25,0))
plot(g_obs_3NodesExample, rescale=T,layout=MixedLPCM_MotivationPlot_positions_3NodesExample,edge.curved=0.0,edge.arrow.size=0.75,
     edge.lty=2-(E(g_obs_3NodesExample)$weight>0),edge.width=(E(g_obs_3NodesExample)$weight+1)*1.5)
par(mfrow=c(1,1),mai = c(1.02, 0.82, 0.82, 0.42))
set.seed(NULL)
```

The more complex 5-node one can be reproduced via:

``` r
set.seed(16)
# Plot the motivation example of the new LPCM with 5 nodes
MixedLPCM_MotivationPlot_GraphMatrix_5NodesExample <- matrix(0,25,25)
MixedLPCM_MotivationPlot_GraphMatrix_5NodesExample[1,c(10,14,18,22)] <- rep(1,4)
MixedLPCM_MotivationPlot_GraphMatrix_5NodesExample[2,c(6,15,19,23)] <- rep(1,4)
MixedLPCM_MotivationPlot_GraphMatrix_5NodesExample[3,c(7,11,20,24)] <- rep(1,4)
MixedLPCM_MotivationPlot_GraphMatrix_5NodesExample[4,c(8,12,16,25)] <- rep(1,4)
MixedLPCM_MotivationPlot_GraphMatrix_5NodesExample[5,c(9,13,17,21)] <- rep(1,4)
library(mvtnorm)
MixedLPCM_MotivationPlot_positions_5NodesExample <- rbind(c(0,1),c(-1*sinpi(0.4),1*cospi(0.4)),c(-1*sinpi(0.8),1*cospi(0.8)),c(-1*sinpi(1.2),1*cospi(1.2)),c(-1*sinpi(1.6),1*cospi(1.6)),
                                                          mvtnorm::rmvnorm(4,c(0,1),(1/9)*diag(2)), # generate overt positions with precision 9 for gorup 1
                                                          mvtnorm::rmvnorm(4,c(-1*sinpi(0.4),1*cospi(0.4)),(1/15)*diag(2)),# with precision 7 for gorup 2
                                                          mvtnorm::rmvnorm(4,c(-1*sinpi(0.8),1*cospi(0.8)),(1/15)*diag(2)),# with precision 5 for gorup 3
                                                          mvtnorm::rmvnorm(4,c(-1*sinpi(1.2),1*cospi(1.2)),(1/15)*diag(2)),# with precision 3 for gorup 4
                                                          mvtnorm::rmvnorm(4,c(-1*sinpi(1.6),1*cospi(1.6)),(1/15)*diag(2)))# with precision 1 for gorup 5
library("igraph")
g_obs_5NodesExample <- graph_from_adjacency_matrix(MixedLPCM_MotivationPlot_GraphMatrix_5NodesExample,mode = "directed")
# Node sizes are proportional to their social impact, \tilde{\gamma}_i
V(g_obs_5NodesExample)$size <- c(rep(15,5),rep(10,20))
# Additional graphical settings
V(g_obs_5NodesExample)$frame.color <- "black"
# V(g_obs_5NodesExample)$label <- c(1:5,rep("",20))
V(g_obs_5NodesExample)$label <- c(TeX(r"($u_1$)"),TeX(r"($u_2$)"),TeX(r"($u_3$)"),TeX(r"($u_4$)"),TeX(r"($u_5$)"),
                                  TeX(r"($v_{1 \leftarrow 2}$)"),TeX(r"($v_{1 \leftarrow 3}$)"),TeX(r"($v_{1 \leftarrow 4}$)"),TeX(r"($v_{1 \leftarrow 5}$)"),
                                  TeX(r"($v_{2 \leftarrow 1}$)"),TeX(r"($v_{2 \leftarrow 3}$)"),TeX(r"($v_{2 \leftarrow 4}$)"),TeX(r"($v_{2 \leftarrow 5}$)"),
                                  TeX(r"($v_{3 \leftarrow 1}$)"),TeX(r"($v_{3 \leftarrow 2}$)"),TeX(r"($v_{3 \leftarrow 4}$)"),TeX(r"($v_{3 \leftarrow 5}$)"),
                                  TeX(r"($v_{4 \leftarrow 1}$)"),TeX(r"($v_{4 \leftarrow 2}$)"),TeX(r"($v_{4 \leftarrow 3}$)"),TeX(r"($v_{4 \leftarrow 5}$)"),
                                  TeX(r"($v_{5 \leftarrow 1}$)"),TeX(r"($v_{5 \leftarrow 2}$)"),TeX(r"($v_{5 \leftarrow 3}$)"),TeX(r"($v_{5 \leftarrow 4}$)"))
V(g_obs_5NodesExample)$label.cex <- c(rep(1,5),rep(0.65,20))


# Node colors indicate the inferred clustering
Customized_colors <- My_colors[1:10]
V(g_obs_5NodesExample)$color <- adjustcolor(c(Customized_colors[6:10],rep(Customized_colors[1],4),rep(Customized_colors[2],4),rep(Customized_colors[3],4),rep(Customized_colors[4],4),rep(Customized_colors[5],4)), alpha.f = .7)

# Assign interactions (0 or 1) to the network
E(g_obs_5NodesExample)$weight <- c(0); beta_ex <- 1.5
exp_beta_ex_minus_dist2_list <- c()
for (i in 1:nrow(get.edgelist(g_obs_5NodesExample))){
  exp_beta_ex_minus_dist2 <- exp(beta_ex-dist(MixedLPCM_MotivationPlot_positions_5NodesExample[get.edgelist(g_obs_5NodesExample)[i,],])^2)
  exp_beta_ex_minus_dist2_list <- c(exp_beta_ex_minus_dist2_list,exp_beta_ex_minus_dist2)
  E(g_obs_5NodesExample)$weight[i] <- round(exp_beta_ex_minus_dist2) # # specify the edge weights by round(PoisMean)
}
exp_beta_ex_minus_dist2_list
E(g_obs_5NodesExample)$weight
E(g_obs_5NodesExample)$color <- colorRampPalette(brewer.pal(9,"Greys")[c(3,6)])(max(E(g_obs_5NodesExample)$weight)+1)[E(g_obs_5NodesExample)$weight+1]

# Make the position plot
par(mfrow=c(1,1),mai = c(0.05, 0.05, 0.05, 0.05),mgp = c(1,0.25,0))
plot(g_obs_5NodesExample, rescale=T,layout=MixedLPCM_MotivationPlot_positions_5NodesExample,edge.curved=0.0,edge.arrow.size=0.5,
     edge.lty=2-(E(g_obs_5NodesExample)$weight>0),edge.width=(E(g_obs_5NodesExample)$weight+1)*1.25)
par(mfrow=c(1,1),mai = c(1.02, 0.82, 0.82, 0.42))
set.seed(NULL)
```

## 2. Simulation Study 1

### 2.1 Simulation Study 1 Scenario 1

We begin by simulating from a Mixed-LPCM with model settings provided in the simulation study 1 section of the paper.

``` r
# Simulation study 1 scenario 1

# Simulate from the MLPCM and store the network

# SS1Sce1_MixedLPCM <- Simulation_Directed_MixedLPCM(beta=1,mu=rbind(c(1.25,1.25),c(1.25,-1.25),c(-1.25,-1.25),c(-1.25,1.25)),tau=c(8,6,4,2),
#                                                    gamma=rep(10,100),d=2,
#                                                    z=c(rep(1,25),rep(2,25),rep(3,25),rep(4,25)),seed=NULL)
# write.csv(SS1Sce1_MixedLPCM$Y,"SS1Sce1_MixedLPCM_Y.csv", row.names = FALSE)
# write.csv(SS1Sce1_MixedLPCM$z,"SS1Sce1_MixedLPCM_z.csv", row.names = FALSE)
# write.csv(SS1Sce1_MixedLPCM$U,"SS1Sce1_MixedLPCM_U.csv", row.names = FALSE)
# write.csv(SS1Sce1_MixedLPCM$V,"SS1Sce1_MixedLPCM_V.csv", row.names = FALSE)

# Load the data

SS1Sce1_MixedLPCM <- list(Y = as.matrix(read.csv("SS1Sce1_MixedLPCM_Y.csv",header = TRUE)),
                              z = c(as.matrix(read.csv("SS1Sce1_MixedLPCM_z.csv",header = TRUE))),
                              U = as.matrix(read.csv("SS1Sce1_MixedLPCM_U.csv",header = TRUE)))
SS1Sce1_MixedLPCM$Z <- t(t(matrix(SS1Sce1_MixedLPCM$z,length(SS1Sce1_MixedLPCM$z),max(SS1Sce1_MixedLPCM$z)))==(1:max(SS1Sce1_MixedLPCM$z)))*1
SS1Sce1_MixedLPCM$V <- array(as.matrix(read.csv("SS1Sce1_MixedLPCM_V.csv",header = TRUE)),
                                 dim=c(dim(t(SS1Sce1_MixedLPCM$U)),nrow(SS1Sce1_MixedLPCM$Y)))
colnames(SS1Sce1_MixedLPCM$Y) <- colnames(SS1Sce1_MixedLPCM$U) <- colnames(SS1Sce1_MixedLPCM$V) <- NULL

# Reference values
SS1Sce1_MixedLPCM_beta <- 1
SS1Sce1_MixedLPCM_mu <- rbind(c(1.25,1.25),c(1.25,-1.25),c(-1.25,-1.25),c(-1.25,1.25))
SS1Sce1_MixedLPCM_tau <- c(8,6,4,2)
SS1Sce1_MixedLPCM_gamma <- rep(10,100)
```

#### 2.1.1 Fitting the MLPCM in Simulation Study 1 Scenario 1

Following the code below provides the model fitting results illustrated in the paper.

``` r
# Apply VB on the simulated MixedLPCM K = 2 with default prior settings
time_start <- Sys.time()
SS1Sce1_MixedLPCM_out_K2 <- VB_Directed_MixedLPCM(Y=SS1Sce1_MixedLPCM$Y,K=2,a=10,b=1,eta=1,rho2=1,omega2=1,xi=1,psi=1,delta=1,tol=1e-2)
time_end <- Sys.time()
SS1Sce1_MixedLPCM_out_K2$time <- time_end-time_start
SS1Sce1_MixedLPCM_out_K2$hat_z <- apply(SS1Sce1_MixedLPCM_out_K2$Pi_t,1,which.max)
SS1Sce1_MixedLPCM_out_K2$LS_hat_z <- LabelSwitching_SG2003(z=SS1Sce1_MixedLPCM_out_K2$hat_z)$z
# Store the PICL criteria value
SS1Sce1_MixedLPCM_out_K2$PICL <- 
  PICL_Directed_MixedLPCM(Y=SS1Sce1_MixedLPCM$Y,
                          hat_U=SS1Sce1_MixedLPCM_out_K2$U_t,
                          hat_V=SS1Sce1_MixedLPCM_out_K2$V_t,
                          hat_z=SS1Sce1_MixedLPCM_out_K2$LS_hat_z,
                          K=2,a=10,b=1,omega2=1,delta=1,tol=.Machine$double.xmin)

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

# Apply VB on the simulated MixedLPCM K = 3 with default prior settings
time_start <- Sys.time()
SS1Sce1_MixedLPCM_out_K3 <- VB_Directed_MixedLPCM(Y=SS1Sce1_MixedLPCM$Y,K=3,a=10,b=1,eta=1,rho2=1,omega2=1,xi=1,psi=1,delta=1,tol=1e-2)
time_end <- Sys.time()
SS1Sce1_MixedLPCM_out_K3$time <- time_end-time_start
SS1Sce1_MixedLPCM_out_K3$hat_z <- apply(SS1Sce1_MixedLPCM_out_K3$Pi_t,1,which.max)
SS1Sce1_MixedLPCM_out_K3$LS_hat_z <- LabelSwitching_SG2003(z=SS1Sce1_MixedLPCM_out_K3$hat_z)$z
# Store the PICL criteria value
SS1Sce1_MixedLPCM_out_K3$PICL <- 
  PICL_Directed_MixedLPCM(Y=SS1Sce1_MixedLPCM$Y,
                          hat_U=SS1Sce1_MixedLPCM_out_K3$U_t,
                          hat_V=SS1Sce1_MixedLPCM_out_K3$V_t,
                          hat_z=SS1Sce1_MixedLPCM_out_K3$LS_hat_z,
                          K=3,a=10,b=1,omega2=1,delta=1,tol=.Machine$double.xmin)

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

# Apply VB on the simulated MixedLPCM K = 4 with default prior settings
time_start <- Sys.time()
SS1Sce1_MixedLPCM_out_K4 <- VB_Directed_MixedLPCM(Y=SS1Sce1_MixedLPCM$Y,K=4,a=10,b=1,eta=1,rho2=1,omega2=1,xi=1,psi=1,delta=1,tol=1e-2)
time_end <- Sys.time()
SS1Sce1_MixedLPCM_out_K4$time <- time_end-time_start
SS1Sce1_MixedLPCM_out_K4$hat_z <- apply(SS1Sce1_MixedLPCM_out_K4$Pi_t,1,which.max)
SS1Sce1_MixedLPCM_out_K4$LS_hat_z <- LabelSwitching_SG2003(z=SS1Sce1_MixedLPCM_out_K4$hat_z)$z
# Store the PICL criteria value
SS1Sce1_MixedLPCM_out_K4$PICL <- 
  PICL_Directed_MixedLPCM(Y=SS1Sce1_MixedLPCM$Y,
                          hat_U=SS1Sce1_MixedLPCM_out_K4$U_t,
                          hat_V=SS1Sce1_MixedLPCM_out_K4$V_t,
                          hat_z=SS1Sce1_MixedLPCM_out_K4$LS_hat_z,
                          K=4,a=10,b=1,omega2=1,delta=1,tol=.Machine$double.xmin)

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

# Apply VB on the simulated MixedLPCM K = 4 with default prior settings but with smaller tol
time_start <- Sys.time()
SS1Sce1_MixedLPCM_out_K4_tol1em3 <- VB_Directed_MixedLPCM(Y=SS1Sce1_MixedLPCM$Y,K=4,a=10,b=1,eta=1,rho2=1,omega2=1,xi=1,psi=1,delta=1,tol=1e-3)
time_end <- Sys.time()
SS1Sce1_MixedLPCM_out_K4_tol1em3$time <- time_end-time_start
SS1Sce1_MixedLPCM_out_K4_tol1em3$hat_z <- apply(SS1Sce1_MixedLPCM_out_K4_tol1em3$Pi_t,1,which.max)
SS1Sce1_MixedLPCM_out_K4_tol1em3$LS_hat_z <- LabelSwitching_SG2003(z=SS1Sce1_MixedLPCM_out_K4_tol1em3$hat_z)$z
# Store the PICL criteria value
SS1Sce1_MixedLPCM_out_K4_tol1em3$PICL <- 
  PICL_Directed_MixedLPCM(Y=SS1Sce1_MixedLPCM$Y,
                          hat_U=SS1Sce1_MixedLPCM_out_K4_tol1em3$U_t,
                          hat_V=SS1Sce1_MixedLPCM_out_K4_tol1em3$V_t,
                          hat_z=SS1Sce1_MixedLPCM_out_K4_tol1em3$LS_hat_z,
                          K=4,a=10,b=1,omega2=1,delta=1,tol=.Machine$double.xmin)

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

# Apply VB on the simulated MixedLPCM K = 5 with default prior settings
time_start <- Sys.time()
SS1Sce1_MixedLPCM_out_K5 <- VB_Directed_MixedLPCM(Y=SS1Sce1_MixedLPCM$Y,K=5,a=10,b=1,eta=1,rho2=1,omega2=1,xi=1,psi=1,delta=1,tol=1e-2)
time_end <- Sys.time()
SS1Sce1_MixedLPCM_out_K5$time <- time_end-time_start
SS1Sce1_MixedLPCM_out_K5$hat_z <- apply(SS1Sce1_MixedLPCM_out_K5$Pi_t,1,which.max)
SS1Sce1_MixedLPCM_out_K5$LS_hat_z <- LabelSwitching_SG2003(z=SS1Sce1_MixedLPCM_out_K5$hat_z)$z
# Store the PICL criteria value
SS1Sce1_MixedLPCM_out_K5$PICL <- 
  PICL_Directed_MixedLPCM(Y=SS1Sce1_MixedLPCM$Y,
                          hat_U=SS1Sce1_MixedLPCM_out_K5$U_t,
                          hat_V=SS1Sce1_MixedLPCM_out_K5$V_t,
                          hat_z=SS1Sce1_MixedLPCM_out_K5$LS_hat_z,
                          K=5,a=10,b=1,omega2=1,delta=1,tol=.Machine$double.xmin)

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

# Apply VB on the simulated MixedLPCM K = 6 with default prior settings
time_start <- Sys.time()
SS1Sce1_MixedLPCM_out_K6 <- VB_Directed_MixedLPCM(Y=SS1Sce1_MixedLPCM$Y,K=6,a=10,b=1,eta=1,rho2=1,omega2=1,xi=1,psi=1,delta=1,tol=1e-2)
time_end <- Sys.time()
SS1Sce1_MixedLPCM_out_K6$time <- time_end-time_start
SS1Sce1_MixedLPCM_out_K6$hat_z <- apply(SS1Sce1_MixedLPCM_out_K6$Pi_t,1,which.max)
SS1Sce1_MixedLPCM_out_K6$LS_hat_z <- LabelSwitching_SG2003(z=SS1Sce1_MixedLPCM_out_K6$hat_z)$z
# Store the PICL criteria value
SS1Sce1_MixedLPCM_out_K6$PICL <- 
  PICL_Directed_MixedLPCM(Y=SS1Sce1_MixedLPCM$Y,
                          hat_U=SS1Sce1_MixedLPCM_out_K6$U_t,
                          hat_V=SS1Sce1_MixedLPCM_out_K6$V_t,
                          hat_z=SS1Sce1_MixedLPCM_out_K6$LS_hat_z,
                          K=6,a=10,b=1,omega2=1,delta=1,tol=.Machine$double.xmin)
```

#### 2.1.2 Fitting the PoisLPCM in Simulation Study 1 Scenario 1

The following code is for the Pois-LPCM fitting in simulation study 1 scenario 1.

``` r
# SS1 Apply PoisLPCM VB on the simulated MixedLPCM network K = 2 with default prior settings
time_start <- Sys.time()
SS1Sce1_PoisLPCM_out_K2 <- VB_Directed_PoisLPCM(Y=SS1Sce1_MixedLPCM$Y,K=2,eta=1,rho2=1,omega2=1,xi=1,psi=1,delta=1,tol=1e-2)
time_end <- Sys.time()
SS1Sce1_PoisLPCM_out_K2$time <- time_end-time_start
SS1Sce1_PoisLPCM_out_K2$hat_z <- apply(SS1Sce1_PoisLPCM_out_K2$Pi_t,1,which.max)
SS1Sce1_PoisLPCM_out_K2$LS_hat_z <- LabelSwitching_SG2003(z=SS1Sce1_PoisLPCM_out_K2$hat_z)$z
# Store the PICL criteria value
SS1Sce1_PoisLPCM_out_K2$PICL <-
  PICL_Directed_PoisLPCM(Y=SS1Sce1_MixedLPCM$Y,
                         hat_U=SS1Sce1_PoisLPCM_out_K2$U_t,
                         hat_z=SS1Sce1_PoisLPCM_out_K2$LS_hat_z,
                         K=2,omega2=1,delta=1,tol=.Machine$double.xmin)

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

# SS1 Apply PoisLPCM VB on the simulated MixedLPCM network K = 3 with default prior settings
time_start <- Sys.time()
SS1Sce1_PoisLPCM_out_K3 <- VB_Directed_PoisLPCM(Y=SS1Sce1_MixedLPCM$Y,K=3,eta=1,rho2=1,omega2=1,xi=1,psi=1,delta=1,tol=1e-2)
time_end <- Sys.time()
SS1Sce1_PoisLPCM_out_K3$time <- time_end-time_start
SS1Sce1_PoisLPCM_out_K3$hat_z <- apply(SS1Sce1_PoisLPCM_out_K3$Pi_t,1,which.max)
SS1Sce1_PoisLPCM_out_K3$LS_hat_z <- LabelSwitching_SG2003(z=SS1Sce1_PoisLPCM_out_K3$hat_z)$z
# Store the PICL criteria value
SS1Sce1_PoisLPCM_out_K3$PICL <-
  PICL_Directed_PoisLPCM(Y=SS1Sce1_MixedLPCM$Y,
                         hat_U=SS1Sce1_PoisLPCM_out_K3$U_t,
                         hat_z=SS1Sce1_PoisLPCM_out_K3$LS_hat_z,
                         K=3,omega2=1,delta=1,tol=.Machine$double.xmin)

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------


# SS1 Apply PoisLPCM VB on the simulated MixedLPCM network K = 4 with default prior settings
time_start <- Sys.time()
SS1Sce1_PoisLPCM_out_K4 <- VB_Directed_PoisLPCM(Y=SS1Sce1_MixedLPCM$Y,K=4,eta=1,rho2=1,omega2=1,xi=1,psi=1,delta=1,tol=1e-2)
time_end <- Sys.time()
SS1Sce1_PoisLPCM_out_K4$time <- time_end-time_start
SS1Sce1_PoisLPCM_out_K4$hat_z <- apply(SS1Sce1_PoisLPCM_out_K4$Pi_t,1,which.max)
SS1Sce1_PoisLPCM_out_K4$LS_hat_z <- LabelSwitching_SG2003(z=SS1Sce1_PoisLPCM_out_K4$hat_z)$z
# Store the PICL criteria value
SS1Sce1_PoisLPCM_out_K4$PICL <-
  PICL_Directed_PoisLPCM(Y=SS1Sce1_MixedLPCM$Y,
                         hat_U=SS1Sce1_PoisLPCM_out_K4$U_t,
                         hat_z=SS1Sce1_PoisLPCM_out_K4$LS_hat_z,
                         K=4,omega2=1,delta=1,tol=.Machine$double.xmin)

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

# SS1 Apply PoisLPCM VB on the simulated MixedLPCM network K = 5 with default prior settings
time_start <- Sys.time()
SS1Sce1_PoisLPCM_out_K5 <- VB_Directed_PoisLPCM(Y=SS1Sce1_MixedLPCM$Y,K=5,eta=1,rho2=1,omega2=1,xi=1,psi=1,delta=1,tol=1e-2)
time_end <- Sys.time()
SS1Sce1_PoisLPCM_out_K5$time <- time_end-time_start
SS1Sce1_PoisLPCM_out_K5$hat_z <- apply(SS1Sce1_PoisLPCM_out_K5$Pi_t,1,which.max)
SS1Sce1_PoisLPCM_out_K5$LS_hat_z <- LabelSwitching_SG2003(z=SS1Sce1_PoisLPCM_out_K5$hat_z)$z
# Store the PICL criteria value
SS1Sce1_PoisLPCM_out_K5$PICL <-
  PICL_Directed_PoisLPCM(Y=SS1Sce1_MixedLPCM$Y,
                         hat_U=SS1Sce1_PoisLPCM_out_K5$U_t,
                         hat_z=SS1Sce1_PoisLPCM_out_K5$LS_hat_z,
                         K=5,omega2=1,delta=1,tol=.Machine$double.xmin)

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

# SS1 Apply PoisLPCM VB on the simulated MixedLPCM network K = 6 with default prior settings
time_start <- Sys.time()
SS1Sce1_PoisLPCM_out_K6 <- VB_Directed_PoisLPCM(Y=SS1Sce1_MixedLPCM$Y,K=6,eta=1,rho2=1,omega2=1,xi=1,psi=1,delta=1,tol=1e-2)
time_end <- Sys.time()
SS1Sce1_PoisLPCM_out_K6$time <- time_end-time_start
SS1Sce1_PoisLPCM_out_K6$hat_z <- apply(SS1Sce1_PoisLPCM_out_K6$Pi_t,1,which.max)
SS1Sce1_PoisLPCM_out_K6$LS_hat_z <- LabelSwitching_SG2003(z=SS1Sce1_PoisLPCM_out_K6$hat_z)$z
# Store the PICL criteria value
SS1Sce1_PoisLPCM_out_K6$PICL <-
  PICL_Directed_PoisLPCM(Y=SS1Sce1_MixedLPCM$Y,
                         hat_U=SS1Sce1_PoisLPCM_out_K6$U_t,
                         hat_z=SS1Sce1_PoisLPCM_out_K6$LS_hat_z,
                         K=6,omega2=1,delta=1,tol=.Machine$double.xmin)
```


#### 2.1.3 Summarizing the Output in Simulation Study 1 Scenario 1

Taking the MLPCM K=4 case as an example here, the following code provides the summarizing process after the model fitting.

``` r
# # Apply MixedLPCM VB on the simulated MixedLPCM K = 4 with default prior settings
# time_start <- Sys.time()
# SS1Sce1_MixedLPCM_out_K4 <- VB_Directed_MixedLPCM(Y=SS1Sce1_MixedLPCM$Y,K=4,a=10,b=1,eta=1,rho2=1,omega2=1,xi=1,psi=1,delta=1,tol=1e-2)
# time_end <- Sys.time()
# SS1Sce1_MixedLPCM_out_K4$time <- time_end-time_start
# SS1Sce1_MixedLPCM_out_K4$hat_z <- apply(SS1Sce1_MixedLPCM_out_K4$Pi_t,1,which.max)
# SS1Sce1_MixedLPCM_out_K4$LS_hat_z <- LabelSwitching_SG2003(z=SS1Sce1_MixedLPCM_out_K4$hat_z)$z
# # Store the PICL criteria value
# SS1Sce1_MixedLPCM_out_K4$PICL <- 
#   PICL_Directed_MixedLPCM(Y=SS1Sce1_MixedLPCM$Y,
#                           hat_U=SS1Sce1_MixedLPCM_out_K4$U_t,
#                           hat_V=SS1Sce1_MixedLPCM_out_K4$V_t,
#                           hat_z=SS1Sce1_MixedLPCM_out_K4$LS_hat_z,
#                           K=4,a=10,b=1,omega2=1,delta=1,tol=.Machine$double.xmin)
# SS1Sce1_MixedLPCM_out_K4$PICL # 14724

# Check time, iterations and clustering
SS1Sce1_MixedLPCM_out_K4$time # Time difference of 15.26885 mins
SS1Sce1_MixedLPCM_out_K4$NumIt # 3753
table(apply(SS1Sce1_MixedLPCM_out_K4$Pi_t,1,which.max),SS1Sce1_MixedLPCM$z) # Check clustering
library("mcclust")
vi.dist(SS1Sce1_MixedLPCM_out_K4$hat_z,SS1Sce1_MixedLPCM$z) # 0.1217233
SS1Sce1_MixedLPCM_out_K4$ELBO_list[length(SS1Sce1_MixedLPCM_out_K4$ELBO_list)] # -11971.87
# Check estimated LPM intercept eta_t for the model parameter beta
SS1Sce1_MixedLPCM_out_K4$eta_t # 0.9655156
SS1Sce1_MixedLPCM_out_K4$rho2_t #0.0002965032
SS1Sce1_MixedLPCM_beta # 1

# Compare the pairwise distance squares
hist(abs(apply((array(t(SS1Sce1_MixedLPCM$U),dim=c(dim(t(SS1Sce1_MixedLPCM$U)),nrow(SS1Sce1_MixedLPCM$U)))-SS1Sce1_MixedLPCM$V)^2,2:3,sum)-
           apply((array(t(SS1Sce1_MixedLPCM_out_K4$U_t),dim=c(dim(t(SS1Sce1_MixedLPCM_out_K4$U_t)),nrow(SS1Sce1_MixedLPCM_out_K4$U_t)))-SS1Sce1_MixedLPCM_out_K4$V_t)^2,2:3,sum)),
     breaks=1000,main="Pairwise DistSquare AbsError",xlab="")
# Mean absolute error
mean(c(abs(apply((array(t(SS1Sce1_MixedLPCM$U),dim=c(dim(t(SS1Sce1_MixedLPCM$U)),nrow(SS1Sce1_MixedLPCM$U)))-SS1Sce1_MixedLPCM$V)^2,2:3,sum)-
             apply((array(t(SS1Sce1_MixedLPCM_out_K4$U_t),dim=c(dim(t(SS1Sce1_MixedLPCM_out_K4$U_t)),nrow(SS1Sce1_MixedLPCM_out_K4$U_t)))-SS1Sce1_MixedLPCM_out_K4$V_t)^2,2:3,sum)))[
               !is.na(c(abs(apply((array(t(SS1Sce1_MixedLPCM$U),dim=c(dim(t(SS1Sce1_MixedLPCM$U)),nrow(SS1Sce1_MixedLPCM$U)))-SS1Sce1_MixedLPCM$V)^2,2:3,sum)-
                              apply((array(t(SS1Sce1_MixedLPCM_out_K4$U_t),dim=c(dim(t(SS1Sce1_MixedLPCM_out_K4$U_t)),nrow(SS1Sce1_MixedLPCM_out_K4$U_t)))-SS1Sce1_MixedLPCM_out_K4$V_t)^2,2:3,sum))))])
# 1.576048
sd(c(abs(apply((array(t(SS1Sce1_MixedLPCM$U),dim=c(dim(t(SS1Sce1_MixedLPCM$U)),nrow(SS1Sce1_MixedLPCM$U)))-SS1Sce1_MixedLPCM$V)^2,2:3,sum)-
           apply((array(t(SS1Sce1_MixedLPCM_out_K4$U_t),dim=c(dim(t(SS1Sce1_MixedLPCM_out_K4$U_t)),nrow(SS1Sce1_MixedLPCM_out_K4$U_t)))-SS1Sce1_MixedLPCM_out_K4$V_t)^2,2:3,sum)))[
             !is.na(c(abs(apply((array(t(SS1Sce1_MixedLPCM$U),dim=c(dim(t(SS1Sce1_MixedLPCM$U)),nrow(SS1Sce1_MixedLPCM$U)))-SS1Sce1_MixedLPCM$V)^2,2:3,sum)-
                            apply((array(t(SS1Sce1_MixedLPCM_out_K4$U_t),dim=c(dim(t(SS1Sce1_MixedLPCM_out_K4$U_t)),nrow(SS1Sce1_MixedLPCM_out_K4$U_t)))-SS1Sce1_MixedLPCM_out_K4$V_t)^2,2:3,sum))))])
# 1.602655

# Check variational MVN variance sigma2_t for U_t
hist(SS1Sce1_MixedLPCM_out_K4$sigma2_t,breaks=100,main=TeX(r"(R1 $\tilde{\sigma}^2$)",bold=TRUE),xlab="")
mean(SS1Sce1_MixedLPCM_out_K4$sigma2_t) # 0.001005315
sd(SS1Sce1_MixedLPCM_out_K4$sigma2_t) # 8.659065e-05
# Check variational MVN variance omega2_t for mu_t
SS1Sce1_MixedLPCM_out_K4$omega2_t
SS1Sce1_MixedLPCM_out_K4$omega2_t[c(4,2,3,1)] # 0.004427848 0.007414585 0.013428057 0.019938164

# Check the estimated overt latent positions
SS1Sce1_MixedLPCM_obs_gamma <- c() # calculate the gamma_j based on the observed V*
SS1Sce1_MixedLPCM_out_K4_V_t_precision <- c() # obtain the precision (1/var) of the {v_t_ij:i=1,2,\dots,N; i != j} for each covert latent position u_j
for (i in 1:nrow(SS1Sce1_MixedLPCM$Y)){
  SS1Sce1_MixedLPCM_obs_gamma <- c(SS1Sce1_MixedLPCM_obs_gamma,1/mean(apply(SS1Sce1_MixedLPCM$V[,-i,i],1,var)))
  SS1Sce1_MixedLPCM_out_K4_V_t_precision <- c(SS1Sce1_MixedLPCM_out_K4_V_t_precision,1/mean(apply(SS1Sce1_MixedLPCM_out_K4$V_t[,-i,i],1,var)))
}
hist(1/SS1Sce1_MixedLPCM_obs_gamma,breaks=100,main=TeX(r"(Obs var$\{v^*_{j\leftarrow i}:i;i\neq j\}_j$ var)",bold=TRUE),xlab="")
hist(1/SS1Sce1_MixedLPCM_out_K4_V_t_precision,breaks=100,main=TeX(r"(R1 $\{\tilde{v}_{j\leftarrow i}:i;i\neq j\}_j$ var)",bold=TRUE),xlab="")
# Check the variational MVN var varphi2_t for the variational overt latent positions
hist(SS1Sce1_MixedLPCM_out_K4$varphi2_t,breaks=100,main=TeX(r"(R1 $\tilde{\varphi}^2$)",bold=TRUE),xlab="")
mean(SS1Sce1_MixedLPCM_out_K4$varphi2_t[!is.na(SS1Sce1_MixedLPCM_out_K4$varphi2_t)]) # 0.1010904
sd(SS1Sce1_MixedLPCM_out_K4$varphi2_t[!is.na(SS1Sce1_MixedLPCM_out_K4$varphi2_t)]) # 0.01315514

# Check estimated overt latent positions' precision, gamma
SS1Sce1_MixedLPCM_out_K4$a_t
SS1Sce1_MixedLPCM_out_K4$b_t
SS1Sce1_MixedLPCM_out_K4$a_t/SS1Sce1_MixedLPCM_out_K4$b_t
mean(SS1Sce1_MixedLPCM_out_K4$b_t/SS1Sce1_MixedLPCM_out_K4$a_t) # 0.1040671
sd(SS1Sce1_MixedLPCM_out_K4$b_t/SS1Sce1_MixedLPCM_out_K4$a_t) # 0.009399977
SS1Sce1_MixedLPCM_out_K4$a_t/(SS1Sce1_MixedLPCM_out_K4$b_t^2)
SS1Sce1_MixedLPCM_gamma

# Check estimated Dirichlet parameters, delta_t
SS1Sce1_MixedLPCM_out_K4$delta_t
SS1Sce1_MixedLPCM_out_K4$delta_t[c(4,2,3,1)] # 25.94473 26.63933 25.63417 25.78198
# Check estimated covert latent positions' group precision, tau_t
SS1Sce1_MixedLPCM_out_K4$xi_t
SS1Sce1_MixedLPCM_out_K4$psi_t
SS1Sce1_MixedLPCM_out_K4$xi_t/SS1Sce1_MixedLPCM_out_K4$psi_t
SS1Sce1_MixedLPCM_out_K4$xi_t/(SS1Sce1_MixedLPCM_out_K4$psi_t^2)
SS1Sce1_MixedLPCM_out_K4$xi_t[c(4,2,3,1)]
SS1Sce1_MixedLPCM_out_K4$psi_t[c(4,2,3,1)]
(SS1Sce1_MixedLPCM_out_K4$xi_t/SS1Sce1_MixedLPCM_out_K4$psi_t)[c(4,2,3,1)] # 9.013482 5.221144 2.982509 1.983509
(SS1Sce1_MixedLPCM_out_K4$xi_t/(SS1Sce1_MixedLPCM_out_K4$psi_t^2))[c(4,2,3,1)] # 3.1313876 1.0233158 0.3470113 0.1525998
SS1Sce1_MixedLPCM_tau # 8 6 4 2 # reference precision

#---------------------------------------------------------------------------------------------------------------
# plot the latent positions
g_obs <- graph_from_adjacency_matrix(SS1Sce1_MixedLPCM$Y,mode = "directed",weighted=TRUE)
# Node sizes are proportional to their social impact, \tilde{\gamma}_i
V(g_obs)$size <- 20*(1/(SS1Sce1_MixedLPCM_out_K4$a_t/SS1Sce1_MixedLPCM_out_K4$b_t))/max(1/(SS1Sce1_MixedLPCM_out_K4$a_t/SS1Sce1_MixedLPCM_out_K4$b_t))
# Additional graphical settings
V(g_obs)$frame.color <- c("black","red")[c((SS1Sce1_MixedLPCM_gamma==1)+1)]
V(g_obs)$label <- ""
# Node colors indicate the inferred clustering
Customized_colors <- My_colors[6:9]
# Edge colors proportional to edge weights
E(g_obs)$color <- colorRampPalette(brewer.pal(9,"Greys")[c(3,9)])(max(E(g_obs)$weight))[E(g_obs)$weight]
# Make the position plot
SS1Sce1_MixedLPCM_out_K4$Pi_t_pie <- list()
for (i in 1:nrow(SS1Sce1_MixedLPCM_out_K4$Pi_t)){
  SS1Sce1_MixedLPCM_out_K4$Pi_t_pie[[i]] <- SS1Sce1_MixedLPCM_out_K4$Pi_t[i,]
}
par(mfrow=c(1,1),mai = c(0.05, 0.05, 0.05, 0.05),mgp = c(1,0.25,0))
plot(g_obs,rescale=F, layout=Procrustes(SS1Sce1_MixedLPCM_out_K4$U_t,SS1Sce1_MixedLPCM$U, translate = TRUE ,dilate = FALSE)$X.new,
     edge.curved=0.2,vertex.shape="pie",
     vertex.pie=SS1Sce1_MixedLPCM_out_K4$Pi_t_pie,
     vertex.pie.color=list(adjustcolor(Customized_colors[c(4,2,3,1)], alpha.f = .5)),
     edge.arrow.size=0.055,xlim=c(-3,3),ylim=c(-3,3),edge.width=E(g_obs)$weight*0.25)
par(mfrow=c(1,1),mai = c(1.02, 0.82, 0.82, 0.42))
```

#### 2.1.4 Multiple implementations with Random initializations in Simulation Study 1 Scenario 1

Recall here that we perform a robustness checking experiment by fitting the MixedLPCM in simulation study 1 scenario 1 for 10 independent runs with moderate random initializations around the default settings.

``` r
# Apply MixedLPCM VB on the simulated MixdLPCM K = 4 with default prior settings and moderate random initializations
set.seed(123)
SS1Sce1_MixedLPCM_out_K4_sigma2_t_RandomInit <- sample(c(0.2,0.4,0.6,0.8,1,1.25,1.5,1.75,2))
SS1Sce1_MixedLPCM_out_K4_varphi2_t_RandomInit <- sample(c(0.2,0.4,0.6,0.8,1,1.25,1.5,1.75,2))
SS1Sce1_MixedLPCM_out_K4_eta_t_RandomInit <- rnorm(10,1,1)
SS1Sce1_MixedLPCM_out_K4_rho2_t_RandomInit <- sample(c(0.2,0.4,0.6,0.8,1,1.25,1.5,1.75,2))
SS1Sce1_MixedLPCM_out_K4_omega2_t_RandomInit <- sample(c(0.2,0.4,0.6,0.8,1,1.25,1.5,1.75,2))
SS1Sce1_MixedLPCM_out_K4_xi_t_RandomInit <- sample(c(0.2,0.4,0.6,0.8,1,1.25,1.5,1.75,2))
SS1Sce1_MixedLPCM_out_K4_psi_t_RandomInit <- sample(c(0.2,0.4,0.6,0.8,1,1.25,1.5,1.75,2))
SS1Sce1_MixedLPCM_out_K4_a_t_RandomInit <- sample(seq(2,20,2))
SS1Sce1_MixedLPCM_out_K4_b_t_RandomInit <- sample(c(0.2,0.4,0.6,0.8,1,1.25,1.5,1.75,2))
SS1Sce1_MixedLPCM_out_K4_delta_t_RandomInit <- sample(c(0.2,0.4,0.6,0.8,1,1.25,1.5,1.75,2))
set.seed(NULL)

SS1Sce1_MixedLPCM_out_K4_RandomInit <- list()
SS1Sce1_MixedLPCM_out_K4_N <- nrow(SS1Sce1_MixedLPCM$Y)
SS1Sce1_MixedLPCM_out_K4_K <- 4
for (R in 1:10){
  cat("Round", R,"\n")
  # Initialize the states
  SS1Sce1_MixedLPCM_out_K4_sigma2_t <- rep(SS1Sce1_MixedLPCM_out_K4_sigma2_t_RandomInit[R],SS1Sce1_MixedLPCM_out_K4_N)
  SS1Sce1_MixedLPCM_out_K4_varphi2_t <- matrix(1,SS1Sce1_MixedLPCM_out_K4_N,SS1Sce1_MixedLPCM_out_K4_N)
  diag(SS1Sce1_MixedLPCM_out_K4_varphi2_t) <- NA
  SS1Sce1_MixedLPCM_out_K4_eta_t <- SS1Sce1_MixedLPCM_out_K4_eta_t_RandomInit[R]
  SS1Sce1_MixedLPCM_out_K4_rho2_t <- SS1Sce1_MixedLPCM_out_K4_rho2_t_RandomInit[R]
  SS1Sce1_MixedLPCM_out_K4_omega2_t <- rep(SS1Sce1_MixedLPCM_out_K4_omega2_t_RandomInit[R],SS1Sce1_MixedLPCM_out_K4_K)
  SS1Sce1_MixedLPCM_out_K4_xi_t <- rep(SS1Sce1_MixedLPCM_out_K4_xi_t_RandomInit[R],SS1Sce1_MixedLPCM_out_K4_K)
  SS1Sce1_MixedLPCM_out_K4_psi_t <- rep(SS1Sce1_MixedLPCM_out_K4_psi_t_RandomInit[R],SS1Sce1_MixedLPCM_out_K4_K)
  SS1Sce1_MixedLPCM_out_K4_a_t <- rep(SS1Sce1_MixedLPCM_out_K4_a_t_RandomInit[R],SS1Sce1_MixedLPCM_out_K4_N)
  SS1Sce1_MixedLPCM_out_K4_b_t <- rep(SS1Sce1_MixedLPCM_out_K4_b_t_RandomInit[R],SS1Sce1_MixedLPCM_out_K4_N)
  SS1Sce1_MixedLPCM_out_K4_delta_t <- rep(SS1Sce1_MixedLPCM_out_K4_delta_t_RandomInit[R],SS1Sce1_MixedLPCM_out_K4_K)
  # Implementation
  time_start <- Sys.time()
  SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]] <- 
    VB_Directed_MixedLPCM(Y=SS1Sce1_MixedLPCM$Y,K=4,a=10,b=1,eta=1,rho2=1,omega2=1,xi=1,psi=1,delta=1,tol=1e-2,
                          sigma2_t=SS1Sce1_MixedLPCM_out_K4_sigma2_t,
                          varphi2_t=SS1Sce1_MixedLPCM_out_K4_varphi2_t,
                          eta_t=SS1Sce1_MixedLPCM_out_K4_eta_t,
                          rho2_t=SS1Sce1_MixedLPCM_out_K4_rho2_t,
                          omega2_t=SS1Sce1_MixedLPCM_out_K4_omega2_t,
                          xi_t=SS1Sce1_MixedLPCM_out_K4_xi_t,
                          psi_t=SS1Sce1_MixedLPCM_out_K4_psi_t,
                          a_t=SS1Sce1_MixedLPCM_out_K4_a_t,
                          b_t=SS1Sce1_MixedLPCM_out_K4_b_t,
                          delta_t=SS1Sce1_MixedLPCM_out_K4_delta_t,)
  time_end <- Sys.time()
  SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$time <- time_end-time_start
  SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$hat_z <- apply(SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$Pi_t,1,which.max)
  SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$LS_hat_z <- LabelSwitching_SG2003(z=SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$hat_z)$z
  # Store the PICL criteria value
  SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$PICL <- 
    PICL_Directed_MixedLPCM(Y=SS1Sce1_MixedLPCM$Y,
                            hat_U=SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$U_t,
                            hat_V=SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$V_t,
                            hat_z=SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$LS_hat_z,
                            K=4,a=10,b=1,omega2=1,delta=1,tol=.Machine$double.xmin)
}

SS1Sce1_MixedLPCM_out_K4_RandomInit_PICL_List <- c()
SS1Sce1_MixedLPCM_out_K4_RandomInit_time_List <- c()
SS1Sce1_MixedLPCM_out_K4_RandomInit_NumIt_List <- c()
library("mcclust")
SS1Sce1_MixedLPCM_out_K4_RandomInit_VIdist_List <- c()
SS1Sce1_MixedLPCM_out_K4_RandomInit_ELBO_list <- c()
SS1Sce1_MixedLPCM_out_K4_RandomInit_eta_t_list <- c()
SS1Sce1_MixedLPCM_out_K4_RandomInit_rho2_t_list <- c()
SS1Sce1_MixedLPCM_out_K4_RandomInit_Esigma2_t_list <- c()
SS1Sce1_MixedLPCM_out_K4_RandomInit_sdsigma2_t_list <- c()
SS1Sce1_MixedLPCM_out_K4_RandomInit_Evarphi2_t_list <- c()
SS1Sce1_MixedLPCM_out_K4_RandomInit_sdvarphi2_t_list <- c()
SS1Sce1_MixedLPCM_out_K4_RandomInit_MAEdistsquare_list <- c()
SS1Sce1_MixedLPCM_out_K4_RandomInit_sdAEdistsquare_list <- c()
SS1Sce1_MixedLPCM_out_K4_RandomInit_E_a_t_div_b_t_list <- c()
SS1Sce1_MixedLPCM_out_K4_RandomInit_sd_a_t_div_b_t_list <- c()
SS1Sce1_MixedLPCM_out_K4_RandomInit_E_b_t_div_a_t_list <- c()
SS1Sce1_MixedLPCM_out_K4_RandomInit_sd_b_t_div_a_t_list <- c()
for (R in 1:10){
  # diag(SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$varphi2_t) <- NA
  SS1Sce1_MixedLPCM_out_K4_RandomInit_PICL_List <-
    c(SS1Sce1_MixedLPCM_out_K4_RandomInit_PICL_List,SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$PICL$PICL)
  SS1Sce1_MixedLPCM_out_K4_RandomInit_time_List <- 
    c(SS1Sce1_MixedLPCM_out_K4_RandomInit_time_List,SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$time)
  SS1Sce1_MixedLPCM_out_K4_RandomInit_NumIt_List <- 
    c(SS1Sce1_MixedLPCM_out_K4_RandomInit_NumIt_List,SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$NumIt)
  SS1Sce1_MixedLPCM_out_K4_RandomInit_VIdist_List <- 
    c(SS1Sce1_MixedLPCM_out_K4_RandomInit_VIdist_List,vi.dist(SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$hat_z,SS1Sce1_MixedLPCM$z))
  SS1Sce1_MixedLPCM_out_K4_RandomInit_ELBO_list <- 
    c(SS1Sce1_MixedLPCM_out_K4_RandomInit_ELBO_list,SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$ELBO_list[length(SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$ELBO_list)])
  SS1Sce1_MixedLPCM_out_K4_RandomInit_eta_t_list <- 
    c(SS1Sce1_MixedLPCM_out_K4_RandomInit_eta_t_list,SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$eta_t)
  SS1Sce1_MixedLPCM_out_K4_RandomInit_rho2_t_list <- 
    c(SS1Sce1_MixedLPCM_out_K4_RandomInit_rho2_t_list,SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$rho2_t)
  SS1Sce1_MixedLPCM_out_K4_RandomInit_Esigma2_t_list <- 
    c(SS1Sce1_MixedLPCM_out_K4_RandomInit_Esigma2_t_list,mean(SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$sigma2_t))
  SS1Sce1_MixedLPCM_out_K4_RandomInit_sdsigma2_t_list <- 
    c(SS1Sce1_MixedLPCM_out_K4_RandomInit_sdsigma2_t_list,sd(SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$sigma2_t))
  SS1Sce1_MixedLPCM_out_K4_RandomInit_Evarphi2_t_list <- 
    c(SS1Sce1_MixedLPCM_out_K4_RandomInit_Evarphi2_t_list,mean(SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$varphi2_t,na.rm=TRUE))
  SS1Sce1_MixedLPCM_out_K4_RandomInit_sdvarphi2_t_list <- 
    c(SS1Sce1_MixedLPCM_out_K4_RandomInit_sdvarphi2_t_list,sd(SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$varphi2_t,na.rm=TRUE))
  SS1Sce1_MixedLPCM_out_K4_RandomInit_MAEdistsquare_list <- 
    c(SS1Sce1_MixedLPCM_out_K4_RandomInit_MAEdistsquare_list,
      mean(c(abs(apply((array(t(SS1Sce1_MixedLPCM$U),dim=c(dim(t(SS1Sce1_MixedLPCM$U)),nrow(SS1Sce1_MixedLPCM$U)))-SS1Sce1_MixedLPCM$V)^2,2:3,sum)-
                   apply((array(t(SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$U_t),dim=c(dim(t(SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$U_t)),
                                                                                          nrow(SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$U_t)))-
                            SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$V_t)^2,2:3,sum))),na.rm=TRUE))
  SS1Sce1_MixedLPCM_out_K4_RandomInit_sdAEdistsquare_list <- 
    c(SS1Sce1_MixedLPCM_out_K4_RandomInit_sdAEdistsquare_list,
      sd(c(abs(apply((array(t(SS1Sce1_MixedLPCM$U),dim=c(dim(t(SS1Sce1_MixedLPCM$U)),nrow(SS1Sce1_MixedLPCM$U)))-SS1Sce1_MixedLPCM$V)^2,2:3,sum)-
                 apply((array(t(SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$U_t),dim=c(dim(t(SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$U_t)),
                                                                                        nrow(SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$U_t)))-
                          SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$V_t)^2,2:3,sum))),na.rm=TRUE))
  SS1Sce1_MixedLPCM_out_K4_RandomInit_E_a_t_div_b_t_list <- 
    c(SS1Sce1_MixedLPCM_out_K4_RandomInit_E_a_t_div_b_t_list,mean(SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$a_t/SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$b_t))
  SS1Sce1_MixedLPCM_out_K4_RandomInit_sd_a_t_div_b_t_list <- 
    c(SS1Sce1_MixedLPCM_out_K4_RandomInit_sd_a_t_div_b_t_list,sd(SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$a_t/SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$b_t))
  
  SS1Sce1_MixedLPCM_out_K4_RandomInit_E_b_t_div_a_t_list <- 
    c(SS1Sce1_MixedLPCM_out_K4_RandomInit_E_b_t_div_a_t_list,mean(SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$b_t/SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$a_t))
  SS1Sce1_MixedLPCM_out_K4_RandomInit_sd_b_t_div_a_t_list <- 
    c(SS1Sce1_MixedLPCM_out_K4_RandomInit_sd_b_t_div_a_t_list,sd(SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$b_t/SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$a_t))
}



# Check PICL
SS1Sce1_MixedLPCM_out_K4_RandomInit_PICL_List
boxplot(SS1Sce1_MixedLPCM_out_K4_RandomInit_PICL_List)
max(SS1Sce1_MixedLPCM_out_K4_RandomInit_PICL_List) # 14977.25
median(SS1Sce1_MixedLPCM_out_K4_RandomInit_PICL_List) # 14752.55
min(SS1Sce1_MixedLPCM_out_K4_RandomInit_PICL_List) # 11175.59
sd(SS1Sce1_MixedLPCM_out_K4_RandomInit_PICL_List) # 1140.57

# Check time
SS1Sce1_MixedLPCM_out_K4_RandomInit_time_List[c(7,9,10)] <- SS1Sce1_MixedLPCM_out_K4_RandomInit_time_List[c(7,9,10)]*60
SS1Sce1_MixedLPCM_out_K4_RandomInit_time_List

# Check number of iterations
SS1Sce1_MixedLPCM_out_K4_RandomInit_NumIt_List

# Check vi.dist between hat_z and z*
SS1Sce1_MixedLPCM_out_K4_RandomInit_VIdist_List
max(SS1Sce1_MixedLPCM_out_K4_RandomInit_VIdist_List[-1]) # 0.2428694
median(SS1Sce1_MixedLPCM_out_K4_RandomInit_VIdist_List[-1]) # 0.1217233
min(SS1Sce1_MixedLPCM_out_K4_RandomInit_VIdist_List[-1]) # 0.1217233
sd(SS1Sce1_MixedLPCM_out_K4_RandomInit_VIdist_List[-1]) # 0.04038203

table(apply(SS1Sce1_MixedLPCM_out_K4_RandomInit[[1]]$Pi_t,1,which.max),SS1Sce1_MixedLPCM$z) # Check clustering
table(apply(SS1Sce1_MixedLPCM_out_K4_RandomInit[[7]]$Pi_t,1,which.max),SS1Sce1_MixedLPCM$z) # Check clustering

# Check distinct hat_z
table(apply(SS1Sce1_MixedLPCM_out_K4_RandomInit[[1]]$Pi_t,1,which.max),SS1Sce1_MixedLPCM$z) # Check clustering
table(apply(SS1Sce1_MixedLPCM_out_K4_RandomInit[[2]]$Pi_t,1,which.max),SS1Sce1_MixedLPCM$z) # Check clustering
table(apply(SS1Sce1_MixedLPCM_out_K4_RandomInit[[7]]$Pi_t,1,which.max),SS1Sce1_MixedLPCM$z) # Check clustering
table(apply(SS1Sce1_MixedLPCM_out_K4_RandomInit[[9]]$Pi_t,1,which.max),SS1Sce1_MixedLPCM$z) # Check clustering
table(apply(SS1Sce1_MixedLPCM_out_K4$Pi_t,1,which.max),SS1Sce1_MixedLPCM$z) # Check clustering

# Check ELBO
SS1Sce1_MixedLPCM_out_K4_RandomInit_ELBO_list
max(SS1Sce1_MixedLPCM_out_K4_RandomInit_ELBO_list) # -11969.67
median(SS1Sce1_MixedLPCM_out_K4_RandomInit_ELBO_list) # -11996.97
min(SS1Sce1_MixedLPCM_out_K4_RandomInit_ELBO_list) # -12760.46
sd(SS1Sce1_MixedLPCM_out_K4_RandomInit_ELBO_list) # 315.2172
SS1Sce1_MixedLPCM_out_K4$ELBO_list[length(SS1Sce1_MixedLPCM_out_K4$ELBO_list)] # -11971.87

# Check eta_t
SS1Sce1_MixedLPCM_out_K4_RandomInit_eta_t_list
max(SS1Sce1_MixedLPCM_out_K4_RandomInit_eta_t_list[-1]) # 1.079989
median(SS1Sce1_MixedLPCM_out_K4_RandomInit_eta_t_list[-1]) # 0.8781608
min(SS1Sce1_MixedLPCM_out_K4_RandomInit_eta_t_list[-1]) # 0.8209029
sd(SS1Sce1_MixedLPCM_out_K4_RandomInit_eta_t_list[-1]) # 0.09156729

# Check rho2_t
SS1Sce1_MixedLPCM_out_K4_RandomInit_rho2_t_list
max(SS1Sce1_MixedLPCM_out_K4_RandomInit_rho2_t_list[-1]) # 0.001096851
median(SS1Sce1_MixedLPCM_out_K4_RandomInit_rho2_t_list[-1]) # 0.0003048629
min(SS1Sce1_MixedLPCM_out_K4_RandomInit_rho2_t_list[-1]) # 0.0002831225
sd(SS1Sce1_MixedLPCM_out_K4_RandomInit_rho2_t_list[-1]) # 0.0002641932

# Check mean pairwise distance between mu_t and mu*
SS1Sce1_MixedLPCM_out_K4_RandomInit_muMeanDist_list <- c(NA)
SS1Sce1_MixedLPCM_out_K4_RandomInit_musdDist_list <- c(NA)
for (R in 2:10){
  SS1Sce1_MixedLPCM_out_K4_RandomInit_muMeanDist_list <- 
    c(SS1Sce1_MixedLPCM_out_K4_RandomInit_muMeanDist_list,
      mean(sqrt(rowSums((Procrustes(
        t(LabelSwitching_SG2003(z=SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$hat_z,
                                matrix_dbyK=t(SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$mu_t))$matrix_dbyK),
        SS1Sce1_MixedLPCM_mu, translate = TRUE ,dilate = FALSE)$X.new-
          SS1Sce1_MixedLPCM_mu)^2))) # mean pairwise distance between mu_t and mu*
    )
  SS1Sce1_MixedLPCM_out_K4_RandomInit_musdDist_list <- 
    c(SS1Sce1_MixedLPCM_out_K4_RandomInit_musdDist_list,
      sd(sqrt(rowSums((Procrustes(
        t(LabelSwitching_SG2003(z=SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$hat_z,
                                matrix_dbyK=t(SS1Sce1_MixedLPCM_out_K4_RandomInit[[R]]$mu_t))$matrix_dbyK),
        SS1Sce1_MixedLPCM_mu, translate = TRUE ,dilate = FALSE)$X.new-
          SS1Sce1_MixedLPCM_mu)^2))) # mean pairwise distance between mu_t and mu*
    )
}
SS1Sce1_MixedLPCM_out_K4_RandomInit_muMeanDist_list
max(SS1Sce1_MixedLPCM_out_K4_RandomInit_muMeanDist_list[-1],na.rm=TRUE) # 0.3222257
median(SS1Sce1_MixedLPCM_out_K4_RandomInit_muMeanDist_list[-1],na.rm=TRUE) # 0.1406596
min(SS1Sce1_MixedLPCM_out_K4_RandomInit_muMeanDist_list[-1],na.rm=TRUE) # 0.1126088
sd(SS1Sce1_MixedLPCM_out_K4_RandomInit_muMeanDist_list[-1],na.rm=TRUE) # 0.06731361

SS1Sce1_MixedLPCM_out_K4_RandomInit_musdDist_list
max(SS1Sce1_MixedLPCM_out_K4_RandomInit_musdDist_list[-1],na.rm=TRUE) # 0.1367614
median(SS1Sce1_MixedLPCM_out_K4_RandomInit_musdDist_list[-1],na.rm=TRUE) # 0.07822845
min(SS1Sce1_MixedLPCM_out_K4_RandomInit_musdDist_list[-1],na.rm=TRUE) # 0.03587664
sd(SS1Sce1_MixedLPCM_out_K4_RandomInit_musdDist_list[-1],na.rm=TRUE) # 0.03094033

# Check E(sigma2_t) and sd(sigma2_t)
SS1Sce1_MixedLPCM_out_K4_RandomInit_Esigma2_t_list
max(SS1Sce1_MixedLPCM_out_K4_RandomInit_Esigma2_t_list[-1]) # 0.001248794
median(SS1Sce1_MixedLPCM_out_K4_RandomInit_Esigma2_t_list[-1]) # 0.001007179
min(SS1Sce1_MixedLPCM_out_K4_RandomInit_Esigma2_t_list[-1]) # 0.0009881262
sd(SS1Sce1_MixedLPCM_out_K4_RandomInit_Esigma2_t_list[-1]) # 8.219956e-05

SS1Sce1_MixedLPCM_out_K4_RandomInit_sdsigma2_t_list
max(SS1Sce1_MixedLPCM_out_K4_RandomInit_sdsigma2_t_list[-1]) # 0.0005377791
median(SS1Sce1_MixedLPCM_out_K4_RandomInit_sdsigma2_t_list[-1]) # 0.0001143457
min(SS1Sce1_MixedLPCM_out_K4_RandomInit_sdsigma2_t_list[-1]) # 7.914545e-05
sd(SS1Sce1_MixedLPCM_out_K4_RandomInit_sdsigma2_t_list[-1]) # 0.0001458901

# Check E(varphi2_t) and sd(varphi2_t)
SS1Sce1_MixedLPCM_out_K4_RandomInit_Evarphi2_t_list
max(SS1Sce1_MixedLPCM_out_K4_RandomInit_Evarphi2_t_list[-1]) # 0.1276918
median(SS1Sce1_MixedLPCM_out_K4_RandomInit_Evarphi2_t_list[-1]) # 0.1012873
min(SS1Sce1_MixedLPCM_out_K4_RandomInit_Evarphi2_t_list[-1]) # 0.09964086
sd(SS1Sce1_MixedLPCM_out_K4_RandomInit_Evarphi2_t_list[-1]) # 0.008824039

SS1Sce1_MixedLPCM_out_K4_RandomInit_sdvarphi2_t_list
max(SS1Sce1_MixedLPCM_out_K4_RandomInit_sdvarphi2_t_list[-1]) # 0.07087906
median(SS1Sce1_MixedLPCM_out_K4_RandomInit_sdvarphi2_t_list[-1]) # 0.01447806
min(SS1Sce1_MixedLPCM_out_K4_RandomInit_sdvarphi2_t_list[-1]) # 0.01257286
sd(SS1Sce1_MixedLPCM_out_K4_RandomInit_sdvarphi2_t_list[-1]) # 0.01886141

# Check mean and sd absolute error of pairwise distance squares
SS1Sce1_MixedLPCM_out_K4_RandomInit_MAEdistsquare_list
max(SS1Sce1_MixedLPCM_out_K4_RandomInit_MAEdistsquare_list[-1]) # 2.886559
median(SS1Sce1_MixedLPCM_out_K4_RandomInit_MAEdistsquare_list[-1]) # 1.598545
min(SS1Sce1_MixedLPCM_out_K4_RandomInit_MAEdistsquare_list[-1]) # 1.519087
sd(SS1Sce1_MixedLPCM_out_K4_RandomInit_MAEdistsquare_list[-1]) # 0.458675

SS1Sce1_MixedLPCM_out_K4_RandomInit_sdAEdistsquare_list
max(SS1Sce1_MixedLPCM_out_K4_RandomInit_sdAEdistsquare_list[-1]) # 2.986887
median(SS1Sce1_MixedLPCM_out_K4_RandomInit_sdAEdistsquare_list[-1]) # 1.65491
min(SS1Sce1_MixedLPCM_out_K4_RandomInit_sdAEdistsquare_list[-1]) # 1.553111
sd(SS1Sce1_MixedLPCM_out_K4_RandomInit_sdAEdistsquare_list[-1]) # 0.4734222

# Check mean and sd of a_t/b_t
SS1Sce1_MixedLPCM_out_K4_RandomInit_E_a_t_div_b_t_list
max(SS1Sce1_MixedLPCM_out_K4_RandomInit_E_a_t_div_b_t_list[-1]) # 9.826578
median(SS1Sce1_MixedLPCM_out_K4_RandomInit_E_a_t_div_b_t_list[-1]) # 9.69322
min(SS1Sce1_MixedLPCM_out_K4_RandomInit_E_a_t_div_b_t_list[-1]) # 8.292436
sd(SS1Sce1_MixedLPCM_out_K4_RandomInit_E_a_t_div_b_t_list[-1]) # 0.4556187

SS1Sce1_MixedLPCM_out_K4_RandomInit_sd_a_t_div_b_t_list
max(SS1Sce1_MixedLPCM_out_K4_RandomInit_sd_a_t_div_b_t_list[-1]) # 1.927929
median(SS1Sce1_MixedLPCM_out_K4_RandomInit_sd_a_t_div_b_t_list[-1]) # 0.8806734
min(SS1Sce1_MixedLPCM_out_K4_RandomInit_sd_a_t_div_b_t_list[-1]) # 0.7603967
sd(SS1Sce1_MixedLPCM_out_K4_RandomInit_sd_a_t_div_b_t_list[-1]) # 0.3627962

# Check mean and sd of b_t/a_t
SS1Sce1_MixedLPCM_out_K4_RandomInit_E_b_t_div_a_t_list
max(SS1Sce1_MixedLPCM_out_K4_RandomInit_E_b_t_div_a_t_list[-1]) # 0.1346866
median(SS1Sce1_MixedLPCM_out_K4_RandomInit_E_b_t_div_a_t_list[-1]) # 0.1042271
min(SS1Sce1_MixedLPCM_out_K4_RandomInit_E_b_t_div_a_t_list[-1]) # 0.1024162
sd(SS1Sce1_MixedLPCM_out_K4_RandomInit_E_b_t_div_a_t_list[-1]) # 0.0101824

SS1Sce1_MixedLPCM_out_K4_RandomInit_sd_b_t_div_a_t_list
max(SS1Sce1_MixedLPCM_out_K4_RandomInit_sd_b_t_div_a_t_list[-1]) # 0.07533027
median(SS1Sce1_MixedLPCM_out_K4_RandomInit_sd_b_t_div_a_t_list[-1]) # 0.01131912
min(SS1Sce1_MixedLPCM_out_K4_RandomInit_sd_b_t_div_a_t_list[-1]) # 0.008581765
sd(SS1Sce1_MixedLPCM_out_K4_RandomInit_sd_b_t_div_a_t_list[-1]) # 0.0215327
```

### 2.2 Simulation Study 1 Scenario 2

Here, similar to SS1 scenario 1, we begin with the simulation process:

``` r
# Simulation study 1 scenario 2

# Simulate from the PoisLPCM and store the data

# SS1Sce2_PoisLPCM <- Simulation_Directed_PoisLPCM(beta=1,mu=rbind(c(1.25,1.25),c(1.25,-1.25),c(-1.25,-1.25),c(-1.25,1.25)),tau=c(8,6,4,2),
#                                                      d=2,z=c(rep(1,25),rep(2,25),rep(3,25),rep(4,25)),seed=NULL)
# write.csv(SS1Sce2_PoisLPCM$Y,"SS1Sce2_PoisLPCM_Y.csv", row.names = FALSE)
# write.csv(SS1Sce2_PoisLPCM$z,"SS1Sce2_PoisLPCM_z.csv", row.names = FALSE)
# write.csv(SS1Sce2_PoisLPCM$U,"SS1Sce2_PoisLPCM_U.csv", row.names = FALSE)

# Load the data

SS1Sce2_MixedLPCM <- list(Y = as.matrix(read.csv("SS1Sce2_MixedLPCM_Y.csv",header = TRUE)),
                          z = c(as.matrix(read.csv("SS1Sce2_MixedLPCM_z.csv",header = TRUE))),
                          U = as.matrix(read.csv("SS1Sce2_MixedLPCM_U.csv",header = TRUE)))
SS1Sce2_MixedLPCM$Z <- t(t(matrix(SS1Sce2_MixedLPCM$z,length(SS1Sce2_MixedLPCM$z),max(SS1Sce2_MixedLPCM$z)))==(1:max(SS1Sce2_MixedLPCM$z)))*1
SS1Sce2_MixedLPCM$V <- array(as.matrix(read.csv("SS1Sce2_MixedLPCM_V.csv",header = TRUE)),
                             dim=c(dim(t(SS1Sce2_MixedLPCM$U)),nrow(SS1Sce2_MixedLPCM$Y)))
colnames(SS1Sce2_MixedLPCM$Y) <- colnames(SS1Sce2_MixedLPCM$U) <- colnames(SS1Sce2_MixedLPCM$V) <- NULL

SS1Sce2_MixedLPCM_beta <- 1
SS1Sce2_MixedLPCM_mu <- rbind(c(1.25,1.25),c(1.25,-1.25),c(-1.25,-1.25),c(-1.25,1.25))
SS1Sce2_MixedLPCM_tau <- c(8,6,4,2)
SS1Sce2_MixedLPCM_gamma <- rep(10,100)
```

The model fitting code and summarizing code is very similar to that shown in the previous SS1 Sce1 sections, so we refer the readers to previous sections for more details. In pratice, simply replacing all "SS1Sce1" by "SS1Sce2" should work.

## 3. Simulation Study 2

Similar as before, we begin with the simulation process:

``` r
# Simulation study 2

# Simulate from the MLPCM and store the data

# SS2_MixedLPCM <- Simulation_Directed_MixedLPCM(beta=1,mu=rbind(c(1.5,1.5),c(1.5,-1.5),c(-1.5,-1),c(-1.5,1)),tau=c(8,6,4,2),
#                                                        gamma=rep(c(rep(1,5),rep(10,20)),4),d=2,
#                                                        z=c(rep(1,25),rep(2,25),rep(3,25),rep(4,25)),seed=NULL)

# write.csv(SS2_MixedLPCM$Y,"SS2_MixedLPCM_Y.csv", row.names = FALSE)
# write.csv(SS2_MixedLPCM$z,"SS2_MixedLPCM_z.csv", row.names = FALSE)
# write.csv(SS2_MixedLPCM$U,"SS2_MixedLPCM_U.csv", row.names = FALSE)
# write.csv(SS2_MixedLPCM$V,"SS2_MixedLPCM_V.csv", row.names = FALSE)

# Load the data

SS2_MixedLPCM <- list(Y = as.matrix(read.csv("SS2_MixedLPCM_Y.csv",header = TRUE)),
                              z = c(as.matrix(read.csv("SS2_MixedLPCM_z.csv",header = TRUE))),
                              U = as.matrix(read.csv("SS2_MixedLPCM_U.csv",header = TRUE)))
SS2_MixedLPCM$Z <- t(t(matrix(SS2_MixedLPCM$z,length(SS2_MixedLPCM$z),max(SS2_MixedLPCM$z)))==(1:max(SS2_MixedLPCM$z)))*1
SS2_MixedLPCM$V <- array(as.matrix(read.csv("SS2_MixedLPCM_V.csv",header = TRUE)),
                                 dim=c(dim(t(SS2_MixedLPCM$U)),nrow(SS2_MixedLPCM$Y)))
colnames(SS2_MixedLPCM$Y) <- colnames(SS2_MixedLPCM$U) <- colnames(SS2_MixedLPCM$V) <- NULL

SS2_MixedLPCM_beta <- 1
SS2_MixedLPCM_mu <- rbind(c(1.5,1.5),c(1.5,-1.5),c(-1.5,-1),c(-1.5,1))
SS2_MixedLPCM_tau <- c(8,6,4,2)
SS2_MixedLPCM_gamma <- rep(c(rep(1,5),rep(10,20)),4)
```

Here, for the model fitting, we provide an example code for $K=4$ case here.

``` r
# Apply MixedLPCM VB on the simulated MixedLPCM K = 4 with default prior settings
time_start <- Sys.time()
SS2_MixedLPCM_out_K4 <- VB_Directed_MixedLPCM(Y=SS2_MixedLPCM$Y,K=4,a=10,b=1,eta=1,rho2=1,omega2=1,xi=1,psi=1,delta=1,tol=1e-2)
time_end <- Sys.time()
SS2_MixedLPCM_out_K4$time <- time_end-time_start
SS2_MixedLPCM_out_K4$hat_z <- apply(SS2_MixedLPCM_out_K4$Pi_t,1,which.max)
SS2_MixedLPCM_out_K4$LS_hat_z <- LabelSwitching_SG2003(z=SS2_MixedLPCM_out_K4$hat_z)$z
# Store the PICL criteria value
SS2_MixedLPCM_out_K4$PICL <-
  PICL_Directed_MixedLPCM(Y=SS2_MixedLPCM$Y,
                          hat_U=SS2_MixedLPCM_out_K4$U_t,
                          hat_V=SS2_MixedLPCM_out_K4$V_t,
                          hat_z=SS2_MixedLPCM_out_K4$LS_hat_z,
                          K=4,a=10,b=1,omega2=1,delta=1,tol=.Machine$double.xmin)

#---------------------------------------------------------------------------------------------------------------
# plot the latent positions
g_obs <- graph_from_adjacency_matrix(SS2_MixedLPCM$Y,mode = "directed",weighted=TRUE)
# Node sizes are proportional to their social impact, \tilde{\gamma}_i
V(g_obs)$size <- 75*(1/(SS2_MixedLPCM_out_K4$a_t/SS2_MixedLPCM_out_K4$b_t))/
  max(1/(SS2_MixedLPCM_out_K4$a_t/SS2_MixedLPCM_out_K4$b_t))
# Additional graphical settings
V(g_obs)$frame.color <- c("black","red")[c((SS2_MixedLPCM_gamma==1)+1)]
V(g_obs)$label <- ""
# Node colors indicate the inferred clustering
Customized_colors <- My_colors[6:9]
V(g_obs)$color <- adjustcolor(Customized_colors[SS2_MixedLPCM_out_K4$LS_hat_z], alpha.f = .5)
E(g_obs)$color <- colorRampPalette(brewer.pal(9,"Greys")[c(3,9)])(max(E(g_obs)$weight))[E(g_obs)$weight]
# Make the position plot
SS2_MixedLPCM_out_K4$Pi_t_pie <- list()
for (i in 1:nrow(SS2_MixedLPCM_out_K4$Pi_t)){
  SS2_MixedLPCM_out_K4$Pi_t_pie[[i]] <- SS2_MixedLPCM_out_K4$Pi_t[i,]
}
par(mfrow=c(1,1),mai = c(0.3, 0.3, 0.3, 0.05),mgp = c(1,0.45,0))
plot(g_obs,rescale=F, layout=Procrustes(SS2_MixedLPCM_out_K4$U_t,SS2_MixedLPCM$U, translate = TRUE ,dilate = FALSE)$X.new,
     edge.curved=0.2,vertex.shape="pie",
     vertex.pie=SS2_MixedLPCM_out_K4$Pi_t_pie,
     vertex.pie.color=list(adjustcolor(Customized_colors[c(1,3,4,2)], alpha.f = .5)),#mark.groups=c(3,7,4,9,14,15), mark.col="#f0f0f0", mark.border=NA
     edge.arrow.size=0.055,xlim=c(-3.50,2.50),ylim=c(-2.75,2.75),edge.width=E(g_obs)$weight*0.25)
axis(1)
axis(2)
box()
par(mfrow=c(1,1),mai = c(1.02, 0.82, 0.82, 0.42))
```

Recall here that we also experiment on a case where we have a better guess in the prior settings, i.e. the adjusted prior case. In our real life, when we get a new network data, we might usually already have some prior information of which person might be more influential and which one might be less. So this experiment aims to mimic this situation to show that better prior information does improve the performance.

``` r
# Apply MixedLPCM VB on the simulated MixedLPCM K = 4 with Adjusted prior settings

time_start <- Sys.time()
SS2_MixedLPCM_out_K4_AdjustedPrior <- VB_Directed_MixedLPCM_VecPrior_a(Y=SS2_MixedLPCM$Y,K=4,a=SS2_MixedLPCM_gamma,b=1,eta=1,rho2=1,omega2=3,xi=1,psi=1,delta=1,tol=1e-2)
time_end <- Sys.time()
SS2_MixedLPCM_out_K4_AdjustedPrior$time <- time_end-time_start
SS2_MixedLPCM_out_K4_AdjustedPrior$hat_z <- apply(SS2_MixedLPCM_out_K4_AdjustedPrior$Pi_t,1,which.max)
SS2_MixedLPCM_out_K4_AdjustedPrior$LS_hat_z <- LabelSwitching_SG2003(z=SS2_MixedLPCM_out_K4_AdjustedPrior$hat_z)$z
# Store the PICL criteria value
SS2_MixedLPCM_out_K4_AdjustedPrior$PICL <-
  PICL_Directed_MixedLPCM(Y=SS2_MixedLPCM$Y,
                          hat_U=SS2_MixedLPCM_out_K4_AdjustedPrior$U_t,
                          hat_V=SS2_MixedLPCM_out_K4_AdjustedPrior$V_t,
                          hat_z=SS2_MixedLPCM_out_K4_AdjustedPrior$LS_hat_z,
                          K=4,a=SS2_MixedLPCM_gamma,b=1,omega2=3,delta=1,tol=.Machine$double.xmin)
```

## 4. Real Data Application

In this real data application, we focus on the International Arms Transfers Database, which
is publicly available from the Stockholm International Peace Research Institute (SIPRI).
The database includes all transfers of major conventional arms between countries, regions,
subregions or non-state actors from 1950 to the most recent full calendar year. We focus
on a subset of this dataset, corresponding to the ego-centric network of the United States
for 2024. 

Without loss of generality, we define that a directed weighted edge exists
from country $i$ to country $j$ if country $i$ ordered arms transfers from country $j$, and
the interaction weight corresponds to the volume of such arms transfers in 2024 that is
recorded in millions of SIPRI Trend-Indicator Values (TIVs). In other words, the edge
direction follows the direction of the order, thus it is opposite to how the arms are moving.
The SIPRI TIV is a novel common unit developed by SIPRI aiming to regularize the data
of deliveries of different weapons and to identify general trends of these arms transfers.

We begin by listing all the countries in 2024 based on the most recent SIPRI Military Expenditure Database. There are 163 countries which have records of military expenditures from 2001 to 2023, and, among these 163 countries, there are 110 countries which have records of arms transfers in 2024. To obtain the egocentric network for the United States, we extract the interaction data in which the United States are involved. After
this step, we remove the United States and all the countries that only have arms transfers with the United States. The remaining network consists of all the interactions that the neighbors of the United States have between themselves, excluding the US themselves and any isolated node. This effectively corresponds to the network as it is seen from the point of view of the United States.

The code for obtaining the resulting network is shown below:

``` r
# # Extract the latest Military expenditures to choose the countries
# 
Military_Spending <- read.csv("military-spending-sipri.csv",header = TRUE) # Extract the original data
Military_Spending <- Military_Spending[Military_Spending[,2] != "",] # keep the country-level data and remove the region-level data
Military_Spending_countries <- levels(as.factor(Military_Spending[,1]))
Military_Spending_Latest <- data.frame()
for (i in 1:length(Military_Spending_countries)){
  Military_Spending_Country_i <- Military_Spending[Military_Spending[,1]==Military_Spending_countries[i],]
  Military_Spending_Latest <- rbind(Military_Spending_Latest,
                                    Military_Spending_Country_i[nrow(Military_Spending_Country_i),])
}
Military_Spending_Latest
Military_Spending_Latest <- Military_Spending_Latest[Military_Spending_Latest[,3]>2000,] # remove countries of which the data only available before year 2000
Military_Spending_Latest <- Military_Spending_Latest[Military_Spending_Latest[,1]!="World",] # remove the world-level data
nrow(Military_Spending_Latest)
order(Military_Spending_Latest[,4],decreasing = TRUE) # rank the Military_Spending in descending order
Military_Spending_Latest$Ranking <- rank(-Military_Spending_Latest[,4], ties.method= "min")
Military_Spending_Latest[order(Military_Spending_Latest[,4],decreasing = TRUE),]
Military_Spending_Latest_Ranking <- Military_Spending_Latest[order(Military_Spending_Latest[,4],decreasing = TRUE),] # store the ranking
rownames(Military_Spending_Latest_Ranking) <- NULL
Military_Spending_Latest_Ranking


#--------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------

# Arms Transfers between 163 countries 2024

# Read the entire file as plain text lines
Arms_transfer_lines <- readLines("163_Countries_TIV_2024.csv")

# Find lines where the split lines start
Arms_transfer_split_points <- grep("^Volume of transfers of major arms", Arms_transfer_lines)
# Add end boundaries
Arms_transfer_split_points <- c(Arms_transfer_split_points, length(Arms_transfer_lines) + 1)

# Extract tables
Arms_transfer_tables <- lapply(seq_along(Arms_transfer_split_points)[-length(Arms_transfer_split_points)], function(i) {
  start <- Arms_transfer_split_points[i] + 9
  end <- Arms_transfer_split_points[i + 1] - 1  # skip the metadata block
  read.csv(text = paste(Arms_transfer_lines[start:end], collapse = "\n"))
})

# Extract all countries which have exports/imports
Arms_transfer_countries <- c()
for (i in 1:length(Arms_transfer_tables)){
  Arms_transfer_countries <- c(Arms_transfer_countries,
                               Arms_transfer_tables[[i]][,1][-length(Arms_transfer_tables[[i]][,1])], # extract import countries
                               sub(".*from (.*)", "\\1", Arms_transfer_tables[[i]][,1][length(Arms_transfer_tables[[i]][,1])])) # extract the export country
}
Arms_transfer_countries <- levels(as.factor(Arms_transfer_countries))

# Construct the adj matrix for the Arms transfers of countries which have imports/exports
Arms_transfer_adj <- Arms_transfer_log2adj <- matrix(0,length(Arms_transfer_countries),length(Arms_transfer_countries))
# res1 <- c()
for (i in 1:length(Arms_transfer_tables)){
  # Find the row that the exporting country at
  column_indicator <- which.max(Arms_transfer_countries==sub(".*from (.*)", "\\1", Arms_transfer_tables[[i]][,1][length(Arms_transfer_tables[[i]][,1])]))
  # res1 <- c(res1,nrow(Arms_transfer_tables[[i]])-1)
  if (length(column_indicator)>0){ # skip the empty tables
    for (j in 1:(nrow(Arms_transfer_tables[[i]])-1)){
      # Find the column that the importing country at
      row_indicator <- which.max(Arms_transfer_countries==Arms_transfer_tables[[i]][j,1])
      # Assign the TIV value to the entry
      if (is.na(Arms_transfer_tables[[i]][j,2])){ # if NA, assign 0
        Arms_transfer_adj[row_indicator,column_indicator] <- Arms_transfer_log2adj[row_indicator,column_indicator] <- 0
      }else{
        if (Arms_transfer_tables[[i]][j,2]!=0){
          Arms_transfer_adj[row_indicator,column_indicator] <- Arms_transfer_tables[[i]][j,2]
          Arms_transfer_log2adj[row_indicator,column_indicator] <- round(log2(Arms_transfer_tables[[i]][j,2]+1)) # log adj: round(log2(x+1))
        }else if (Arms_transfer_tables[[i]][j,2]==0){
          Arms_transfer_adj[row_indicator,column_indicator] <- 1 # 0 values correspond to the TIV between 0 to 0.5; so we assign 1 to it to distinguish with no import/export
          Arms_transfer_log2adj[row_indicator,column_indicator] <- round(log2(2)) # we assign round(log2(2)) in log adj
        }
      }
    }
  }
}

# Countries ranking based on exports+imports
length(Arms_transfer_countries) # check number of countries have exports/imports in 2024
(rowSums(Arms_transfer_adj)+colSums(Arms_transfer_adj))/(2*sum(Arms_transfer_adj)) # check exports+imports importance score for each country
order((rowSums(Arms_transfer_adj)+colSums(Arms_transfer_adj))/(2*sum(Arms_transfer_adj)),decreasing = TRUE) # show positions of rankings
sort((rowSums(Arms_transfer_adj)+colSums(Arms_transfer_adj))/(2*sum(Arms_transfer_adj)),decreasing = TRUE) # rank the scores
Arms_transfer_countries[order((rowSums(Arms_transfer_adj)+colSums(Arms_transfer_adj))/(2*sum(Arms_transfer_adj)),decreasing = TRUE)] # corresponding ranked countries
cumsum(sort((rowSums(Arms_transfer_adj)+colSums(Arms_transfer_adj))/(2*sum(Arms_transfer_adj)),decreasing = TRUE)) # cumulative sum of the ranked scores

# # -----------------------------------------------------------------------------------------------------------
# Only keep those countries which have interactions with the US
which(Arms_transfer_countries=="United States")
Arms_transfer_countries_with_US_Index <- ((Arms_transfer_log2adj[107,]>0)+(Arms_transfer_log2adj[,107]>0)>0)
Arms_transfer_countries_with_US_Index[107] <- TRUE
Arms_transfer_removed_countries <- Arms_transfer_countries[Arms_transfer_countries_with_US_Index==FALSE]
Arms_transfer_countries <- Arms_transfer_countries[Arms_transfer_countries_with_US_Index]
Arms_transfer_log2adj <- Arms_transfer_log2adj[Arms_transfer_countries_with_US_Index,Arms_transfer_countries_with_US_Index]
Arms_transfer_adj <- Arms_transfer_adj[Arms_transfer_countries_with_US_Index,Arms_transfer_countries_with_US_Index]

# remove the United States
Arms_transfer_adj <- Arms_transfer_adj[-which(Arms_transfer_countries=="United States"),-which(Arms_transfer_countries=="United States")]
Arms_transfer_log2adj <- Arms_transfer_log2adj[-which(Arms_transfer_countries=="United States"),-which(Arms_transfer_countries=="United States")]
Arms_transfer_countries <- Arms_transfer_countries[-which(Arms_transfer_countries=="United States")]
Arms_transfer_removed_countries <- sort(c(Arms_transfer_removed_countries,"United States"))
# remove other countries which no longer has arms transports
Arms_transfer_removed_countries_after_US_Index <- which((rowSums(Arms_transfer_adj)+colSums(Arms_transfer_adj))==0)
Arms_transfer_adj <- Arms_transfer_adj[-Arms_transfer_removed_countries_after_US_Index,-Arms_transfer_removed_countries_after_US_Index]
Arms_transfer_log2adj <- Arms_transfer_log2adj[-Arms_transfer_removed_countries_after_US_Index,-Arms_transfer_removed_countries_after_US_Index]
Arms_transfer_removed_countries <- sort(c(Arms_transfer_removed_countries,Arms_transfer_countries[Arms_transfer_removed_countries_after_US_Index]))
Arms_transfer_countries <- Arms_transfer_countries[-Arms_transfer_removed_countries_after_US_Index]

# Plot the adj
Arms_transfer_adj_dataframe <- as.data.frame(Arms_transfer_adj)
rownames(Arms_transfer_adj_dataframe) <- colnames(Arms_transfer_adj_dataframe) <- 1:nrow(Arms_transfer_adj)
library("pheatmap")
library("RColorBrewer")
Arms_transfer_adj_heatmap <-
  pheatmap(Arms_transfer_adj_dataframe,
           color=colorRampPalette(brewer.pal(9,"YlOrRd"))(max(Arms_transfer_adj)+1),
           cluster_cols = FALSE,cluster_rows= FALSE,show_rownames=FALSE,show_colnames=FALSE,border_color=FALSE,legend=TRUE)
Arms_transfer_adj_heatmap_draw <- cowplot::ggdraw(Arms_transfer_adj_heatmap[[4]])#+ theme(plot.background =element_rect(fill=brewer.pal(9,"Greys")[3]))
print(Arms_transfer_adj_heatmap_draw)

Arms_transfer_log2adj_dataframe <- as.data.frame(Arms_transfer_log2adj)
rownames(Arms_transfer_log2adj_dataframe) <- colnames(Arms_transfer_log2adj_dataframe) <- 1:nrow(Arms_transfer_log2adj)
library("pheatmap")
library("RColorBrewer")
Arms_transfer_log2adj_heatmap <-
  pheatmap(Arms_transfer_log2adj_dataframe,
           color=colorRampPalette(brewer.pal(9,"YlOrRd"))(max(Arms_transfer_log2adj)+1),
           cluster_cols = FALSE,cluster_rows= FALSE,show_rownames=FALSE,show_colnames=FALSE,border_color=FALSE,legend=TRUE)
Arms_transfer_log2adj_heatmap_draw <- cowplot::ggdraw(Arms_transfer_log2adj_heatmap[[4]])#+ theme(plot.background =element_rect(fill=brewer.pal(9,"Greys")[3]))
print(Arms_transfer_log2adj_heatmap_draw)

library("corrplot")
colnames(Arms_transfer_log2adj) <- rownames(Arms_transfer_log2adj) <- Arms_transfer_countries
# colnames(Arms_transfer_log2adj) <- rownames(Arms_transfer_log2adj) <- 1:nrow(Arms_transfer_log2adj)
corrplot(Arms_transfer_log2adj,col=colorRampPalette(brewer.pal(9,"YlOrRd"))(max(Arms_transfer_log2adj)), method = "square", is.corr = FALSE, addgrid.col = NA, diag = F,
         mar = c(0, 0, 1, 0), tl.cex = 0.3, tl.col = NA, cl.cex = 0.75, cl.align.text = "l", main = "")

#--------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------


# Match the resulting countries with the initial list of countries
# And rename the countries if they have different names in Arms_transfer database and in Military_Spending database

Arms_transfer_countries
Military_Spending_Latest_Ranking
match(Arms_transfer_countries,Military_Spending_Latest_Ranking[,1])
match(c("Democratic Republic of Congo","Turkey","United Arab Emirates","Vietnam"),Military_Spending_Latest_Ranking[,1])
Military_Spending_Latest_Ranking[,1][
  match(c("Democratic Republic of Congo","Turkey","United Arab Emirates","Vietnam"),Military_Spending_Latest_Ranking[,1])] <-
  c("DR Congo","Turkiye","UAE","Viet Nam") # Rename the countries which are the same
Military_Spending_Latest_Ranking
match(Arms_transfer_countries,Military_Spending_Latest_Ranking[,1])
Military_Spending_Latest_Ranking[match(Arms_transfer_countries,Military_Spending_Latest_Ranking[,1]),]
```

The example model fitting code for the K=3 case is provided below:

``` r
# Apply MixedLPCM VB on the directed weighted Arms_transfer dataset K = 3 with the prior settings that encourage separation of different groups
time_start <- Sys.time()
RDA_Arms_transfer_MixedLPCM_out_K3 <- VB_Directed_MixedLPCM(Y=Arms_transfer_log2adj,K=3,a=10,b=1,eta=1,rho2=1,omega2=3,xi=3,psi=1,delta=3,tol=1e-2)
time_end <- Sys.time()
RDA_Arms_transfer_MixedLPCM_out_K3$time <- time_end-time_start
RDA_Arms_transfer_MixedLPCM_out_K3$hat_z <- apply(RDA_Arms_transfer_MixedLPCM_out_K3$Pi_t,1,which.max)
RDA_Arms_transfer_MixedLPCM_out_K3$LS_hat_z <- LabelSwitching_SG2003(z=RDA_Arms_transfer_MixedLPCM_out_K3$hat_z)$z
# Store the PICL criteria value
RDA_Arms_transfer_MixedLPCM_out_K3$PICL <-
  PICL_Directed_MixedLPCM(Y=Arms_transfer_log2adj,
                          hat_U=RDA_Arms_transfer_MixedLPCM_out_K3$U_t,
                          hat_V=RDA_Arms_transfer_MixedLPCM_out_K3$V_t,
                          hat_z=RDA_Arms_transfer_MixedLPCM_out_K3$LS_hat_z,
                          K=3,a=10,b=1,omega2=3,delta=3,tol=.Machine$double.xmin)
RDA_Arms_transfer_MixedLPCM_out_K3$PICL$PICL # 172.9893

RDA_Arms_transfer_MixedLPCM_out_K3$time # Time difference of 9.919651 mins
RDA_Arms_transfer_MixedLPCM_out_K3$NumIt # 6739
table(apply(RDA_Arms_transfer_MixedLPCM_out_K3$Pi_t,1,which.max)) # Check clustering
RDA_Arms_transfer_MixedLPCM_out_K3$Pi_t
RDA_Arms_transfer_MixedLPCM_out_K3$hat_z <- apply(RDA_Arms_transfer_MixedLPCM_out_K3$Pi_t,1,which.max)#LabelSwitching_SG2003(z=apply(RDA_Arms_transfer_MixedLPCM_out_K3$Pi_t,1,which.max))$z
RDA_Arms_transfer_MixedLPCM_out_K3$LS_hat_z <- LabelSwitching_SG2003(z=RDA_Arms_transfer_MixedLPCM_out_K3$hat_z)$z
RDA_Arms_transfer_MixedLPCM_out_K3$ELBO_list[length(RDA_Arms_transfer_MixedLPCM_out_K3$ELBO_list)]
# Check estimated LPM intercept eta_t for the model parameter beta
RDA_Arms_transfer_MixedLPCM_out_K3$eta_t # 1.32703
RDA_Arms_transfer_MixedLPCM_out_K3$rho2_t # 0.001276286
# Check estimated Dirichlet parameters, delta_t
RDA_Arms_transfer_MixedLPCM_out_K3$delta_t
# Check estimated covert latent positions' group precision, tau_t
RDA_Arms_transfer_MixedLPCM_out_K3$xi_t
RDA_Arms_transfer_MixedLPCM_out_K3$psi_t
RDA_Arms_transfer_MixedLPCM_out_K3$xi_t/RDA_Arms_transfer_MixedLPCM_out_K3$psi_t
RDA_Arms_transfer_MixedLPCM_out_K3$xi_t/(RDA_Arms_transfer_MixedLPCM_out_K3$psi_t^2)
RDA_Arms_transfer_MixedLPCM_out_K3$xi_t
RDA_Arms_transfer_MixedLPCM_out_K3$psi_t
(RDA_Arms_transfer_MixedLPCM_out_K3$xi_t/RDA_Arms_transfer_MixedLPCM_out_K3$psi_t)
(RDA_Arms_transfer_MixedLPCM_out_K3$xi_t/(RDA_Arms_transfer_MixedLPCM_out_K3$psi_t^2))
# Check estimated covert latent positions U_t and group centers mu_t
library("latex2exp")
plot(RDA_Arms_transfer_MixedLPCM_out_K3$U_t,xlim=c(-10,10),ylim=c(-10,10),col=My_colors[RDA_Arms_transfer_MixedLPCM_out_K3$LS_hat_z+5],xlab="",ylab="", main=TeX(r"(R1 $\tilde{U}$)",bold=TRUE))
points(RDA_Arms_transfer_MixedLPCM_out_K3$mu_t,col=My_colors[10],pch=2)
# Check variational MVN variance sigma2_t for U_t
hist(RDA_Arms_transfer_MixedLPCM_out_K3$sigma2_t,breaks=100,main=TeX(r"(R1 $\tilde{\sigma}^2$)",bold=TRUE),xlab="")
# Check variational MVN variance omega2_t for mu_t
RDA_Arms_transfer_MixedLPCM_out_K3$omega2_t

# Check estimated overt latent positions' precision, gamma
RDA_Arms_transfer_MixedLPCM_out_K3$a_t
RDA_Arms_transfer_MixedLPCM_out_K3$b_t
RDA_Arms_transfer_MixedLPCM_out_K3$a_t/RDA_Arms_transfer_MixedLPCM_out_K3$b_t
RDA_Arms_transfer_MixedLPCM_out_K3$a_t/(RDA_Arms_transfer_MixedLPCM_out_K3$b_t^2)

# Check the estimated overt latent positions
RDA_Arms_transfer_MixedLPCM_out_K3_V_t_precision <- c() # obtain the precision (1/var) of the {v_t_ij:i=1,2,\dots,N; i != j} for each covert latent position u_j
for (i in 1:nrow(Arms_transfer_adj)){
  RDA_Arms_transfer_MixedLPCM_out_K3_V_t_precision <- c(RDA_Arms_transfer_MixedLPCM_out_K3_V_t_precision,1/mean(apply(RDA_Arms_transfer_MixedLPCM_out_K3$V_t[,-i,i],1,var)))
}
hist(1/RDA_Arms_transfer_MixedLPCM_out_K3_V_t_precision,breaks=100,main=TeX(r"(R1 $\{\tilde{v}_{j\leftarrow i}:i;i\neq j\}_j$ var)",bold=TRUE),xlab="")
# Check the variational MVN var varphi2_t for the variational overt latent positions
hist(RDA_Arms_transfer_MixedLPCM_out_K3$varphi2_t,breaks=100,main=TeX(r"(R1 $\tilde{\varphi}^2$)",bold=TRUE),xlab="")
mean(RDA_Arms_transfer_MixedLPCM_out_K3$varphi2_t[!is.na(RDA_Arms_transfer_MixedLPCM_out_K3$varphi2_t)]) # 0.1948025
# 
#--------------------------------------------------------------------------------------------------------------------------------------------
# 2024 Egocentric ArmsTransfers K=3 plots
library("igraph")
# plot the inferred latent positions with inferred clustering
g_obs <- graph_from_adjacency_matrix(Arms_transfer_log2adj,mode = "directed",weighted=TRUE)
# Node sizes are proportional to their social impact, \tilde{\gamma}_i
V(g_obs)$size <- 1/(RDA_Arms_transfer_MixedLPCM_out_K3$a_t/RDA_Arms_transfer_MixedLPCM_out_K3$b_t)
V(g_obs)$size <- 125*V(g_obs)$size
# Additional graphical settings
V(g_obs)$frame.color <- "black"
V(g_obs)$label <- ""
# Node colors indicate the inferred clustering
Customized_colors <- My_colors[6:17][c(3,1,2)]
V(g_obs)$color <- adjustcolor(Customized_colors[RDA_Arms_transfer_MixedLPCM_out_K3$hat_z], alpha.f = .5)
# Edge colors proportional to edge weights
E(g_obs)$color <- colorRampPalette(brewer.pal(9,"Greys")[c(3,9)])(max(E(g_obs)$weight))[E(g_obs)$weight]

RDA_Arms_transfer_MixedLPCM_out_K3$Pi_t_pie <- list()
for (i in 1:nrow(RDA_Arms_transfer_MixedLPCM_out_K3$Pi_t)){
  RDA_Arms_transfer_MixedLPCM_out_K3$Pi_t_pie[[i]] <- RDA_Arms_transfer_MixedLPCM_out_K3$Pi_t[i,]#[2:1]
}
par(mfrow=c(2,1),mai = c(0.3, 0.3, 0.3, 0.05),mgp = c(1,0.45,0))
plot(g_obs,rescale=F, layout=RDA_Arms_transfer_MixedLPCM_out_K3$U_t,edge.curved=0.1,
     vertex.shape="pie",vertex.pie=RDA_Arms_transfer_MixedLPCM_out_K3$Pi_t_pie,
     main=TeX(r"(2024 Egocentric ArmsTransfers K=3 with ${\tilde{\Pi}}$)",bold=TRUE),
     vertex.pie.color=list(adjustcolor(Customized_colors, alpha.f = .5)),
     edge.arrow.size=0.25,xlim=c(-4.25,4.25),ylim=c(-3.5,4),edge.width=E(g_obs)$weight*0.15)
axis(1)
axis(2)
box()

V(g_obs)$label <- Arms_transfer_countries
V(g_obs)$label.cex <- 0.75
plot(g_obs,rescale=F, layout=RDA_Arms_transfer_MixedLPCM_out_K3$U_t,edge.curved=0.1,
     main=TeX(r"(2024 Egocentric ArmsTransfers K=3 with Point Estimate Clustering)",bold=TRUE),
     edge.arrow.size=0.25,xlim=c(-4.25,4.25),ylim=c(-3.5,4),edge.width=E(g_obs)$weight*0.15)
axis(1)
axis(2)
box()

par(mfrow=c(1,1),mai = c(1.02, 0.82, 0.82, 0.42),mgp = c(3,1,0))
```

And here ends this tutorial for the Mixed-LPCM paper experiments.
