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

#### 2.1.1 Fitting the MLPCM in Simulation Study 1 Scenario 1

#### 2.1.2 Fitting the PoisLPCM in Simulation Study 1 Scenario 1

#### 2.1.3 Summarizing the Output in Simulation Study 1 Scenario 1

#### 2.1.4 Multiple implementations with Random initializations in Simulation Study 1 Scenario 1

### 2.2 Simulation Study 1 Scenario 2

## 3. Simulation Study 2

## 4. Real Data Application

