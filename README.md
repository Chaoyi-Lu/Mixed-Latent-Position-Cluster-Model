# Mixed-Latent-Position-Cluster-Model

This repository provides all the R code as well as the data for the simulation studies and ral data application of the paper [**Mixed Latent Position Cluster Model**].

This paper develops a novel mixed latent position cluster model especially for directed social networks to characterize not only weighted interactions, clustering, latent space visualizations, but also possibly significant discrepancy of outgoing and incoming interactions between different individuals. 
The MLPCM extends the mixed-membership idea to latent positions leading to the so-called mixed latent positions, and can also be treated as using a novel approach to incorporate random effects to the LPCM.

The inference is based on the proposed variational Bayes approaches within which the ELBO can be accurately derived.
A novel partially integrated classification log-likelihood model selection criterion is also developed and accommodated to our model in order to find the optimal number of clusters.

The file [`Application.md`] is the main tutorial notes for the simulation studies and the real data applications.
The file [`Functions.R`] contains all the [`R`] functions we need to use.
The [`.csv`] files contain the simulated data and the real data that we implement on for the simulation studies and the real data application.
