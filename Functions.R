#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------

################ MixedLPCM ################

#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------

# Simulation from the new MixedLPCM

Simulation_Directed_MixedLPCM <- function(beta,mu,tau,gamma,d=2,z=NULL,Pi=NULL,seed=NULL){
  # beta: LPM intercept parameter
  # mu: K X d matrix with each row denoting a mu_k; MVN center for the latent latent position u
  # tau: 1 X K vector with each entry denoting tau_k; MVN precision for the latent latent position u
  # gamma: 1 X N vector with each entry denoting gamma_i; MVN precision for the latent position v
  # Note that c(1:N) can be treated as either a N X 1 matrix or a 1 X N matrix or a 1 X N vector
  # d: dimension of latent positions; d=2 by default
  # z: 1 X N vector; the specified memberships; Here, z is possible to contain empty clusters
  # Pi: 1 X K vector; Clustering probability to each group
  
  set.seed(seed)
  
  K <- nrow(mu)
  #----------------------------------------------------------------------------------------------------------------------------------------
  if (!is.null(z)) {
    N <- length(z)
    Z <- t(t(matrix(z,N,K))==(1:K))*1 # transform the membership 1 X N vector z to a N X K matrix Z; empty groups included
  }else{
    Z <- t(rmultinom(N,1,Pi)) # the rmultinom(N,1,Pi) generated matrix is a K X N matrix
    z <- c(Z%*%c(1:K))
  }
  #----------------------------------------------------------------------------------------------------------------------------------------
  library(mvtnorm)
  n_k <- table(factor(z,levels=1:K)) # empty groups might exist
  # Generate covert latent positions u_i for each i
  U <- matrix(0,N,d)
  for (k in 1:K){
    if (n_k[[k]]!=0){U[which(z==k),] <- mvtnorm::rmvnorm(n_k[[k]],mu[k,],(1/tau[k])*diag(d))} # mvtnorm::rmvnorm() outputs a n_k X d matrix
  }
  #----------------------------------------------------------------------------------------------------------------------------------------
  # Generate overt latent positions v_{j<--i} for each edge i,j
  V <- array(0,dim=c(d,N,N))
  for (j in 1:N){ # the ijth entry of V is v_{j<--i} corresponding to the latent position of node j for y_ij
    V[,,j] <- t(mvtnorm::rmvnorm(N,U[j,],(1/gamma[j])*diag(d))) # mvtnorm::rmvnorm() outputs a N X d matrix
    V[,j,j] <- rep(NA,d)
  }
  dij2 <- apply((array(t(U),dim=c(d,N,N))-V)^2,2:3,sum)
  diag(dij2) <- 0
  #----------------------------------------------------------------------------------------------------------------------------------------
  # Simulate adj Y
  Y <- matrix(rpois(N^2,exp(beta-dij2)),N,N)
  diag(Y) <- 0
  
  return(list(U=U, V=V, Z=Z, z=z, Y=Y))
}

#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------

# VB inference from the new MixedLPCM

VB_Directed_MixedLPCM <- 
  function(Y,d=2, a=1,b=1,eta=1,rho2=1,omega2=1,xi=1,psi=1,K=10,delta=1, # delta=rep(1e-3,K), # data, dimension and prior parameters
           U_t=NA,sigma2_t=NA,V_t=NA,varphi2_t=NA,Pi_t=NA,eta_t=NA,rho2_t=NA,mu_t=NA,omega2_t=NA,xi_t=NA,psi_t=NA,a_t=NA,b_t=NA,delta_t=NA, # initialize VB parameters
           tol=1e-2,epsilon_init=1,seed=NULL){ # VB algorithm settings
    
    # Y: the N X N observed directed binary network adjacency matrix
    # d: the dimension assumed for latent position space
    # a,b: prior parameters for each \gamma_i; default setting Ga(1,1)
    # eta,rho2: prior parameters for \beta; default setting N(0,1)
    # omega2: prior parameter for each \bm{\mu}_k; default setting MVN_d(0,1)
    # xi,psi: prior parameters for each \tau_k; default setting Ga(1,1)
    # K: an upper bound of the number of non-empty groups; default setting K=20; 
    # extra groups will be made redundant automatically by very small \tilde{\pi_ik}
    # delta: a 1 X K vector or a single value if a symmetric Dirichlet is considered; prior parameters for \bm{\Pi}; 
    # default setting Dirichlet(rep(1e-3,K)); Shrinkage prior required that each entry of delta should be < 1/K (Rousseau and Mengersen, 2011)
    
    # tol: tolerance level; if the ELBO^{(t+1)}-ELBO^{(t)} < tol, then terminate the implementation; default setting tol=1e-5
    # epsilon_init: initial epsilon set for each the gradient ascent step; default setting epsilon=1
    # seed: the seed for random number generator
    
    set.seed(seed)
    N <- nrow(Y)
    Check.Replace.Smallest.Number <- function(x) max(x,.Machine$double.xmin)
    Check.Replace.Smallest.Number <- Vectorize(Check.Replace.Smallest.Number)
    # It's possible to have very small numbers during the inference and continuing the process would make the corresponding model parameters become too small so that the machine cannot represent and a zero instead appears
    # This will bring the logarithm of such a parameter being -Inf. So we limit the smallest number to .Machine$double.xmin
    
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    # Initialization
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    
    # U_t: \tilde{U}, a N X d matrix
    # sigma2_t: \tilde{sigma2}, a 1 X N vector
    # V_t: \tilde{V}, a d X N X N array
    # varphi2_t: \tilde{varphi}, a N X N matrix
    # Pi_t: \tilde{Pi}, a N X K matrix; each row is the 1 X K vector \tilde{\Pi}_i
    # eta_t and rho2_t: \tilde{eta} and \tilde{rho2}, two normal distribution parameters
    # mu_t: \tilde{mu}, a K X d matrix
    # omega2_t: \tilde{omega2}, a 1 X K vector
    # xi_t: \tilde{xi}, a 1 X K vector
    # psi_t: \tilde{psi}, a 1 X K vector
    # a_t: \tilde{a}, a 1 X N vector
    # b_t: \tilde{b}, a 1 X N vector
    # delta_t: \tilde{delta}, a 1 X K vector
    
    # We instead denote ( U_tp,sigma2_tp,V_tp,varphi2_tp,Pi_tp,eta_tp,rho2_tp,mu_tp,omega2_tp,xi_tp,psi_tp,a_tp,b_tp,delta_tp ) as the new update
    
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize covert latent positions U_t by Classical (Metric) Multidimensional Scaling (MDS) on asymmetric distance matrix
    library(ergmito) # for geodesic
    Y_GM <- geodesic(Y,force=TRUE)[[1]] # geodesic matrix (GM) of Y, the length/the number of edges of the shortest paths between each pair of nodes
    Y_GM[is.na(Y_GM)] <- 0 # remove NA as 0
    if (is.na(U_t[1])){U_t <- cmdscale(Y_GM,d)} # obtain the initial latent positions by MDS; cmdscale(,d) returns a N X d matrix
    # Note here that MDS may bring the same position for multiple nodes
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize overt latent positions V_t by setting all { v_{j<--i} }^N_{i=1,i!=j} be the same as the u_j for j = 1,2,...,N
    if (is.na(V_t[1])){
      V_t <- array(0,dim=c(d,N,N))
      for (j in 1:N){
        V_t[,,j] <- matrix(U_t[j,],d,N)
        V_t[,j,j] <- rep(NA,d)
      }
    }
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize the group MVN distribution center mu_t by first applying kmeans() on the covert latent positions to obtain the initial grouping,
    # and then let each \tilde{mu}_k be the mean of the latent positions contained in group k
    z_init <- kmeans(U_t,centers=K,nstart = 1000)$cluster # check 1000 different start and pick the best one
    if (is.na(mu_t[1])){
      mu_t <- matrix(0,K,d)
      for (k in 1:K){mu_t[k,] <- apply(matrix(U_t[z_init==k,],ncol=d),2,mean)}
    }
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize sigma2_t, a_t and b_t
    if (is.na(sigma2_t[1])){sigma2_t <- rep(1,N)}
    if (is.na(a_t[1])){a_t <- rep(a,N)}
    if (is.na(b_t[1])){b_t <- rep(b,N)}
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize varphi2_t
    if (is.na(varphi2_t[2])){ # Check whether varphi2_t is NA or a matrix with diag=NA (so the 2nd element is instead checked)
      varphi2_t <- matrix(1,N,N)
      diag(varphi2_t) <- NA
    }
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize omega2_t, xi_t, psi_t and delta_t
    if (is.na(omega2_t[1])){omega2_t <- rep(omega2,K)}
    if (is.na(xi_t[1])){xi_t <- rep(xi,K)}
    if (is.na(psi_t[1])){psi_t <- rep(psi,K)}
    if (is.na(delta_t[1])){delta_t <- rep(delta,K)}
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize eta_t and rho2_t
    if (is.na(eta_t[1])){eta_t <- eta}
    if (is.na(rho2_t[1])){rho2_t <- rho2}
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize Pi_t by evaluating the analytical solution based on all other initial settings
    library(proxy) # for proxy::dist; pairwise distance between two different set of vectors
    delta_t0 <- sum(delta_t)
    if (is.na(Pi_t)){
      # We first define the A_{jk} which is a N X K matrix
      Ajk <- t(matrix(0.5*d*(digamma(xi_t)-log(psi_t))+digamma(delta_t)-digamma(delta_t0),K,N))-
        t(matrix(0.5*xi_t/psi_t,K,N))*(proxy::dist(U_t, mu_t, method="Euclidean")^2+(matrix(sigma2_t,N,K)+t(matrix(omega2_t,K,N)))*d)
      exp_Ajk <- exp(Ajk)
      Pi_t <- exp_Ajk/matrix(rowSums(exp_Ajk),N,K) # a N X K matrix with each row being \Pi_i
      Pi_t <- apply(Pi_t,1:2,Check.Replace.Smallest.Number)
    }
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize the epsilon for gradient ascent steps
    epsilon_U <- rep(epsilon_init,N)
    epsilon_V <- matrix(epsilon_init,N,N); diag(epsilon_V) <- NA
    epsilon_beta <- epsilon_Pi <- epsilon_init
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Evaluate the ELBO at the initial state
    # Calculate the 1st sum over i,j = 1,2,...,N and i != j
    DistSquare_u_t_i_minus_v_t_ij_NNmatrix <- apply((array(t(U_t),dim=c(d,N,N))-V_t)^2,2:3,sum) # NNmatrix means N X N matrix
    sigma2_t_i_plus_varphi2_t_ij_NNmatrix <- matrix(sigma2_t,N,N)+varphi2_t
    The_1_plus_2var_t_NNmatrix <- 1+2*sigma2_t_i_plus_varphi2_t_ij_NNmatrix
    a_t_j_NNmatrix <- t(matrix(a_t,N,N))
    b_t_j_NNmatrix <- t(matrix(b_t,N,N))
    ELBO_t_Sum_ij <- Y*(eta_t-DistSquare_u_t_i_minus_v_t_ij_NNmatrix-d*sigma2_t_i_plus_varphi2_t_ij_NNmatrix)-
      exp(eta_t+rho2_t/2-DistSquare_u_t_i_minus_v_t_ij_NNmatrix/The_1_plus_2var_t_NNmatrix)/(The_1_plus_2var_t_NNmatrix^(0.5*d))+
      0.5*d*(digamma(a_t_j_NNmatrix)-log(b_t_j_NNmatrix)+log(varphi2_t))-
      0.5*a_t_j_NNmatrix/b_t_j_NNmatrix*(apply((V_t-array(do.call(rbind, replicate(N, t(U_t), simplify=FALSE)),dim=c(d,N,N)))^2,2:3,sum)+d*(t(matrix(sigma2_t,N,N))+varphi2_t))
    diag(ELBO_t_Sum_ij) <- 0
    # Calculate the 2nd sum over i = 1,2,...,N and k = 1,2,...,K
    xi_t_k_NKmatrix <- t(matrix(xi_t,K,N))
    psi_t_k_NKmatrix <- t(matrix(psi_t,K,N))
    ELBO_t_Sum_ik <-
      Pi_t*(0.5*d*(digamma(xi_t_k_NKmatrix)-log(psi_t_k_NKmatrix))-0.5*xi_t_k_NKmatrix/psi_t_k_NKmatrix* # (apply((array(t(U_t),dim=c(d,N,K))-array(do.call(rbind, replicate(N, t(mu_t), simplify=FALSE)),dim=c(d,N,K)))^2,2:3,sum) # array() way to calculate ||u_i-mu_k||^2
              (proxy::dist(U_t, mu_t, method="Euclidean")^2+d*(matrix(sigma2_t,N,K)+t(matrix(omega2_t,K,N))))+
              digamma(t(matrix(delta_t,K,N)))-digamma(delta_t0)-log(Pi_t))
    # Calculate the 3rd sum over k = 1,2,...,K
    ELBO_t_Sum_k <- 
      (delta-delta_t)*(digamma(delta_t)-digamma(delta_t0))+lgamma(delta_t)-0.5/omega2*(rowSums(mu_t^2)+d*omega2_t)+0.5*d*log(omega2_t)+
      (xi-xi_t)*digamma(xi_t)-xi*log(psi_t)-psi*xi_t/psi_t+lgamma(xi_t)+xi_t
    # Calculate the 4th sum over i = 1,2,...,N
    ELBO_t_Sum_i <- (a-a_t)*digamma(a_t)-a*log(b_t)-b*a_t/b_t+lgamma(a_t)+a_t+0.5*d*log(sigma2_t)
    # Calculate the ELBO with the rest single terms
    ELBO_t <- sum(ELBO_t_Sum_ij)+sum(ELBO_t_Sum_ik)+sum(ELBO_t_Sum_k)+sum(ELBO_t_Sum_i)-
      0.5/rho2*((eta_t-eta)^2+rho2_t)+0.5*log(rho2_t)-lgamma(delta_t0)
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize the ELBO list
    ELBO_list <- c(ELBO_t)
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize the update state
    U_tp<-U_t; sigma2_tp<-sigma2_t; V_tp<-V_t; varphi2_tp<-varphi2_t; Pi_tp<-Pi_t; a_tp<-a_t; b_tp<-b_t; eta_tp<-eta_t; rho2_tp<-rho2_t
    mu_tp<-mu_t; omega2_tp<-omega2_t; xi_tp<-xi_t; psi_tp<-psi_t; delta_tp<-delta_t
    
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    # Inference steps
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    t <- 0; STOP <- FALSE
    while (STOP == FALSE){
      t <- t + 1
      if ((t%%100) == 0){cat("Iteration:", t,"\n")}
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      # Step 1: update U_t and sigma2_t by NGA
      for (i in 1:N){
        Increase_U <- FALSE
        # Define the terms which are calculated for more than twice in order to reduce the computation burden
        u_t_i_dNminus1matrix <- matrix(U_t[i,],d,N-1) # a d X N-1 matrix with each column being the \tilde{u}_i
        u_t_i_minus_v_t_ij_dNminus1matrix <- u_t_i_dNminus1matrix-V_t[,i,-i] # matrix(V_t[,i,-i],nrow=d)
        DistSquare_u_t_i_minus_v_t_ij_1Nminus1vector <- colSums(u_t_i_minus_v_t_ij_dNminus1matrix^2)
        The_1_plus_2sigma2_t_i_plus_2varphi2_t_ij_1Nminus1vector <- 1+2*(sigma2_t[i]+varphi2_t[i,-i])
        DistSqaure_over_1_plus_2var_t_1Nminus1vector <- DistSquare_u_t_i_minus_v_t_ij_1Nminus1vector/The_1_plus_2sigma2_t_i_plus_2varphi2_t_ij_1Nminus1vector
        The_exp_term_t_1Nminus1vector <- exp(eta_t+rho2_t/2-DistSqaure_over_1_plus_2var_t_1Nminus1vector)
        The_big_shared_term_t_U_1Nminus1vector <- The_exp_term_t_1Nminus1vector/(The_1_plus_2sigma2_t_i_plus_2varphi2_t_ij_1Nminus1vector^(0.5*d+1))
        v_t_ji_minus_u_t_i_dNminus1matrix <- V_t[,-i,i]-u_t_i_dNminus1matrix
        t_mu_t <- t(mu_t)
        u_t_i_minus_mu_t_k_dKmatrix <- matrix(U_t[i,],d,K)-t_mu_t
        pi_t_ik_times_xi_t_k_over_psi_t_k_1Kvector <- Pi_t[i,]*xi_t/psi_t
        # Calculate the ELBO relevant to the current state U_t[i,] and sigma2_t[i]
        ELBO_t_U <- sum(-Y[i,-i]*(DistSquare_u_t_i_minus_v_t_ij_1Nminus1vector+d*sigma2_t[i])-The_exp_term_t_1Nminus1vector/(The_1_plus_2sigma2_t_i_plus_2varphi2_t_ij_1Nminus1vector^(0.5*d))-
                          0.5*a_t[i]/b_t[i]*colSums(v_t_ji_minus_u_t_i_dNminus1matrix^2))+
          0.5*d*(log(sigma2_t[i])-(N-1)*a_t[i]/b_t[i]*sigma2_t[i]) - 0.5*sum(pi_t_ik_times_xi_t_k_over_psi_t_k_1Kvector*(colSums(u_t_i_minus_mu_t_k_dKmatrix^2)+d*sigma2_t[i]))
        # Set the while loop to check whether the ELBO increases or not
        while (Increase_U == FALSE){
          # Update U_tp
          U_tp[i,] <- U_t[i,]+epsilon_U[i]*sigma2_t[i]*
            (rowSums(t(replicate(d,2*(-Y[i,-i]+The_big_shared_term_t_U_1Nminus1vector)))*u_t_i_minus_v_t_ij_dNminus1matrix + a_t[i]/b_t[i]*v_t_ji_minus_u_t_i_dNminus1matrix)-
               rowSums(u_t_i_minus_mu_t_k_dKmatrix*t(matrix(replicate(d,pi_t_ik_times_xi_t_k_over_psi_t_k_1Kvector),K,d)))) # matrix(,K,d) here aims to deal with the K=1 case
          # Update sigma2_tp
          sigma2_tp[i] <- sigma2_t[i]*exp(epsilon_U[i]*2/d*sigma2_t[i]*
                                            (0.5*d*(1/sigma2_t[i]-(N-1)*a_t[i]/b_t[i]-sum(pi_t_ik_times_xi_t_k_over_psi_t_k_1Kvector))-
                                               sum(d*Y[i,-i]+The_big_shared_term_t_U_1Nminus1vector*(2*DistSqaure_over_1_plus_2var_t_1Nminus1vector-d))))
          if (is.nan(sigma2_tp[i])){sigma2_tp[i] <- sigma2_t[i]} # it's possible that some initial state would first bring very very very small sigma2_t leading to NaN; so in this case we lower-bound the value
          # Determine whether the ELBO increases or not; if no increase, half the epsilon and repeat again
          # Calculate the ELBO relevant to the update state U_tp[i,] and sigma2_tp[i]
          u_tp_i_dNminus1matrix <- matrix(U_tp[i,],d,N-1)
          DistSquare_u_tp_i_minus_v_t_ij_1Nminus1vector <- colSums((u_tp_i_dNminus1matrix-V_t[,i,-i])^2)
          The_1_plus_2sigma2_tp_i_plus_2varphi2_t_ij_1Nminus1vector <- 1+2*(sigma2_tp[i]+varphi2_t[i,-i])
          ELBO_tp_U <- sum(-Y[i,-i]*(DistSquare_u_tp_i_minus_v_t_ij_1Nminus1vector+d*sigma2_tp[i])-
                             exp(eta_t+rho2_t/2-DistSquare_u_tp_i_minus_v_t_ij_1Nminus1vector/The_1_plus_2sigma2_tp_i_plus_2varphi2_t_ij_1Nminus1vector)/(The_1_plus_2sigma2_tp_i_plus_2varphi2_t_ij_1Nminus1vector^(0.5*d))-
                             0.5*a_t[i]/b_t[i]*colSums((V_t[,-i,i]-u_tp_i_dNminus1matrix)^2))+
            0.5*d*(log(sigma2_tp[i])-(N-1)*a_t[i]/b_t[i]*sigma2_tp[i]) - 0.5*sum(pi_t_ik_times_xi_t_k_over_psi_t_k_1Kvector*(colSums((matrix(U_tp[i,],d,K)-t_mu_t)^2)+d*sigma2_tp[i]))
          # Determine increase or not
          if (ELBO_tp_U>=ELBO_t_U & !is.nan(ELBO_tp_U)){Increase_U<-TRUE}else{epsilon_U[i] <- epsilon_U[i]/2}
        }
      }
      
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      # Step 2: update V_t and varphi2_t by NGA
      for (i in 1:N){
        for (j in 1:N){
          if (i != j){
            Increase_V <- FALSE
            u_tp_i_minus_v_t_ij <- U_tp[i,]-V_t[,i,j]
            DistSquare_u_tp_i_minus_v_t_ij <- sum(u_tp_i_minus_v_t_ij^2)
            The_1_plus_2sigma2_tp_i_plus_2varphi2_t_ij <- 1+2*(sigma2_tp[i]+varphi2_t[i,j])
            DistSqaure_over_1_plus_2var_t <- DistSquare_u_tp_i_minus_v_t_ij/The_1_plus_2sigma2_tp_i_plus_2varphi2_t_ij
            The_exp_term_t <- exp(eta_t+rho2_t/2-DistSqaure_over_1_plus_2var_t)
            The_big_shared_term_t_V <- The_exp_term_t/(The_1_plus_2sigma2_tp_i_plus_2varphi2_t_ij^(0.5*d+1))
            v_t_ij_minus_u_tp_j <- V_t[,i,j]-U_tp[j,]
            # Calculate the ELBO relevant to the current state V_t[,i,j] and varphi2_t[i,j]
            ELBO_t_V <- -Y[i,j]*(DistSquare_u_tp_i_minus_v_t_ij+d*varphi2_t[i,j])-
              The_exp_term_t/(The_1_plus_2sigma2_tp_i_plus_2varphi2_t_ij^(0.5*d))-
              0.5*a_t[j]/b_t[j]*(sum(v_t_ij_minus_u_tp_j^2)+d*varphi2_t[i,j])+0.5*d*log(varphi2_t[i,j])
            while (Increase_V == FALSE){
              # Update V_tp
              V_tp[,i,j] <- V_t[,i,j]+epsilon_V[i,j]*varphi2_t[i,j]*
                (2*(Y[i,j]-The_big_shared_term_t_V)*u_tp_i_minus_v_t_ij-a_t[j]/b_t[j]*v_t_ij_minus_u_tp_j)
              # Update varphi2_tp
              varphi2_tp[i,j] <- varphi2_t[i,j]*exp(epsilon_V[i,j]*2/d*varphi2_t[i,j]*(
                0.5*d*(-2*Y[i,j]-a_t[j]/b_t[j]+1/varphi2_t[i,j])-The_big_shared_term_t_V*(2*DistSqaure_over_1_plus_2var_t-d)))
              # Determine whether the ELBO increases or not; if no increase, half the epsilon and repeat again
              # Calculate the ELBO relevant to the update state V_tp[,i,j] and varphi2_tp[i,j]
              DistSquare_u_tp_i_minus_v_tp_ij <- sum((U_tp[i,]-V_tp[,i,j])^2)
              The_1_plus_2sigma2_tp_i_plus_2varphi2_tp_ij <- 1+2*(sigma2_tp[i]+varphi2_tp[i,j])
              ELBO_tp_V <- -Y[i,j]*(DistSquare_u_tp_i_minus_v_tp_ij+d*varphi2_tp[i,j])-
                exp(eta_t+rho2_t/2-DistSquare_u_tp_i_minus_v_tp_ij/The_1_plus_2sigma2_tp_i_plus_2varphi2_tp_ij)/(The_1_plus_2sigma2_tp_i_plus_2varphi2_tp_ij^(0.5*d))-
                0.5*a_t[j]/b_t[j]*(sum((V_tp[,i,j]-U_tp[j,])^2)+d*varphi2_tp[i,j])+0.5*d*log(varphi2_tp[i,j])
              # Determine increase or not
              if (ELBO_tp_V>=ELBO_t_V){Increase_V<-TRUE}else{epsilon_V[i,j] <- epsilon_V[i,j]/2}
            }
          }
        }
      }
      
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      # Step 3: update Pi_t by AS
      Ajk <- t(matrix(0.5*d*(digamma(xi_t)-log(psi_t))+digamma(delta_t)-digamma(delta_t0),K,N))-
        t(matrix(0.5*xi_t/psi_t,K,N))*(proxy::dist(U_tp, mu_t, method="Euclidean")^2+(matrix(sigma2_tp,N,K)+t(matrix(omega2_t,K,N)))*d)
      exp_Ajk <- exp(Ajk)
      Pi_tp <- exp_Ajk/matrix(rowSums(exp_Ajk),N,K)
      Pi_tp <- apply(Pi_tp,1:2,Check.Replace.Smallest.Number)
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      # Step 4: update a_t and b_t by AS
      a_tp <- rep(a+0.5*d*(N-1),N)
      DistSquare_v_tp_ij_minus_u_tp_j_NNmatrix <- apply((V_tp-array(do.call(rbind, replicate(N, t(U_tp), simplify=FALSE)),dim=c(d,N,N)))^2,2:3,sum)
      diag(DistSquare_v_tp_ij_minus_u_tp_j_NNmatrix) <- 0
      varphi2_tp_ij_plus_sigma2_tp_i_NNmatrix <- varphi2_tp+t(matrix(sigma2_tp,N,N))
      diag(varphi2_tp_ij_plus_sigma2_tp_i_NNmatrix) <- 0
      b_tp <- b+0.5*(colSums(DistSquare_v_tp_ij_minus_u_tp_j_NNmatrix)+d*colSums(varphi2_tp_ij_plus_sigma2_tp_i_NNmatrix))
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      # Step 5: update eta_t and rho2_t by NGA
      Increase_beta <- FALSE
      Sum_Yij <- Y; diag(Sum_Yij) <- 0; Sum_Yij <- sum(Sum_Yij)
      The_1_plus_2sigma2_tp_i_plus_2varphi2_tp_ij_NNmatrix <- 1+2*(matrix(sigma2_tp,N,N)+varphi2_tp)
      The_1_plus_2var_tp_power0.5d_NNmatrix <- The_1_plus_2sigma2_tp_i_plus_2varphi2_tp_ij_NNmatrix^(0.5*d)
      DistSquare_over_1_plus_2var_tp_NNmatrix <- apply((array(t(U_tp),dim=c(d,N,N))-V_tp)^2,2:3,sum)/The_1_plus_2sigma2_tp_i_plus_2varphi2_tp_ij_NNmatrix
      The_exp_term_t_NNmatrix <- exp(eta_t+rho2_t/2-DistSquare_over_1_plus_2var_tp_NNmatrix)
      The_big_shared_term_t_beta <- The_exp_term_t_NNmatrix/The_1_plus_2var_tp_power0.5d_NNmatrix
      diag(The_big_shared_term_t_beta) <- 0 # since diagonal entries are NAs; recall: diagonal entries of V and varphi are NAs
      Sum_big_shared_term_t_beta <- sum(The_big_shared_term_t_beta)
      # Calculate the ELBO relevant to the current state eta_t and rho2_t
      ELBO_t_beta <- 0.5*(log(rho2_t)-1/rho2*((eta_t-eta)^2+rho2_t))+eta_t*Sum_Yij-Sum_big_shared_term_t_beta
      while (Increase_beta == FALSE){
        # Update eta_tp
        eta_tp <- eta_t+epsilon_beta*rho2_t*(Sum_Yij-(eta_t-eta)/rho2-Sum_big_shared_term_t_beta)
        # Update rho2_tp
        rho2_tp_inside_exp <- epsilon_beta*2*rho2_t*0.5*(1/rho2_t-1/rho2-Sum_big_shared_term_t_beta)
        rho2_tp <- rho2_t*exp(rho2_tp_inside_exp) # so the log(rho2_tp)=log(rho2_t)+rho2_tp_inside_exp in the ELBO; this aims to avoid log(exp(-1000))=-Inf in R
        # Determine whether the ELBO increases or not; if no increase, half the epsilon and repeat again
        # Calculate the ELBO relevant to the update state eta_tp and rho2_tp
        The_exp_over_1_plus_2var_power0.5d_tp_NNmatrix <- exp(eta_tp+rho2_tp/2-DistSquare_over_1_plus_2var_tp_NNmatrix)/The_1_plus_2var_tp_power0.5d_NNmatrix
        diag(The_exp_over_1_plus_2var_power0.5d_tp_NNmatrix) <- 0
        ELBO_tp_beta <- 0.5*((rho2_tp_inside_exp+log(rho2_t))-1/rho2*((eta_tp-eta)^2+rho2_tp))+eta_tp*Sum_Yij-sum(The_exp_over_1_plus_2var_power0.5d_tp_NNmatrix)
        # Determine increase or not
        if (ELBO_tp_beta>=ELBO_t_beta){Increase_beta<-TRUE}else{epsilon_beta <- epsilon_beta/2}
      }
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      # Step 6: update mu_t by AS
      colSums_Pi_tp <- colSums(Pi_tp)
      mu_tp <- matrix(omega2*xi_t/(omega2*xi_t*colSums_Pi_tp+psi_t),K,d)*
        t(apply(array(t(U_tp),dim=c(d,N,K))*array(t(matrix(do.call(rbind, replicate(d, t(Pi_tp), simplify=FALSE)),K,N*d)),dim=c(d,N,K)),c(1,3),sum))
      # t(matrix(do.call(rbind, replicate(d, t(Pi_tp), simplify=FALSE)),K,N*d)) # this aims to construct a N*d X K matrix where, 
      # e.g. when d=2 we have that the 1,1th, 2,1th entries are the same and are Pi_tp[1,1]; similarly, the 3,1th, 4,1th entries are the same and are Pi_tp[2,1]; 
      # the 1,2th, 2,2th entries are the same and are Pi_tp[1,2]; the 3,2th, 4,2th entries are the same and are Pi_tp[2,2], and so on
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      # Step 7: update omega2_t by AS
      omega2_tp <- omega2*psi_t/(psi_t+omega2*xi_t*colSums_Pi_tp)
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      # Step 8: update xi_t,psi_t by AS
      xi_tp <- xi +0.5*d*colSums_Pi_tp
      psi_tp <- psi + 0.5*colSums(Pi_tp*(proxy::dist(U_tp, mu_tp, method="Euclidean")^2+(matrix(sigma2_tp,N,K)+t(matrix(omega2_tp,K,N)))*d))
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      # Step 9: update delta_t by SGA
      Increase_Pi <- FALSE
      # Calculate the ELBO relevant to the current state delta_t
      ELBO_t_Pi <- -lgamma(delta_t0) + sum(Pi_tp*(t(matrix(digamma(delta_t),K,N)))) - digamma(delta_t0)*N + sum((delta-delta_t)*(digamma(delta_t)-digamma(delta_t0))+lgamma(delta_t))
      while (Increase_Pi == FALSE){
        # Update all delta_t simultaneously
        delta_tp <- delta_t*exp(epsilon_Pi*delta_t*(trigamma(delta_t)*(delta-delta_t+colSums_Pi_tp)-trigamma(delta_t0)*(sum(delta-delta_t)+N)))
        delta_tp0 <- sum(delta_tp)
        # Calculate the ELBO relevant to the update state delta_tp
        ELBO_tp_Pi <- -lgamma(delta_tp0) + sum(Pi_tp*(t(matrix(digamma(delta_tp),K,N)))) - digamma(delta_tp0)*N + sum((delta-delta_tp)*(digamma(delta_tp)-digamma(delta_tp0))+lgamma(delta_tp))
        # Determine increase or not
        if (ELBO_tp_Pi>=ELBO_t_Pi & (lgamma(delta_tp0)+1)!=lgamma(delta_tp0) #& ELBO_tp_Pi!=0
            ){Increase_Pi<-TRUE}else{epsilon_Pi <- epsilon_Pi/2}
        # (lgamma(delta_tp0)+1)!=lgamma(delta_tp0) and ELBO_tp_Pi!=0 case here corresponds to the problem where the the calculations exceeds the capability of R programming, e.g., (lgamma(1.428950e+21)+1)==lgamma(1.428950e+21)
      }
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      #--------------------------------------------------------------------------------------------------------------------------------------------
      # Evaluate the ELBO at the update state
      # Calculate the 1st sum over i,j = 1,2,...,N and i != j
      DistSquare_u_tp_i_minus_v_tp_ij_NNmatrix <- apply((array(t(U_tp),dim=c(d,N,N))-V_tp)^2,2:3,sum) # NNmatrix means N X N matrix
      sigma2_tp_i_plus_varphi2_tp_ij_NNmatrix <- matrix(sigma2_tp,N,N)+varphi2_tp
      The_1_plus_2var_tp_NNmatrix <- 1+2*sigma2_tp_i_plus_varphi2_tp_ij_NNmatrix
      a_tp_j_NNmatrix <- t(matrix(a_tp,N,N))
      b_tp_j_NNmatrix <- t(matrix(b_tp,N,N))
      ELBO_tp_Sum_ij <- Y*(eta_tp-DistSquare_u_tp_i_minus_v_tp_ij_NNmatrix-d*sigma2_tp_i_plus_varphi2_tp_ij_NNmatrix)-
        exp(eta_tp+rho2_tp/2-DistSquare_u_tp_i_minus_v_tp_ij_NNmatrix/The_1_plus_2var_tp_NNmatrix)/(The_1_plus_2var_tp_NNmatrix^(0.5*d))+
        0.5*d*(digamma(a_tp_j_NNmatrix)-log(b_tp_j_NNmatrix)+log(varphi2_tp))-
        0.5*a_tp_j_NNmatrix/b_tp_j_NNmatrix*(apply((V_tp-array(do.call(rbind, replicate(N, t(U_tp), simplify=FALSE)),dim=c(d,N,N)))^2,2:3,sum)+d*(t(matrix(sigma2_tp,N,N))+varphi2_tp))
      diag(ELBO_tp_Sum_ij) <- 0
      # Calculate the 2nd sum over i = 1,2,...,N and k = 1,2,...,K
      xi_tp_k_NKmatrix <- t(matrix(xi_tp,K,N))
      psi_tp_k_NKmatrix <- t(matrix(psi_tp,K,N))
      ELBO_tp_Sum_ik <-
        Pi_tp*(0.5*d*(digamma(xi_tp_k_NKmatrix)-log(psi_tp_k_NKmatrix))-0.5*xi_tp_k_NKmatrix/psi_tp_k_NKmatrix* # (apply((array(t(U_tp),dim=c(d,N,K))-array(do.call(rbind, replicate(N, t(mu_tp), simplify=FALSE)),dim=c(d,N,K)))^2,2:3,sum) # array() way to calculate ||u_i-mu_k||^2
                 (proxy::dist(U_tp, mu_tp, method="Euclidean")^2+d*(matrix(sigma2_tp,N,K)+t(matrix(omega2_tp,K,N))))+
                 digamma(t(matrix(delta_tp,K,N)))-digamma(delta_tp0)-log(Pi_tp))
      # Calculate the 3rd sum over k = 1,2,...,K
      ELBO_tp_Sum_k <- 
        (delta-delta_tp)*(digamma(delta_tp)-digamma(delta_tp0))+lgamma(delta_tp)-0.5/omega2*(rowSums(mu_tp^2)+d*omega2_tp)+0.5*d*log(omega2_tp)+
        (xi-xi_tp)*digamma(xi_tp)-xi*log(psi_tp)-psi*xi_tp/psi_tp+lgamma(xi_tp)+xi_tp
      # Calculate the 4th sum over i = 1,2,...,N
      ELBO_tp_Sum_i <- (a-a_tp)*digamma(a_tp)-a*log(b_tp)-b*a_tp/b_tp+lgamma(a_tp)+a_tp+0.5*d*log(sigma2_tp)
      # Calculate the ELBO with the rest single terms
      ELBO_tp <- sum(ELBO_tp_Sum_ij)+sum(ELBO_tp_Sum_ik)+sum(ELBO_tp_Sum_i)+sum(ELBO_tp_Sum_k)-
        0.5/rho2*((eta_tp-eta)^2+rho2_tp)+0.5*log(rho2_tp)-lgamma(delta_tp0)
      #--------------------------------------------------------------------------------------------------------------------------------------------
      # If ELBO(theta_tp)-ELBO(theta_t)<tol, then stop the algorithm
      if (ELBO_tp-ELBO_t<=tol){STOP <- TRUE}
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      # Update the current state
      U_t<-U_tp; sigma2_t<-sigma2_tp; V_t<-V_tp; varphi2_t<-varphi2_tp; Pi_t<-Pi_tp; a_t<-a_tp; b_t<-b_tp; eta_t<-eta_tp; rho2_t<-rho2_tp
      mu_t<-mu_tp; omega2_t<-omega2_tp; xi_t<-xi_tp; psi_t<-psi_tp; delta_t<-delta_tp; delta_t0<-delta_tp0
      ELBO_list <- c(ELBO_list,ELBO_tp); ELBO_t <- ELBO_tp
    }
    return(list(NumIt=t,ELBO_list=ELBO_list, U_t=U_t,sigma2_t=sigma2_t, V_t=V_t,varphi2_t=varphi2_t, Pi_t=Pi_t, a_t=a_t,b_t=b_t, eta_t=eta_t,rho2_t=rho2_t,
                mu_t=mu_t,omega2_t=omega2_t, xi_t=xi_t,psi_t=psi_t, delta_t=delta_t,
                epsilon_U=epsilon_U,epsilon_V=epsilon_V,epsilon_beta=epsilon_beta,epsilon_Pi=epsilon_Pi))
  }

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# VB inference from the new MixedLPCM
# But with specified prior parameters a

VB_Directed_MixedLPCM_VecPrior_a <- 
  function(Y,d=2, a=rep(1,nrow(Y)),b=1,eta=1,rho2=1,omega2=1,xi=1,psi=1,K=10,delta=1, # delta=rep(1e-3,K), # data, dimension and prior parameters
           U_t=NA,sigma2_t=NA,V_t=NA,varphi2_t=NA,Pi_t=NA,eta_t=NA,rho2_t=NA,mu_t=NA,omega2_t=NA,xi_t=NA,psi_t=NA,a_t=NA,b_t=NA,delta_t=NA,  # initialize VB parameters
           tol=1e-2,epsilon_init=1,seed=NULL){ # VB algorithm settings
    
    # Y: the N X N observed directed binary network adjacency matrix
    # d: the dimension assumed for latent position space
    # a,b: prior parameters for each \gamma_i; default setting Ga(1,1)
    # eta,rho2: prior parameters for \beta; default setting N(0,1)
    # omega2: prior parameter for each \bm{\mu}_k; default setting MVN_d(0,1)
    # xi,psi: prior parameters for each \tau_k; default setting Ga(1,1)
    # K: an upper bound of the number of non-empty groups; default setting K=20; 
    # extra groups will be made redundant automatically by very small \tilde{\pi_ik}
    # delta: a 1 X K vector or a single value if a symmetric Dirichlet is considered; prior parameters for \bm{\Pi}; 
    # default setting Dirichlet(rep(1e-3,K)); Shrinkage prior required that each entry of delta should be < 1/K (Rousseau and Mengersen, 2011)
    
    # tol: tolerance level; if the ELBO^{(t+1)}-ELBO^{(t)} < tol, then terminate the implementation; default setting tol=1e-5
    # epsilon_init: initial epsilon set for each the gradient ascent step; default setting epsilon=1
    # seed: the seed for random number generator
    
    set.seed(seed)
    N <- nrow(Y)
    Check.Replace.Smallest.Number <- function(x) max(x,.Machine$double.xmin)
    Check.Replace.Smallest.Number <- Vectorize(Check.Replace.Smallest.Number)
    # It's possible to have very small numbers during the inference and continuing the process would make the corresponding model parameters become too small so that the machine cannot represent and a zero instead appears
    # This will bring the logarithm of such a parameter being -Inf. So we limit the smallest number to .Machine$double.xmin
    
    # U_t: \tilde{U}, a N X d matrix
    # sigma2_t: \tilde{sigma2}, a 1 X N vector
    # V_t: \tilde{V}, a d X N X N array
    # varphi2_t: \tilde{varphi}, a N X N matrix
    # Pi_t: \tilde{Pi}, a N X K matrix; each row is the 1 X K vector \tilde{\Pi}_i
    # eta_t and rho2_t: \tilde{eta} and \tilde{rho2}, two normal distribution parameters
    # mu_t: \tilde{mu}, a K X d matrix
    # omega2_t: \tilde{omega2}, a 1 X K vector
    # xi_t: \tilde{xi}, a 1 X K vector
    # psi_t: \tilde{psi}, a 1 X K vector
    # a_t: \tilde{a}, a 1 X N vector
    # b_t: \tilde{b}, a 1 X N vector
    # delta_t: \tilde{delta}, a 1 X K vector
    
    # We instead denote ( U_tp,sigma2_tp,V_tp,varphi2_tp,Pi_tp,eta_tp,rho2_tp,mu_tp,omega2_tp,xi_tp,psi_tp,a_tp,b_tp,delta_tp ) as the new update
    
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    # Initialization
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize covert latent positions U_t by Classical (Metric) Multidimensional Scaling (MDS) on asymmetric distance matrix
    library(ergmito) # for geodesic
    Y_GM <- geodesic(Y,force=TRUE)[[1]] # geodesic matrix (GM) of Y, the length/the number of edges of the shortest paths between each pair of nodes
    Y_GM[is.na(Y_GM)] <- 0 # remove NA as 0
    if (is.na(U_t[1])){U_t <- cmdscale(Y_GM,d)} # obtain the initial latent positions by MDS; cmdscale(,d) returns a N X d matrix
    # Note here that MDS may bring the same position for multiple nodes
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize overt latent positions V_t by setting all { v_{j<--i} }^N_{i=1,i!=j} be the same as the u_j for j = 1,2,...,N
    if (is.na(V_t[1])){
      V_t <- array(0,dim=c(d,N,N))
      for (j in 1:N){
        V_t[,,j] <- matrix(U_t[j,],d,N)
        V_t[,j,j] <- rep(NA,d)
      }
    }
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize the group MVN distribution center mu_t by first applying kmeans() on the covert latent positions to obtain the initial grouping,
    # and then let each \tilde{mu}_k be the mean of the latent positions contained in group k
    z_init <- kmeans(U_t,centers=K,nstart = 1000)$cluster # check 1000 different start and pick the best one
    if (is.na(mu_t[1])){
      mu_t <- matrix(0,K,d)
      for (k in 1:K){mu_t[k,] <- apply(matrix(U_t[z_init==k,],ncol=d),2,mean)}
    }
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize sigma2_t, a_t and b_t
    if (is.na(sigma2_t[1])){sigma2_t <- rep(1,N)}
    if (is.na(a_t[1])){a_t <- a}
    if (is.na(b_t[1])){b_t <- rep(b,N)}
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize varphi2_t
    if (is.na(varphi2_t[2])){ # Check whether varphi2_t is NA or a matrix with diag=NA (so the 2nd element is instead checked)
      varphi2_t <- matrix(1,N,N)
      diag(varphi2_t) <- NA
    }
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize omega2_t, xi_t, psi_t and delta_t
    if (is.na(omega2_t[1])){omega2_t <- rep(omega2,K)}
    if (is.na(xi_t[1])){xi_t <- rep(xi,K)}
    if (is.na(psi_t[1])){psi_t <- rep(psi,K)}
    if (is.na(delta_t[1])){delta_t <- rep(delta,K)}
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize eta_t and rho2_t
    if (is.na(eta_t[1])){eta_t <- eta}
    if (is.na(rho2_t[1])){rho2_t <- rho2}
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize Pi_t by evaluating the analytical solution based on all other initial settings
    library(proxy) # for proxy::dist; pairwise distance between two different set of vectors
    delta_t0 <- sum(delta_t)
    if (is.na(Pi_t)){
      # We first define the A_{jk} which is a N X K matrix
      Ajk <- t(matrix(0.5*d*(digamma(xi_t)-log(psi_t))+digamma(delta_t)-digamma(delta_t0),K,N))-
        t(matrix(0.5*xi_t/psi_t,K,N))*(proxy::dist(U_t, mu_t, method="Euclidean")^2+(matrix(sigma2_t,N,K)+t(matrix(omega2_t,K,N)))*d)
      exp_Ajk <- exp(Ajk)
      Pi_t <- exp_Ajk/matrix(rowSums(exp_Ajk),N,K) # a N X K matrix with each row being \Pi_i
      Pi_t <- apply(Pi_t,1:2,Check.Replace.Smallest.Number)
    }
    
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize the epsilon for gradient ascent steps
    epsilon_U <- rep(epsilon_init,N)
    epsilon_V <- matrix(epsilon_init,N,N); diag(epsilon_V) <- NA
    epsilon_beta <- epsilon_Pi <- epsilon_init
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Evaluate the ELBO at the initial state
    # Calculate the 1st sum over i,j = 1,2,...,N and i != j
    DistSquare_u_t_i_minus_v_t_ij_NNmatrix <- apply((array(t(U_t),dim=c(d,N,N))-V_t)^2,2:3,sum) # NNmatrix means N X N matrix
    sigma2_t_i_plus_varphi2_t_ij_NNmatrix <- matrix(sigma2_t,N,N)+varphi2_t
    The_1_plus_2var_t_NNmatrix <- 1+2*sigma2_t_i_plus_varphi2_t_ij_NNmatrix
    a_t_j_NNmatrix <- t(matrix(a_t,N,N))
    b_t_j_NNmatrix <- t(matrix(b_t,N,N))
    ELBO_t_Sum_ij <- Y*(eta_t-DistSquare_u_t_i_minus_v_t_ij_NNmatrix-d*sigma2_t_i_plus_varphi2_t_ij_NNmatrix)-
      exp(eta_t+rho2_t/2-DistSquare_u_t_i_minus_v_t_ij_NNmatrix/The_1_plus_2var_t_NNmatrix)/(The_1_plus_2var_t_NNmatrix^(0.5*d))+
      0.5*d*(digamma(a_t_j_NNmatrix)-log(b_t_j_NNmatrix)+log(varphi2_t))-
      0.5*a_t_j_NNmatrix/b_t_j_NNmatrix*(apply((V_t-array(do.call(rbind, replicate(N, t(U_t), simplify=FALSE)),dim=c(d,N,N)))^2,2:3,sum)+d*(t(matrix(sigma2_t,N,N))+varphi2_t))
    diag(ELBO_t_Sum_ij) <- 0
    # Calculate the 2nd sum over i = 1,2,...,N and k = 1,2,...,K
    xi_t_k_NKmatrix <- t(matrix(xi_t,K,N))
    psi_t_k_NKmatrix <- t(matrix(psi_t,K,N))
    ELBO_t_Sum_ik <-
      Pi_t*(0.5*d*(digamma(xi_t_k_NKmatrix)-log(psi_t_k_NKmatrix))-0.5*xi_t_k_NKmatrix/psi_t_k_NKmatrix* # (apply((array(t(U_t),dim=c(d,N,K))-array(do.call(rbind, replicate(N, t(mu_t), simplify=FALSE)),dim=c(d,N,K)))^2,2:3,sum) # array() way to calculate ||u_i-mu_k||^2
              (proxy::dist(U_t, mu_t, method="Euclidean")^2+d*(matrix(sigma2_t,N,K)+t(matrix(omega2_t,K,N))))+
              digamma(t(matrix(delta_t,K,N)))-digamma(delta_t0)-log(Pi_t))
    # Calculate the 3rd sum over k = 1,2,...,K
    ELBO_t_Sum_k <- 
      (delta-delta_t)*(digamma(delta_t)-digamma(delta_t0))+lgamma(delta_t)-0.5/omega2*(rowSums(mu_t^2)+d*omega2_t)+0.5*d*log(omega2_t)+
      (xi-xi_t)*digamma(xi_t)-xi*log(psi_t)-psi*xi_t/psi_t+lgamma(xi_t)+xi_t
    # Calculate the 4th sum over i = 1,2,...,N
    ELBO_t_Sum_i <- (a-a_t)*digamma(a_t)-a*log(b_t)-b*a_t/b_t+lgamma(a_t)+a_t+0.5*d*log(sigma2_t)
    # Calculate the ELBO with the rest single terms
    ELBO_t <- sum(ELBO_t_Sum_ij)+sum(ELBO_t_Sum_ik)+sum(ELBO_t_Sum_k)+sum(ELBO_t_Sum_i)-
      0.5/rho2*((eta_t-eta)^2+rho2_t)+0.5*log(rho2_t)-lgamma(delta_t0)
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize the ELBO list
    ELBO_list <- c(ELBO_t)
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize the update state
    U_tp<-U_t; sigma2_tp<-sigma2_t; V_tp<-V_t; varphi2_tp<-varphi2_t; Pi_tp<-Pi_t; a_tp<-a_t; b_tp<-b_t; eta_tp<-eta_t; rho2_tp<-rho2_t
    mu_tp<-mu_t; omega2_tp<-omega2_t; xi_tp<-xi_t; psi_tp<-psi_t; delta_tp<-delta_t
    
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    # Inference steps
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    t <- 0; STOP <- FALSE
    while (STOP == FALSE){
      t <- t + 1
      if ((t%%100) == 0){cat("Iteration:", t,"\n")}
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      # Step 1: update U_t and sigma2_t by NGA
      for (i in 1:N){
        Increase_U <- FALSE
        # Define the terms which are calculated for more than twice in order to reduce the computation burden
        u_t_i_dNminus1matrix <- matrix(U_t[i,],d,N-1) # a d X N-1 matrix with each column being the \tilde{u}_i
        u_t_i_minus_v_t_ij_dNminus1matrix <- u_t_i_dNminus1matrix-V_t[,i,-i] # matrix(V_t[,i,-i],nrow=d)
        DistSquare_u_t_i_minus_v_t_ij_1Nminus1vector <- colSums(u_t_i_minus_v_t_ij_dNminus1matrix^2)
        The_1_plus_2sigma2_t_i_plus_2varphi2_t_ij_1Nminus1vector <- 1+2*(sigma2_t[i]+varphi2_t[i,-i])
        DistSqaure_over_1_plus_2var_t_1Nminus1vector <- DistSquare_u_t_i_minus_v_t_ij_1Nminus1vector/The_1_plus_2sigma2_t_i_plus_2varphi2_t_ij_1Nminus1vector
        The_exp_term_t_1Nminus1vector <- exp(eta_t+rho2_t/2-DistSqaure_over_1_plus_2var_t_1Nminus1vector)
        The_big_shared_term_t_U_1Nminus1vector <- The_exp_term_t_1Nminus1vector/(The_1_plus_2sigma2_t_i_plus_2varphi2_t_ij_1Nminus1vector^(0.5*d+1))
        v_t_ji_minus_u_t_i_dNminus1matrix <- V_t[,-i,i]-u_t_i_dNminus1matrix
        t_mu_t <- t(mu_t)
        u_t_i_minus_mu_t_k_dKmatrix <- matrix(U_t[i,],d,K)-t_mu_t
        pi_t_ik_times_xi_t_k_over_psi_t_k_1Kvector <- Pi_t[i,]*xi_t/psi_t
        # Calculate the ELBO relevant to the current state U_t[i,] and sigma2_t[i]
        ELBO_t_U <- sum(-Y[i,-i]*(DistSquare_u_t_i_minus_v_t_ij_1Nminus1vector+d*sigma2_t[i])-The_exp_term_t_1Nminus1vector/(The_1_plus_2sigma2_t_i_plus_2varphi2_t_ij_1Nminus1vector^(0.5*d))-
                          0.5*a_t[i]/b_t[i]*colSums(v_t_ji_minus_u_t_i_dNminus1matrix^2))+
          0.5*d*(log(sigma2_t[i])-(N-1)*a_t[i]/b_t[i]*sigma2_t[i]) - 0.5*sum(pi_t_ik_times_xi_t_k_over_psi_t_k_1Kvector*(colSums(u_t_i_minus_mu_t_k_dKmatrix^2)+d*sigma2_t[i]))
        # Set the while loop to check whether the ELBO increases or not
        while (Increase_U == FALSE){
          # Update U_tp
          U_tp[i,] <- U_t[i,]+epsilon_U[i]*sigma2_t[i]*
            (rowSums(t(replicate(d,2*(-Y[i,-i]+The_big_shared_term_t_U_1Nminus1vector)))*u_t_i_minus_v_t_ij_dNminus1matrix + a_t[i]/b_t[i]*v_t_ji_minus_u_t_i_dNminus1matrix)-
               rowSums(u_t_i_minus_mu_t_k_dKmatrix*t(matrix(replicate(d,pi_t_ik_times_xi_t_k_over_psi_t_k_1Kvector),K,d)))) # matrix(,K,d) here aims to deal with the K=1 case
          # Update sigma2_tp
          sigma2_tp[i] <- sigma2_t[i]*exp(epsilon_U[i]*2/d*sigma2_t[i]*
                                            (0.5*d*(1/sigma2_t[i]-(N-1)*a_t[i]/b_t[i]-sum(pi_t_ik_times_xi_t_k_over_psi_t_k_1Kvector))-
                                               sum(d*Y[i,-i]+The_big_shared_term_t_U_1Nminus1vector*(2*DistSqaure_over_1_plus_2var_t_1Nminus1vector-d))))
          if (is.nan(sigma2_tp[i])){sigma2_tp[i] <- sigma2_t[i]} # it's possible that some initial state would first bring very very very small sigma2_t leading to NaN; so in this case we lower-bound the value
          # Determine whether the ELBO increases or not; if no increase, half the epsilon and repeat again
          # Calculate the ELBO relevant to the update state U_tp[i,] and sigma2_tp[i]
          u_tp_i_dNminus1matrix <- matrix(U_tp[i,],d,N-1)
          DistSquare_u_tp_i_minus_v_t_ij_1Nminus1vector <- colSums((u_tp_i_dNminus1matrix-V_t[,i,-i])^2)
          The_1_plus_2sigma2_tp_i_plus_2varphi2_t_ij_1Nminus1vector <- 1+2*(sigma2_tp[i]+varphi2_t[i,-i])
          ELBO_tp_U <- sum(-Y[i,-i]*(DistSquare_u_tp_i_minus_v_t_ij_1Nminus1vector+d*sigma2_tp[i])-
                             exp(eta_t+rho2_t/2-DistSquare_u_tp_i_minus_v_t_ij_1Nminus1vector/The_1_plus_2sigma2_tp_i_plus_2varphi2_t_ij_1Nminus1vector)/(The_1_plus_2sigma2_tp_i_plus_2varphi2_t_ij_1Nminus1vector^(0.5*d))-
                             0.5*a_t[i]/b_t[i]*colSums((V_t[,-i,i]-u_tp_i_dNminus1matrix)^2))+
            0.5*d*(log(sigma2_tp[i])-(N-1)*a_t[i]/b_t[i]*sigma2_tp[i]) - 0.5*sum(pi_t_ik_times_xi_t_k_over_psi_t_k_1Kvector*(colSums((matrix(U_tp[i,],d,K)-t_mu_t)^2)+d*sigma2_tp[i]))
          # Determine increase or not
          if (ELBO_tp_U>=ELBO_t_U & !is.nan(ELBO_tp_U)){Increase_U<-TRUE}else{epsilon_U[i] <- epsilon_U[i]/2}
        }
      }
      
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      # Step 2: update V_t and varphi2_t by NGA
      for (i in 1:N){
        for (j in 1:N){
          if (i != j){
            Increase_V <- FALSE
            u_tp_i_minus_v_t_ij <- U_tp[i,]-V_t[,i,j]
            DistSquare_u_tp_i_minus_v_t_ij <- sum(u_tp_i_minus_v_t_ij^2)
            The_1_plus_2sigma2_tp_i_plus_2varphi2_t_ij <- 1+2*(sigma2_tp[i]+varphi2_t[i,j])
            DistSqaure_over_1_plus_2var_t <- DistSquare_u_tp_i_minus_v_t_ij/The_1_plus_2sigma2_tp_i_plus_2varphi2_t_ij
            The_exp_term_t <- exp(eta_t+rho2_t/2-DistSqaure_over_1_plus_2var_t)
            The_big_shared_term_t_V <- The_exp_term_t/(The_1_plus_2sigma2_tp_i_plus_2varphi2_t_ij^(0.5*d+1))
            v_t_ij_minus_u_tp_j <- V_t[,i,j]-U_tp[j,]
            # Calculate the ELBO relevant to the current state V_t[,i,j] and varphi2_t[i,j]
            ELBO_t_V <- -Y[i,j]*(DistSquare_u_tp_i_minus_v_t_ij+d*varphi2_t[i,j])-
              The_exp_term_t/(The_1_plus_2sigma2_tp_i_plus_2varphi2_t_ij^(0.5*d))-
              0.5*a_t[j]/b_t[j]*(sum(v_t_ij_minus_u_tp_j^2)+d*varphi2_t[i,j])+0.5*d*log(varphi2_t[i,j])
            while (Increase_V == FALSE){
              # Update V_tp
              V_tp[,i,j] <- V_t[,i,j]+epsilon_V[i,j]*varphi2_t[i,j]*
                (2*(Y[i,j]-The_big_shared_term_t_V)*u_tp_i_minus_v_t_ij-a_t[j]/b_t[j]*v_t_ij_minus_u_tp_j)
              # Update varphi2_tp
              varphi2_tp[i,j] <- varphi2_t[i,j]*exp(epsilon_V[i,j]*2/d*varphi2_t[i,j]*(
                0.5*d*(-2*Y[i,j]-a_t[j]/b_t[j]+1/varphi2_t[i,j])-The_big_shared_term_t_V*(2*DistSqaure_over_1_plus_2var_t-d)))
              # Determine whether the ELBO increases or not; if no increase, half the epsilon and repeat again
              # Calculate the ELBO relevant to the update state V_tp[,i,j] and varphi2_tp[i,j]
              DistSquare_u_tp_i_minus_v_tp_ij <- sum((U_tp[i,]-V_tp[,i,j])^2)
              The_1_plus_2sigma2_tp_i_plus_2varphi2_tp_ij <- 1+2*(sigma2_tp[i]+varphi2_tp[i,j])
              ELBO_tp_V <- -Y[i,j]*(DistSquare_u_tp_i_minus_v_tp_ij+d*varphi2_tp[i,j])-
                exp(eta_t+rho2_t/2-DistSquare_u_tp_i_minus_v_tp_ij/The_1_plus_2sigma2_tp_i_plus_2varphi2_tp_ij)/(The_1_plus_2sigma2_tp_i_plus_2varphi2_tp_ij^(0.5*d))-
                0.5*a_t[j]/b_t[j]*(sum((V_tp[,i,j]-U_tp[j,])^2)+d*varphi2_tp[i,j])+0.5*d*log(varphi2_tp[i,j])
              # Determine increase or not
              if (ELBO_tp_V>=ELBO_t_V){Increase_V<-TRUE}else{epsilon_V[i,j] <- epsilon_V[i,j]/2}
            }
          }
        }
      }
      
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      # Step 3: update Pi_t by AS
      Ajk <- t(matrix(0.5*d*(digamma(xi_t)-log(psi_t))+digamma(delta_t)-digamma(delta_t0),K,N))-
        t(matrix(0.5*xi_t/psi_t,K,N))*(proxy::dist(U_tp, mu_t, method="Euclidean")^2+(matrix(sigma2_tp,N,K)+t(matrix(omega2_t,K,N)))*d)
      exp_Ajk <- exp(Ajk)
      Pi_tp <- exp_Ajk/matrix(rowSums(exp_Ajk),N,K)
      Pi_tp <- apply(Pi_tp,1:2,Check.Replace.Smallest.Number)
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      # Step 4: update a_t and b_t by AS
      a_tp <- a+0.5*d*(N-1)
      DistSquare_v_tp_ij_minus_u_tp_j_NNmatrix <- apply((V_tp-array(do.call(rbind, replicate(N, t(U_tp), simplify=FALSE)),dim=c(d,N,N)))^2,2:3,sum)
      diag(DistSquare_v_tp_ij_minus_u_tp_j_NNmatrix) <- 0
      varphi2_tp_ij_plus_sigma2_tp_i_NNmatrix <- varphi2_tp+t(matrix(sigma2_tp,N,N))
      diag(varphi2_tp_ij_plus_sigma2_tp_i_NNmatrix) <- 0
      b_tp <- b+0.5*(colSums(DistSquare_v_tp_ij_minus_u_tp_j_NNmatrix)+d*colSums(varphi2_tp_ij_plus_sigma2_tp_i_NNmatrix))
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      # Step 5: update eta_t and rho2_t by NGA
      Increase_beta <- FALSE
      Sum_Yij <- Y; diag(Sum_Yij) <- 0; Sum_Yij <- sum(Sum_Yij)
      The_1_plus_2sigma2_tp_i_plus_2varphi2_tp_ij_NNmatrix <- 1+2*(matrix(sigma2_tp,N,N)+varphi2_tp)
      The_1_plus_2var_tp_power0.5d_NNmatrix <- The_1_plus_2sigma2_tp_i_plus_2varphi2_tp_ij_NNmatrix^(0.5*d)
      DistSquare_over_1_plus_2var_tp_NNmatrix <- apply((array(t(U_tp),dim=c(d,N,N))-V_tp)^2,2:3,sum)/The_1_plus_2sigma2_tp_i_plus_2varphi2_tp_ij_NNmatrix
      The_exp_term_t_NNmatrix <- exp(eta_t+rho2_t/2-DistSquare_over_1_plus_2var_tp_NNmatrix)
      The_big_shared_term_t_beta <- The_exp_term_t_NNmatrix/The_1_plus_2var_tp_power0.5d_NNmatrix
      diag(The_big_shared_term_t_beta) <- 0 # since diagonal entries are NAs; recall: diagonal entries of V and varphi are NAs
      Sum_big_shared_term_t_beta <- sum(The_big_shared_term_t_beta)
      # Calculate the ELBO relevant to the current state eta_t and rho2_t
      ELBO_t_beta <- 0.5*(log(rho2_t)-1/rho2*((eta_t-eta)^2+rho2_t))+eta_t*Sum_Yij-Sum_big_shared_term_t_beta
      while (Increase_beta == FALSE){
        # Update eta_tp
        eta_tp <- eta_t+epsilon_beta*rho2_t*(Sum_Yij-(eta_t-eta)/rho2-Sum_big_shared_term_t_beta)
        # Update rho2_tp
        rho2_tp_inside_exp <- epsilon_beta*2*rho2_t*0.5*(1/rho2_t-1/rho2-Sum_big_shared_term_t_beta)
        rho2_tp <- rho2_t*exp(rho2_tp_inside_exp) # so the log(rho2_tp)=log(rho2_t)+rho2_tp_inside_exp in the ELBO; this aims to avoid log(exp(-1000))=-Inf in R
        # Determine whether the ELBO increases or not; if no increase, half the epsilon and repeat again
        # Calculate the ELBO relevant to the update state eta_tp and rho2_tp
        The_exp_over_1_plus_2var_power0.5d_tp_NNmatrix <- exp(eta_tp+rho2_tp/2-DistSquare_over_1_plus_2var_tp_NNmatrix)/The_1_plus_2var_tp_power0.5d_NNmatrix
        diag(The_exp_over_1_plus_2var_power0.5d_tp_NNmatrix) <- 0
        ELBO_tp_beta <- 0.5*((rho2_tp_inside_exp+log(rho2_t))-1/rho2*((eta_tp-eta)^2+rho2_tp))+eta_tp*Sum_Yij-sum(The_exp_over_1_plus_2var_power0.5d_tp_NNmatrix)
        # print(c(epsilon_beta,ELBO_t_beta,ELBO_tp_beta))
        # Determine increase or not
        if (ELBO_tp_beta>=ELBO_t_beta){Increase_beta<-TRUE}else{epsilon_beta <- epsilon_beta/2}
      }
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      # Step 6: update mu_t by AS
      colSums_Pi_tp <- colSums(Pi_tp)
      mu_tp <- matrix(omega2*xi_t/(omega2*xi_t*colSums_Pi_tp+psi_t),K,d)*
        t(apply(array(t(U_tp),dim=c(d,N,K))*array(t(matrix(do.call(rbind, replicate(d, t(Pi_tp), simplify=FALSE)),K,N*d)),dim=c(d,N,K)),c(1,3),sum))
      # t(matrix(do.call(rbind, replicate(d, t(Pi_tp), simplify=FALSE)),K,N*d)) # this aims to construct a N*d X K matrix where, 
      # e.g. when d=2 we have that the 1,1th, 2,1th entries are the same and are Pi_tp[1,1]; similarly, the 3,1th, 4,1th entries are the same and are Pi_tp[2,1]; 
      # the 1,2th, 2,2th entries are the same and are Pi_tp[1,2]; the 3,2th, 4,2th entries are the same and are Pi_tp[2,2], and so on
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      # Step 7: update omega2_t by AS
      omega2_tp <- omega2*psi_t/(psi_t+omega2*xi_t*colSums_Pi_tp)
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      # Step 8: update xi_t,psi_t by AS
      xi_tp <- xi +0.5*d*colSums_Pi_tp
      psi_tp <- psi + 0.5*colSums(Pi_tp*(proxy::dist(U_tp, mu_tp, method="Euclidean")^2+(matrix(sigma2_tp,N,K)+t(matrix(omega2_tp,K,N)))*d))
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      # Step 9: update delta_t by SGA
      Increase_Pi <- FALSE
      # Calculate the ELBO relevant to the current state delta_t
      ELBO_t_Pi <- -lgamma(delta_t0) + sum(Pi_tp*(t(matrix(digamma(delta_t),K,N)))) - digamma(delta_t0)*N + sum((delta-delta_t)*(digamma(delta_t)-digamma(delta_t0))+lgamma(delta_t))
      while (Increase_Pi == FALSE){
        # Update all delta_t simultaneously
        delta_tp <- delta_t*exp(epsilon_Pi*delta_t*(trigamma(delta_t)*(delta-delta_t+colSums_Pi_tp)-trigamma(delta_t0)*(sum(delta-delta_t)+N)))
        delta_tp0 <- sum(delta_tp)
        # Calculate the ELBO relevant to the update state delta_tp
        ELBO_tp_Pi <- -lgamma(delta_tp0) + sum(Pi_tp*(t(matrix(digamma(delta_tp),K,N)))) - digamma(delta_tp0)*N + sum((delta-delta_tp)*(digamma(delta_tp)-digamma(delta_tp0))+lgamma(delta_tp))
        # Determine increase or not
        if (ELBO_tp_Pi>=ELBO_t_Pi & (lgamma(delta_tp0)+1)!=lgamma(delta_tp0) #& ELBO_tp_Pi!=0
            ){Increase_Pi<-TRUE}else{epsilon_Pi <- epsilon_Pi/2}
        # (lgamma(delta_tp0)+1)!=lgamma(delta_tp0) or ELBO_tp_Pi!=0 case here corresponds to the problem where the the calculations exceeds the capability of R programming, e.g., (lgamma(1.428950e+21)+1)==lgamma(1.428950e+21)
      }
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      #--------------------------------------------------------------------------------------------------------------------------------------------
      # Evaluate the ELBO at the update state
      # Calculate the 1st sum over i,j = 1,2,...,N and i != j
      DistSquare_u_tp_i_minus_v_tp_ij_NNmatrix <- apply((array(t(U_tp),dim=c(d,N,N))-V_tp)^2,2:3,sum) # NNmatrix means N X N matrix
      sigma2_tp_i_plus_varphi2_tp_ij_NNmatrix <- matrix(sigma2_tp,N,N)+varphi2_tp
      The_1_plus_2var_tp_NNmatrix <- 1+2*sigma2_tp_i_plus_varphi2_tp_ij_NNmatrix
      a_tp_j_NNmatrix <- t(matrix(a_tp,N,N))
      b_tp_j_NNmatrix <- t(matrix(b_tp,N,N))
      ELBO_tp_Sum_ij <- Y*(eta_tp-DistSquare_u_tp_i_minus_v_tp_ij_NNmatrix-d*sigma2_tp_i_plus_varphi2_tp_ij_NNmatrix)-
        exp(eta_tp+rho2_tp/2-DistSquare_u_tp_i_minus_v_tp_ij_NNmatrix/The_1_plus_2var_tp_NNmatrix)/(The_1_plus_2var_tp_NNmatrix^(0.5*d))+
        0.5*d*(digamma(a_tp_j_NNmatrix)-log(b_tp_j_NNmatrix)+log(varphi2_tp))-
        0.5*a_tp_j_NNmatrix/b_tp_j_NNmatrix*(apply((V_tp-array(do.call(rbind, replicate(N, t(U_tp), simplify=FALSE)),dim=c(d,N,N)))^2,2:3,sum)+d*(t(matrix(sigma2_tp,N,N))+varphi2_tp))
      diag(ELBO_tp_Sum_ij) <- 0
      # Calculate the 2nd sum over i = 1,2,...,N and k = 1,2,...,K
      xi_tp_k_NKmatrix <- t(matrix(xi_tp,K,N))
      psi_tp_k_NKmatrix <- t(matrix(psi_tp,K,N))
      ELBO_tp_Sum_ik <-
        Pi_tp*(0.5*d*(digamma(xi_tp_k_NKmatrix)-log(psi_tp_k_NKmatrix))-0.5*xi_tp_k_NKmatrix/psi_tp_k_NKmatrix* # (apply((array(t(U_tp),dim=c(d,N,K))-array(do.call(rbind, replicate(N, t(mu_tp), simplify=FALSE)),dim=c(d,N,K)))^2,2:3,sum) # array() way to calculate ||u_i-mu_k||^2
                 (proxy::dist(U_tp, mu_tp, method="Euclidean")^2+d*(matrix(sigma2_tp,N,K)+t(matrix(omega2_tp,K,N))))+
                 digamma(t(matrix(delta_tp,K,N)))-digamma(delta_tp0)-log(Pi_tp))
      # Calculate the 3rd sum over k = 1,2,...,K
      ELBO_tp_Sum_k <- 
        (delta-delta_tp)*(digamma(delta_tp)-digamma(delta_tp0))+lgamma(delta_tp)-0.5/omega2*(rowSums(mu_tp^2)+d*omega2_tp)+0.5*d*log(omega2_tp)+
        (xi-xi_tp)*digamma(xi_tp)-xi*log(psi_tp)-psi*xi_tp/psi_tp+lgamma(xi_tp)+xi_tp
      # Calculate the 4th sum over i = 1,2,...,N
      ELBO_tp_Sum_i <- (a-a_tp)*digamma(a_tp)-a*log(b_tp)-b*a_tp/b_tp+lgamma(a_tp)+a_tp+0.5*d*log(sigma2_tp)
      # Calculate the ELBO with the rest single terms
      ELBO_tp <- sum(ELBO_tp_Sum_ij)+sum(ELBO_tp_Sum_ik)+sum(ELBO_tp_Sum_i)+sum(ELBO_tp_Sum_k)-
        0.5/rho2*((eta_tp-eta)^2+rho2_tp)+0.5*log(rho2_tp)-lgamma(delta_tp0)
      #--------------------------------------------------------------------------------------------------------------------------------------------
      # If ELBO(theta_tp)-ELBO(theta_t)<tol, then stop the algorithm
      if (ELBO_tp-ELBO_t<=tol){STOP <- TRUE}
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      # Update the current state
      U_t<-U_tp; sigma2_t<-sigma2_tp; V_t<-V_tp; varphi2_t<-varphi2_tp; Pi_t<-Pi_tp; a_t<-a_tp; b_t<-b_tp; eta_t<-eta_tp; rho2_t<-rho2_tp
      mu_t<-mu_tp; omega2_t<-omega2_tp; xi_t<-xi_tp; psi_t<-psi_tp; delta_t<-delta_tp; delta_t0<-delta_tp0
      ELBO_list <- c(ELBO_list,ELBO_tp); ELBO_t <- ELBO_tp
    }
    return(list(NumIt=t,ELBO_list=ELBO_list, U_t=U_t,sigma2_t=sigma2_t, V_t=V_t,varphi2_t=varphi2_t, Pi_t=Pi_t, a_t=a_t,b_t=b_t, eta_t=eta_t,rho2_t=rho2_t,
                mu_t=mu_t,omega2_t=omega2_t, xi_t=xi_t,psi_t=psi_t, delta_t=delta_t,
                epsilon_U=epsilon_U,epsilon_V=epsilon_V,epsilon_beta=epsilon_beta,epsilon_Pi=epsilon_Pi))
  }

#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------

# SG2003 Label switching function
LabelSwitching_SG2003 <- 
  function(z=NULL, Z=NULL, vector_1byK=NULL,matrix_dbyK=NULL,matrix_KbyK=NULL,matrix2_KbyK=NULL){
    # z: the membership vector to be label-switched
    # Z: the membership matrix to be label-switched
    # vector_1byK: clustering dependent 1 X K vector to be label-switched
    # matrix_dbyK: clustering dependent d X K matrix to be label-switched
    # matrix_KbyK: clustering dependent K X K matrix to be label-switched
    # matrix2_KbyK: the 2nd clustering dependent K X K matrix to be label-switched
    
    if (!is.null(z)&is.null(Z)){
      N <- length(z); K <- max(z)
      Z <- t(t(matrix(z,N,K))==(1:K))*1 # matrix transformation
    }else if(is.null(z)&!is.null(Z)){
      N <- nrow(Z); K <- ncol(Z)
      z <- c(Z%*%1:K)
    }else if(!is.null(z)&!is.null(Z)){
      N <- length(z); K <- ncol(Z)
    }else{
      stop("Either Z or z should be provided!")
    }
    
    if (length(which(apply(Z,2,sum)==0))>0){ # if empty groups exists
      Z <- matrix(Z[,-which(apply(Z,2,sum)==0)],nrow=N) # remove the empty groups
      K <- ncol(Z)
      z <- c(Z%*%1:K)
    }
    
    # Construct the reassignment rule which is a K X 2 matrix; example: c(3,1) in a row means the initial group 3 is relabeled as group 1
    reassignment_rule <- cbind(order(apply(Z,2,which.max)),1:K)
    # First create a K X N matrix where each column is reassignment_rule[,1] and we apply t() to obtain the transport N X K matrix.
    # Then we compare each column of the transport with z and replace the output "true" by the corresponding element in the reassignment_rule[,2] N X K transport matrix
    # Finally, we apply row sums to obtain the relabeled clustering
    # z <- apply((t(matrix(reassignment_rule[,1],K,N))==z)*t(matrix(reassignment_rule[,2],K,N)),1,sum)
    # z <- c((t(matrix(reassignment_rule[,1],K,N))==z)%*%1:K)
    z <- c((t(matrix(order(apply(Z,2,which.max)),K,N))==z)%*%1:K)
    
    if (!is.null(vector_1byK)){
      vector_1byK <- vector_1byK[reassignment_rule[,1]]
    }
    if (!is.null(matrix_dbyK)){
      matrix_dbyK <- matrix_dbyK[,reassignment_rule[,1]]
    }
    if (!is.null(matrix_KbyK)){
      matrix_KbyK <- matrix_KbyK[reassignment_rule[,1],reassignment_rule[,1]]
    }
    if (!is.null(matrix2_KbyK)){
      matrix2_KbyK <- matrix2_KbyK[reassignment_rule[,1],reassignment_rule[,1]]
    }
    return(list(z=z, Z=Z, vector_1byK=vector_1byK,matrix_dbyK=matrix_dbyK,matrix_KbyK=matrix_KbyK,matrix2_KbyK=matrix2_KbyK))
  }

#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------

# The Partially Integrated Classification log-Likelihood (PICL) criteria for the MixedLPCM model selection

PICL_Directed_MixedLPCM <- function(Y,hat_U,hat_V,hat_z,K,
                                          a,b,omega2,delta,tol=.Machine$double.xmin){
  # Y: N X N adjacency matrix
  # hat_U: N X d estimated covert latent positions
  # hat_V: d X N X N estimated overt latent positions
  # hat_z: 1 X N estimated clustering
  # Note that hat_z should be label-switched
  N <- nrow(Y); d <- ncol(hat_U); diag(Y) <- 0
  hat_z <- LabelSwitching_SG2003(z=hat_z)$z
  table_hat_z <- table(factor(hat_z,levels=1:K))
  # a,b: prior parameters for \gamma_j ~ Ga(a,b)
  # omega2: prior parameters for \mu_k ~ MVN(0,omega2*I)
  # delta: prior parameters for \Pi ~ Dirichlet(delta,...,delta)
  # tol: the tolerance value to stop the gradient ascent step of max w.r.t. \tau_g
  
  # Calculate the 1st BIC approximation of the log(p(Y|hat_U,hat_V)) term
  DistSquare_hat_u_i_minus_hat_v_ij_NNmatrix <- apply((array(t(hat_U),dim=c(d,N,N))-hat_V)^2,2:3,sum)
  diag(DistSquare_hat_u_i_minus_hat_v_ij_NNmatrix) <- 0
  Exp_minus_DistSquare_hat_u_i_minus_hat_v_ij_NNmatrix <- exp(-DistSquare_hat_u_i_minus_hat_v_ij_NNmatrix)
  diag(Exp_minus_DistSquare_hat_u_i_minus_hat_v_ij_NNmatrix) <- 0 # exp(0)=1 in the diagonal
  Sum_Y <- sum(Y)
  
  PICL_Y_term <- -Sum_Y*log(sum(Exp_minus_DistSquare_hat_u_i_minus_hat_v_ij_NNmatrix))-sum(Y*DistSquare_hat_u_i_minus_hat_v_ij_NNmatrix)
  
  # Calculate the 2nd integrated log(p(hat_V|hat_U)) term
  DistSquare_hat_v_ij_minus_hat_u_j_NNmatrix <- apply((hat_V-array(do.call(rbind, replicate(N, t(hat_U), simplify=FALSE)),dim=c(d,N,N)))^2,2:3,sum)
  diag(DistSquare_hat_v_ij_minus_hat_u_j_NNmatrix) <- 0
  
  PICL_hat_V_term <- -sum((a+0.5*d*(N-1))*log(b+0.5*colSums(DistSquare_hat_v_ij_minus_hat_u_j_NNmatrix)))
  
  # Calculate the 3rd partial intergal + BIC approximation of the log(p(hat_U|hat_z)) term
  PICL_hat_U_term <- -0.5*K*log(N)
  max_w.r.t.tau_g_term <- function(tau_g,hat_n_g,d,omega2,DistSquare_Sum_hat_U_g,Sum_DistSquare_hat_U_g){ # hat_U_g := {hat_u_i: hat_z_i=g}
    return(0.5*d*(hat_n_g*log(tau_g)-log(tau_g*hat_n_g*omega2+1))+
      0.5*tau_g^2*omega2/(tau_g*hat_n_g*omega2+1)*DistSquare_Sum_hat_U_g-
      0.5*tau_g*Sum_DistSquare_hat_U_g)
  }
  tau_list <- c()
  for (g in 1:K){
    STOP <- FALSE; epsilon <- 1; tau_g_t <- 1; hat_n_g <- table_hat_z[[g]]
    hat_U_g <- matrix(hat_U[hat_z==g,],ncol=2); DistSquare_Sum_hat_U_g <- sum(colSums(hat_U_g)^2); Sum_DistSquare_hat_U_g <- sum(hat_U_g^2)
    max_term_t <- max_w.r.t.tau_g_term(tau_g=tau_g_t, hat_n_g=hat_n_g, d=d,omega2=omega2,
                                       DistSquare_Sum_hat_U_g=DistSquare_Sum_hat_U_g,Sum_DistSquare_hat_U_g=Sum_DistSquare_hat_U_g)
    while (STOP == FALSE){
      tau_g_tp <- tau_g_t*exp(epsilon*tau_g_t*(0.5*d*hat_n_g*(1/tau_g_t-omega2/(tau_g_t*hat_n_g*omega2+1))+
                                         0.5*tau_g_t*omega2*(tau_g_t*hat_n_g*omega2+2)/((tau_g_t*hat_n_g*omega2+1)^2)*DistSquare_Sum_hat_U_g-
                                         0.5*Sum_DistSquare_hat_U_g))
      max_term_tp <- max_w.r.t.tau_g_term(tau_g=tau_g_tp, hat_n_g=hat_n_g, d=d,omega2=omega2,
                                          DistSquare_Sum_hat_U_g=DistSquare_Sum_hat_U_g,Sum_DistSquare_hat_U_g=Sum_DistSquare_hat_U_g)
      if (max_term_tp>=max_term_t){
        if (max_term_tp-max_term_t<=tol){
          STOP <- TRUE
          PICL_hat_U_term <- PICL_hat_U_term + max_term_tp
          tau_list <- c(tau_list,tau_g_tp)
        }else{tau_g_t <- tau_g_tp;max_term_t <- max_term_tp}
      }else{epsilon <- epsilon/2}
    }
  }
  
  # Calculate the 4th integrated log(p(hat_z|K)) term
  PICL_hat_z_term <- sum(lgamma(table_hat_z+delta))-lgamma(N+K*delta)+lgamma(K*delta)-K*lgamma(delta)
  
  # Calculate the constant term
  if (length(a)==1){ # check whether each individual j has different prior a
    PICL_const_term <- (log(Sum_Y)-1)*Sum_Y-sum(log(factorial(Y)))-0.5*log(N*(N-1))+
      N*(-0.5*d*(N-1)*log(2*pi)+a*log(b)-lgamma(a)+lgamma(a+0.5*d*(N-1)))-
      0.5*d*N*log(2*pi)
  }else{
    PICL_const_term <- (log(Sum_Y)-1)*Sum_Y-sum(log(factorial(Y)))-0.5*log(N*(N-1))+
      N*(-0.5*d*(N-1)*log(2*pi))+sum(a*log(b)-lgamma(a)+lgamma(a+0.5*d*(N-1)))-
      0.5*d*N*log(2*pi)
  }
  
  # Return the final PICL value
  return(list(PICL=PICL_Y_term+PICL_hat_V_term+PICL_hat_U_term+PICL_hat_z_term+PICL_const_term,max_tau=tau_list))
}

#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------

################ Pois-LPCM ################

#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------

# Simulation from the Pois-LPCM

Simulation_Directed_PoisLPCM <- function(beta,mu,tau,d=2,z=NULL,Pi=NULL,seed=NULL){
  # beta: LPM intercept parameter
  # mu: K X d matrix with each row denoting a mu_k; MVN center for the latent latent position u
  # tau: 1 X K vector with each entry denoting tau_k; MVN precision for the latent latent position u
  # Note that c(1:N) can be treated as either a N X 1 matrix or a 1 X N matrix or a 1 X N vector
  # d: dimension of latent positions; d=2 by default
  # z: 1 X N vector; the specified memberships; Here, z is possible to contain empty clusters
  # Pi: 1 X K vector; Clustering probability to each group
  
  set.seed(seed)
  
  K <- nrow(mu)
  #----------------------------------------------------------------------------------------------------------------------------------------
  if (!is.null(z)) {
    N <- length(z)
    Z <- t(t(matrix(z,N,K))==(1:K))*1 # transform the membership 1 X N vector z to a N X K matrix Z; empty groups included
  }else{
    Z <- t(rmultinom(N,1,Pi)) # the rmultinom(N,1,Pi) generated matrix is a K X N matrix
    z <- c(Z%*%c(1:K))
  }
  #----------------------------------------------------------------------------------------------------------------------------------------
  library(mvtnorm)
  n_k <- table(factor(z,levels=1:K)) # empty groups might exist
  # Generate covert latent positions u_i for each i
  U <- matrix(0,N,d)
  for (k in 1:K){
    if (n_k[[k]]!=0){U[which(z==k),] <- mvtnorm::rmvnorm(n_k[[k]],mu[k,],(1/tau[k])*diag(d))} # mvtnorm::rmvnorm() outputs a n_k X d matrix
  }
  #----------------------------------------------------------------------------------------------------------------------------------------
  # Calculate pairwise distance^2 for each pair of nodes i,j
  library(Rfast) # for Rfast::Dist()
  dij2 <- Rfast::Dist(U)^2
  #----------------------------------------------------------------------------------------------------------------------------------------
  # Simulate adj Y
  Y <- matrix(rpois(N^2,exp(beta-dij2)),N,N)
  diag(Y) <- 0
  
  return(list(U=U, Z=Z, z=z, Y=Y))
}

#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------

# VB inference from the PoisLPCM

VB_Directed_PoisLPCM <- 
  function(Y,d=2, eta=1,rho2=1,omega2=1,xi=1,psi=1,K=10,delta=1, # delta=rep(1e-3,K), # data, dimension and prior parameters
           U_t=NA,sigma2_t=NA,Pi_t=NA,eta_t=NA,rho2_t=NA,mu_t=NA,omega2_t=NA,xi_t=NA,psi_t=NA,delta_t=NA, # initialize VB parameters
           tol=1e-2,epsilon_init=1,seed=NULL){ # VB algorithm settings
    
    # Y: the N X N observed directed binary network adjacency matrix
    # d: the dimension assumed for latent position space
    # eta,rho2: prior parameters for \beta; default setting N(0,1)
    # omega2: prior parameter for each \bm{\mu}_k; default setting MVN_d(0,1)
    # xi,psi: prior parameters for each \tau_k; default setting Ga(1,1)
    # K: an upper bound of the number of non-empty groups; default setting K=20; 
    # extra groups will be made redundant automatically by very small \tilde{\pi_ik}
    # delta: a 1 X K vector or a single value if a symmetric Dirichlet is considered; prior parameters for \bm{\Pi}; 
    # default setting Dirichlet(rep(1e-3,K)); Shrinkage prior required that each entry of delta should be < 1/K (Rousseau and Mengersen, 2011)
    
    # tol: tolerance level; if the ELBO^{(t+1)}-ELBO^{(t)} < tol, then terminate the implementation; default setting tol=1e-5
    # epsilon_init: initial epsilon set for each the gradient ascent step; default setting epsilon=1
    # seed: the seed for random number generator
    
    set.seed(seed)
    N <- nrow(Y)
    Check.Replace.Smallest.Number <- function(x) max(x,.Machine$double.xmin)
    Check.Replace.Smallest.Number <- Vectorize(Check.Replace.Smallest.Number)
    # It's possible to have very small numbers during the inference and continuing the process would make the corresponding model parameters become too small so that the machine cannot represent and a zero instead appears
    # This will bring the logarithm of such a parameter being -Inf. So we limit the smallest number to .Machine$double.xmin
    
    # U_t: \tilde{U}, a N X d matrix
    # sigma2_t: \tilde{sigma2}, a 1 X N vector
    # Pi_t: \tilde{Pi}, a N X K matrix; each row is the 1 X K vector \tilde{\Pi}_i
    # eta_t and rho2_t: \tilde{eta} and \tilde{rho2}, two normal distribution parameters
    # mu_t: \tilde{mu}, a K X d matrix
    # omega2_t: \tilde{omega2}, a 1 X K vector
    # xi_t: \tilde{xi}, a 1 X K vector
    # psi_t: \tilde{psi}, a 1 X K vector
    # delta_t: \tilde{delta}, a 1 X K vector
    
    # We instead denote ( U_tp,sigma2_tp,Pi_tp,eta_tp,rho2_tp,mu_tp,omega2_tp,xi_tp,psi_tp,delta_tp ) as the new update
    
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    # Initialization
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize covert latent positions U_t by Classical (Metric) Multidimensional Scaling (MDS) on asymmetric distance matrix
    library(ergmito) # for geodesic
    Y_GM <- geodesic(Y)[[1]] # geodesic matrix (GM) of Y, the length/the number of edges of the shortest paths between each pair of nodes
    Y_GM[is.na(Y_GM)] <- 0 # remove NA as 0
    if (is.na(U_t[1])){U_t <- cmdscale(Y_GM,d)} # obtain the initial latent positions by MDS; cmdscale(,d) returns a N X d matrix
    # Note here that MDS may bring the same position for multiple nodes
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize the group MVN distribution center mu_t by first applying kmeans() on the covert latent positions to obtain the initial grouping,
    # and then let each \tilde{mu}_k be the mean of the latent positions contained in group k
    z_init <- kmeans(U_t,centers=K,nstart = 1000)$cluster # check 1000 different start and pick the best one
    if (is.na(mu_t[1])){
      mu_t <- matrix(0,K,d)
      for (k in 1:K){mu_t[k,] <- apply(matrix(U_t[z_init==k,],ncol=d),2,mean)}
    }
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize sigma2_t, a_t and b_t
    if (is.na(sigma2_t[1])){sigma2_t <- rep(1,N)}
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize omega2_t, xi_t, psi_t and delta_t
    if (is.na(omega2_t[1])){omega2_t <- rep(omega2,K)}
    if (is.na(xi_t[1])){xi_t <- rep(xi,K)}
    if (is.na(psi_t[1])){psi_t <- rep(psi,K)}
    if (is.na(delta_t[1])){delta_t <- rep(delta,K)}
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize eta_t and rho2_t
    if (is.na(eta_t[1])){eta_t <- eta}
    if (is.na(rho2_t[1])){rho2_t <- rho2}
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize Pi_t by evaluating the analytical solution based on all other initial settings
    library(proxy) # for proxy::dist; pairwise distance between two different set of vectors
    delta_t0 <- sum(delta_t)
    if (is.na(Pi_t)){
      # We first define the A_{jk} which is a N X K matrix
      Ajk <- t(matrix(0.5*d*(digamma(xi_t)-log(psi_t))+digamma(delta_t)-digamma(delta_t0),K,N))-
        t(matrix(0.5*xi_t/psi_t,K,N))*(proxy::dist(U_t, mu_t, method="Euclidean")^2+(matrix(sigma2_t,N,K)+t(matrix(omega2_t,K,N)))*d)
      exp_Ajk <- exp(Ajk)
      Pi_t <- exp_Ajk/matrix(rowSums(exp_Ajk),N,K) # a N X K matrix with each row being \Pi_i
      Pi_t <- apply(Pi_t,1:2,Check.Replace.Smallest.Number)
    }
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize the epsilon for gradient ascent steps
    epsilon_U <- rep(epsilon_init,N)
    epsilon_beta <- epsilon_Pi <- epsilon_init
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Evaluate the ELBO at the initial state
    # Calculate the 1st sum over i,j = 1,2,...,N and i != j
    DistSquare_u_t_i_minus_u_t_j_NNmatrix <- proxy::dist(U_t, U_t, method="Euclidean")^2 # NNmatrix means N X N matrix
    sigma2_t_i_plus_sigma2_t_j_NNmatrix <- matrix(sigma2_t,N,N)
    sigma2_t_i_plus_sigma2_t_j_NNmatrix <- sigma2_t_i_plus_sigma2_t_j_NNmatrix + t(sigma2_t_i_plus_sigma2_t_j_NNmatrix)
    The_1_plus_2var_t_NNmatrix <- 1+2*sigma2_t_i_plus_sigma2_t_j_NNmatrix
    ELBO_t_Sum_ij <- Y*(eta_t-DistSquare_u_t_i_minus_u_t_j_NNmatrix-d*sigma2_t_i_plus_sigma2_t_j_NNmatrix)-
      exp(eta_t+rho2_t/2-DistSquare_u_t_i_minus_u_t_j_NNmatrix/The_1_plus_2var_t_NNmatrix)/(The_1_plus_2var_t_NNmatrix^(0.5*d))
    diag(ELBO_t_Sum_ij) <- 0
    # Calculate the 2nd sum over i = 1,2,...,N and k = 1,2,...,K
    xi_t_k_NKmatrix <- t(matrix(xi_t,K,N))
    psi_t_k_NKmatrix <- t(matrix(psi_t,K,N))
    ELBO_t_Sum_ik <-
      Pi_t*(0.5*d*(digamma(xi_t_k_NKmatrix)-log(psi_t_k_NKmatrix))-0.5*xi_t_k_NKmatrix/psi_t_k_NKmatrix* # (apply((array(t(U_t),dim=c(d,N,K))-array(do.call(rbind, replicate(N, t(mu_t), simplify=FALSE)),dim=c(d,N,K)))^2,2:3,sum) # array() way to calculate ||u_i-mu_k||^2
              (proxy::dist(U_t, mu_t, method="Euclidean")^2+d*(matrix(sigma2_t,N,K)+t(matrix(omega2_t,K,N))))+
              digamma(t(matrix(delta_t,K,N)))-digamma(delta_t0)-log(Pi_t))
    # Calculate the 3rd sum over k = 1,2,...,K
    ELBO_t_Sum_k <- 
      (delta-delta_t)*(digamma(delta_t)-digamma(delta_t0))+lgamma(delta_t)-0.5/omega2*(rowSums(mu_t^2)+d*omega2_t)+0.5*d*log(omega2_t)+
      (xi-xi_t)*digamma(xi_t)-xi*log(psi_t)-psi*xi_t/psi_t+lgamma(xi_t)+xi_t
    # Calculate the 4th sum over i = 1,2,...,N
    ELBO_t_Sum_i <- 0.5*d*log(sigma2_t)
    # Calculate the ELBO with the rest single terms
    ELBO_t <- sum(ELBO_t_Sum_ij)+sum(ELBO_t_Sum_ik)+sum(ELBO_t_Sum_k)+sum(ELBO_t_Sum_i)-
      0.5/rho2*((eta_t-eta)^2+rho2_t)+0.5*log(rho2_t)-lgamma(delta_t0)
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize the ELBO list
    ELBO_list <- c(ELBO_t)
    #--------------------------------------------------------------------------------------------------------------------------------------------
    # Initialize the update state
    U_tp<-U_t; sigma2_tp<-sigma2_t; Pi_tp<-Pi_t; eta_tp<-eta_t; rho2_tp<-rho2_t
    mu_tp<-mu_t; omega2_tp<-omega2_t; xi_tp<-xi_t; psi_tp<-psi_t; delta_tp<-delta_t
    
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    # Inference steps
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------------------------------------------------------
    t <- 0; STOP <- FALSE
    while (STOP == FALSE){
      t <- t + 1
      if ((t%%100) == 0){
        cat("Iteration:", t,"\n")
      }
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      # Step 1: update U_t and sigma2_t by NGA
      for (i in 1:N){
        Increase_U <- FALSE
        # Define the terms which are calculated for more than twice in order to reduce the computation burden
        # u_t_i_dNminus1matrix <- matrix(U_t[i,],d,N-1) # a d X N-1 matrix with each column being the \tilde{u}_i
        # u_t_i_minus_v_t_ij_dNminus1matrix <- u_t_i_dNminus1matrix-V_t[,i,-i] # matrix(V_t[,i,-i],nrow=d)
        DistSquare_u_t_i_minus_u_t_j_1Nminus1vector <- c(proxy::dist(matrix(U_t[i,],1,2),U_t[-i,],method="Euclidean")^2)
        The_1_plus_2sigma2_t_i_plus_2sigma2_t_j_1Nminus1vector <- 1+2*(sigma2_t[i]+sigma2_t[-i])
        DistSqaure_over_1_plus_2var_t_1Nminus1vector <- DistSquare_u_t_i_minus_u_t_j_1Nminus1vector/The_1_plus_2sigma2_t_i_plus_2sigma2_t_j_1Nminus1vector
        The_exp_term_t_1Nminus1vector <- exp(eta_t+rho2_t/2-DistSqaure_over_1_plus_2var_t_1Nminus1vector)
        The_big_shared_term_t_U_1Nminus1vector <- The_exp_term_t_1Nminus1vector/(The_1_plus_2sigma2_t_i_plus_2sigma2_t_j_1Nminus1vector^(0.5*d+1))
        # v_t_ji_minus_u_t_i_dNminus1matrix <- V_t[,-i,i]-u_t_i_dNminus1matrix
        t_mu_t <- t(mu_t)
        u_t_i_minus_mu_t_k_dKmatrix <- matrix(U_t[i,],d,K)-t_mu_t
        pi_t_ik_times_xi_t_k_over_psi_t_k_1Kvector <- Pi_t[i,]*xi_t/psi_t
        # Calculate the ELBO relevant to the current state U_t[i,] and sigma2_t[i]
        ELBO_t_U <- sum(-(Y[i,-i]+Y[-i,i])*(DistSquare_u_t_i_minus_u_t_j_1Nminus1vector+d*(sigma2_t[i]+sigma2_t[-i]))-2*The_exp_term_t_1Nminus1vector/(The_1_plus_2sigma2_t_i_plus_2sigma2_t_j_1Nminus1vector^(0.5*d)))+
          0.5*d*log(sigma2_t[i]) - 0.5*sum(pi_t_ik_times_xi_t_k_over_psi_t_k_1Kvector*(colSums(u_t_i_minus_mu_t_k_dKmatrix^2)+d*sigma2_t[i]))
        # Set the while loop to check whether the ELBO increases or not
        while (Increase_U == FALSE){
          # Update U_tp
          U_tp[i,] <- U_t[i,]+epsilon_U[i]*sigma2_t[i]*
            (rowSums(t(replicate(d,2*(-Y[i,-i]-Y[-i,i]+2*The_big_shared_term_t_U_1Nminus1vector)))*(matrix(U_t[i,],d,N-1)-t(U_t[-i,])))-
               rowSums(u_t_i_minus_mu_t_k_dKmatrix*t(matrix(replicate(d,pi_t_ik_times_xi_t_k_over_psi_t_k_1Kvector),K,d)))) # matrix(,K,d) here aims to deal with the K=1 case
          # Update sigma2_tp
          sigma2_tp[i] <- sigma2_t[i]*exp(epsilon_U[i]*2/d*sigma2_t[i]*
                                            (0.5*d*(1/sigma2_t[i]-sum(pi_t_ik_times_xi_t_k_over_psi_t_k_1Kvector))-
                                               sum(d*(Y[i,-i]+Y[-i,i])+2*The_big_shared_term_t_U_1Nminus1vector*(2*DistSqaure_over_1_plus_2var_t_1Nminus1vector-d))))
          if (is.nan(sigma2_tp[i])){sigma2_tp[i] <- sigma2_t[i]} # it's possible that some initial state would first bring very very very small sigma2_t leading to NaN; so in this case we lower-bound the value
          # Determine whether the ELBO increases or not; if no increase, half the epsilon and repeat again
          # Calculate the ELBO relevant to the update state U_tp[i,] and sigma2_tp[i]
          # u_tp_i_dNminus1matrix <- matrix(U_tp[i,],d,N-1)
          DistSquare_u_tp_i_minus_u_tp_j_1Nminus1vector <- c(proxy::dist(matrix(U_tp[i,],1,2),U_tp[-i,],method="Euclidean")^2)
          The_1_plus_2sigma2_tp_i_plus_2sigma2_tp_j_1Nminus1vector <- 1+2*(sigma2_tp[i]+sigma2_tp[-i])
          ELBO_tp_U <- sum(-(Y[i,-i]+Y[-i,i])*(DistSquare_u_tp_i_minus_u_tp_j_1Nminus1vector+d*(sigma2_tp[i]+sigma2_tp[-i]))-
                             2*exp(eta_t+rho2_t/2-DistSquare_u_tp_i_minus_u_tp_j_1Nminus1vector/The_1_plus_2sigma2_tp_i_plus_2sigma2_tp_j_1Nminus1vector)/(The_1_plus_2sigma2_tp_i_plus_2sigma2_tp_j_1Nminus1vector^(0.5*d)))+
            0.5*d*log(sigma2_tp[i]) - 0.5*sum(pi_t_ik_times_xi_t_k_over_psi_t_k_1Kvector*(colSums((matrix(U_tp[i,],d,K)-t_mu_t)^2)+d*sigma2_tp[i]))
          # Determine increase or not
          if (ELBO_tp_U>=ELBO_t_U & !is.nan(ELBO_tp_U)){Increase_U<-TRUE}else{epsilon_U[i] <- epsilon_U[i]/2}
        }
        U_t <- U_tp # use the newest state for the inference
        sigma2_t <- sigma2_tp
      }
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      # Step 2: update Pi_t by AS
      Ajk <- t(matrix(0.5*d*(digamma(xi_t)-log(psi_t))+digamma(delta_t)-digamma(delta_t0),K,N))-
        t(matrix(0.5*xi_t/psi_t,K,N))*(proxy::dist(U_tp, mu_t, method="Euclidean")^2+(matrix(sigma2_tp,N,K)+t(matrix(omega2_t,K,N)))*d)
      exp_Ajk <- exp(Ajk)
      Pi_tp <- exp_Ajk/matrix(rowSums(exp_Ajk),N,K)
      Pi_tp <- apply(Pi_tp,1:2,Check.Replace.Smallest.Number)
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      # Step 3: update eta_t and rho2_t by NGA
      Increase_beta <- FALSE
      Sum_Yij <- Y; diag(Sum_Yij) <- 0; Sum_Yij <- sum(Sum_Yij)
      The_1_plus_2sigma2_tp_i_plus_2sigma2_tp_j_NNmatrix <- matrix(sigma2_tp,N,N)
      The_1_plus_2sigma2_tp_i_plus_2sigma2_tp_j_NNmatrix <- 1+2*(The_1_plus_2sigma2_tp_i_plus_2sigma2_tp_j_NNmatrix+
                                                                   t(The_1_plus_2sigma2_tp_i_plus_2sigma2_tp_j_NNmatrix))
      The_1_plus_2var_tp_power0.5d_NNmatrix <- The_1_plus_2sigma2_tp_i_plus_2sigma2_tp_j_NNmatrix^(0.5*d)
      DistSquare_over_1_plus_2var_tp_NNmatrix <- (proxy::dist(U_tp, U_tp, method="Euclidean")^2)/The_1_plus_2sigma2_tp_i_plus_2sigma2_tp_j_NNmatrix
      The_exp_term_t_NNmatrix <- exp(eta_t+rho2_t/2-DistSquare_over_1_plus_2var_tp_NNmatrix)
      The_big_shared_term_t_beta <- The_exp_term_t_NNmatrix/The_1_plus_2var_tp_power0.5d_NNmatrix
      diag(The_big_shared_term_t_beta) <- 0 # since diagonal entries are NAs; recall: diagonal entries of V and varphi are NAs
      Sum_big_shared_term_t_beta <- sum(The_big_shared_term_t_beta)
      # Calculate the ELBO relevant to the current state eta_t and rho2_t
      ELBO_t_beta <- 0.5*(log(rho2_t)-1/rho2*((eta_t-eta)^2+rho2_t))+eta_t*Sum_Yij-Sum_big_shared_term_t_beta
      while (Increase_beta == FALSE){
        # Update eta_tp
        eta_tp <- eta_t+epsilon_beta*rho2_t*(Sum_Yij-(eta_t-eta)/rho2-Sum_big_shared_term_t_beta)
        # Update rho2_tp
        rho2_tp_inside_exp <- epsilon_beta*2*rho2_t*0.5*(1/rho2_t-1/rho2-Sum_big_shared_term_t_beta)
        rho2_tp <- rho2_t*exp(rho2_tp_inside_exp) # so the log(rho2_tp)=log(rho2_t)+rho2_tp_inside_exp in the ELBO; this aims to avoid log(exp(-1000))=-Inf in R
        # Determine whether the ELBO increases or not; if no increase, half the epsilon and repeat again
        # Calculate the ELBO relevant to the update state eta_tp and rho2_tp
        The_exp_over_1_plus_2var_power0.5d_tp_NNmatrix <- exp(eta_tp+rho2_tp/2-DistSquare_over_1_plus_2var_tp_NNmatrix)/The_1_plus_2var_tp_power0.5d_NNmatrix
        diag(The_exp_over_1_plus_2var_power0.5d_tp_NNmatrix) <- 0
        ELBO_tp_beta <- 0.5*((rho2_tp_inside_exp+log(rho2_t))-1/rho2*((eta_tp-eta)^2+rho2_tp))+eta_tp*Sum_Yij-sum(The_exp_over_1_plus_2var_power0.5d_tp_NNmatrix)
        # Determine increase or not
        if (ELBO_tp_beta>=ELBO_t_beta){Increase_beta<-TRUE}else{epsilon_beta <- epsilon_beta/2}
      }
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      # Step 4: update mu_t by AS
      colSums_Pi_tp <- colSums(Pi_tp)
      mu_tp <- matrix(omega2*xi_t/(omega2*xi_t*colSums_Pi_tp+psi_t),K,d)*
        t(apply(array(t(U_tp),dim=c(d,N,K))*array(t(matrix(do.call(rbind, replicate(d, t(Pi_tp), simplify=FALSE)),K,N*d)),dim=c(d,N,K)),c(1,3),sum))
      # t(matrix(do.call(rbind, replicate(d, t(Pi_tp), simplify=FALSE)),K,N*d)) # this aims to construct a N*d X K matrix where, 
      # e.g. when d=2 we have that the 1,1th, 2,1th entries are the same and are Pi_tp[1,1]; similarly, the 3,1th, 4,1th entries are the same and are Pi_tp[2,1]; 
      # the 1,2th, 2,2th entries are the same and are Pi_tp[1,2]; the 3,2th, 4,2th entries are the same and are Pi_tp[2,2], and so on
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      # Step 5: update omega2_t by AS
      omega2_tp <- omega2*psi_t/(psi_t+omega2*xi_t*colSums_Pi_tp)
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      # Step 6: update xi_t,psi_t by AS
      xi_tp <- xi +0.5*d*colSums_Pi_tp
      psi_tp <- psi + 0.5*colSums(Pi_tp*(proxy::dist(U_tp, mu_tp, method="Euclidean")^2+(matrix(sigma2_tp,N,K)+t(matrix(omega2_tp,K,N)))*d))
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      # Step 7: update delta_t by SGA
      Increase_Pi <- FALSE
      # Calculate the ELBO relevant to the current state delta_t
      ELBO_t_Pi <- -lgamma(delta_t0) + sum(Pi_tp*(t(matrix(digamma(delta_t),K,N)))) - digamma(delta_t0)*N + sum((delta-delta_t)*(digamma(delta_t)-digamma(delta_t0))+lgamma(delta_t))
      while (Increase_Pi == FALSE){
        # Update all delta_t simultaneously
        delta_tp <- delta_t*exp(epsilon_Pi*delta_t*(trigamma(delta_t)*(delta-delta_t+colSums_Pi_tp)-trigamma(delta_t0)*(sum(delta-delta_t)+N)))
        delta_tp0 <- sum(delta_tp)
        # Calculate the ELBO relevant to the update state delta_tp
        ELBO_tp_Pi <- -lgamma(delta_tp0) + sum(Pi_tp*(t(matrix(digamma(delta_tp),K,N)))) - digamma(delta_tp0)*N + sum((delta-delta_tp)*(digamma(delta_tp)-digamma(delta_tp0))+lgamma(delta_tp))
        # Determine increase or not
        if (ELBO_tp_Pi>=ELBO_t_Pi & (lgamma(delta_tp0)+1)!=lgamma(delta_tp0)
        ){Increase_Pi<-TRUE}else{epsilon_Pi <- epsilon_Pi/2}
        # (lgamma(delta_tp0)+1)!=lgamma(delta_tp0) case here corresponds to the problem where the the calculations exceeds the capability of R programming, e.g., (lgamma(1.428950e+21)+1)==lgamma(1.428950e+21)
      }
      
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      #--------------------------------------------------------------------------------------------------------------------------------------------
      # Evaluate the ELBO at the update state
      # Calculate the 1st sum over i,j = 1,2,...,N and i != j
      DistSquare_u_tp_i_minus_u_tp_j_NNmatrix <- proxy::dist(U_tp, U_tp, method="Euclidean")^2 # NNmatrix means N X N matrix
      sigma2_tp_i_plus_sigma2_tp_j_NNmatrix <- matrix(sigma2_tp,N,N)
      sigma2_tp_i_plus_sigma2_tp_j_NNmatrix <- sigma2_tp_i_plus_sigma2_tp_j_NNmatrix + t(sigma2_tp_i_plus_sigma2_tp_j_NNmatrix)
      The_1_plus_2var_tp_NNmatrix <- 1+2*sigma2_tp_i_plus_sigma2_tp_j_NNmatrix
      ELBO_tp_Sum_ij <- Y*(eta_tp-DistSquare_u_tp_i_minus_u_tp_j_NNmatrix-d*sigma2_tp_i_plus_sigma2_tp_j_NNmatrix)-
        exp(eta_tp+rho2_tp/2-DistSquare_u_tp_i_minus_u_tp_j_NNmatrix/The_1_plus_2var_tp_NNmatrix)/(The_1_plus_2var_tp_NNmatrix^(0.5*d))
      diag(ELBO_tp_Sum_ij) <- 0
      # Calculate the 2nd sum over i = 1,2,...,N and k = 1,2,...,K
      xi_tp_k_NKmatrix <- t(matrix(xi_tp,K,N))
      psi_tp_k_NKmatrix <- t(matrix(psi_tp,K,N))
      ELBO_tp_Sum_ik <-
        Pi_tp*(0.5*d*(digamma(xi_tp_k_NKmatrix)-log(psi_tp_k_NKmatrix))-0.5*xi_tp_k_NKmatrix/psi_tp_k_NKmatrix* # (apply((array(t(U_tp),dim=c(d,N,K))-array(do.call(rbind, replicate(N, t(mu_tp), simplify=FALSE)),dim=c(d,N,K)))^2,2:3,sum) # array() way to calculate ||u_i-mu_k||^2
                 (proxy::dist(U_tp, mu_tp, method="Euclidean")^2+d*(matrix(sigma2_tp,N,K)+t(matrix(omega2_tp,K,N))))+
                 digamma(t(matrix(delta_tp,K,N)))-digamma(delta_tp0)-log(Pi_tp))
      # Calculate the 3rd sum over k = 1,2,...,K
      ELBO_tp_Sum_k <- 
        (delta-delta_tp)*(digamma(delta_tp)-digamma(delta_tp0))+lgamma(delta_tp)-0.5/omega2*(rowSums(mu_tp^2)+d*omega2_tp)+0.5*d*log(omega2_tp)+
        (xi-xi_tp)*digamma(xi_tp)-xi*log(psi_tp)-psi*xi_tp/psi_tp+lgamma(xi_tp)+xi_tp
      # Calculate the 4th sum over i = 1,2,...,N
      ELBO_tp_Sum_i <- 0.5*d*log(sigma2_tp)
      # Calculate the ELBO with the rest single terms
      ELBO_tp <- sum(ELBO_tp_Sum_ij)+sum(ELBO_tp_Sum_ik)+sum(ELBO_tp_Sum_i)+sum(ELBO_tp_Sum_k)-
        0.5/rho2*((eta_tp-eta)^2+rho2_tp)+0.5*log(rho2_tp)-lgamma(delta_tp0)
      #--------------------------------------------------------------------------------------------------------------------------------------------
      if (ELBO_tp-ELBO_t<=tol){STOP <- TRUE}
      #-----------------------------------------------------------------------------------------------------------------------------------------------
      # Update the current state
      U_t<-U_tp; sigma2_t<-sigma2_tp; Pi_t<-Pi_tp; eta_t<-eta_tp; rho2_t<-rho2_tp
      mu_t<-mu_tp; omega2_t<-omega2_tp; xi_t<-xi_tp; psi_t<-psi_tp; delta_t<-delta_tp; delta_t0<-delta_tp0
      ELBO_list <- c(ELBO_list,ELBO_tp); ELBO_t <- ELBO_tp
    }
    return(list(NumIt=t,ELBO_list=ELBO_list, U_t=U_t,sigma2_t=sigma2_t, Pi_t=Pi_t, eta_t=eta_t,rho2_t=rho2_t,
                mu_t=mu_t,omega2_t=omega2_t, xi_t=xi_t,psi_t=psi_t, delta_t=delta_t,
                epsilon_U=epsilon_U,epsilon_beta=epsilon_beta,epsilon_Pi=epsilon_Pi))
  }

#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------

# The Partially Integrated Classification log-Likelihood (PICL) criteria for the Pois-LPCM model selection

PICL_Directed_PoisLPCM <- function(Y,hat_U,hat_z,K,omega2,delta,tol=.Machine$double.xmin){
  # Y: N X N adjacency matrix
  # hat_U: N X d estimated covert latent positions
  # hat_z: 1 X N estimated clustering
  # Note that hat_z should be label-switched
  N <- nrow(Y); d <- ncol(hat_U); diag(Y) <- 0
  hat_z <- LabelSwitching_SG2003(z=hat_z)$z
  table_hat_z <- table(factor(hat_z,levels=1:K))
  # omega2: prior parameters for \mu_k ~ MVN(0,omega2*I)
  # delta: prior parameters for \Pi ~ Dirichlet(delta,...,delta)
  # tol: the tolerance value to stop the gradient ascent step of max w.r.t. \tau_g
  
  # Calculate the 1st BIC approximation of the log(p(Y|hat_U)) term
  DistSquare_hat_u_i_minus_hat_u_j_NNmatrix <- proxy::dist(hat_U, hat_U, method="Euclidean")^2
  diag(DistSquare_hat_u_i_minus_hat_u_j_NNmatrix) <- 0
  Exp_minus_DistSquare_hat_u_i_minus_hat_u_j_NNmatrix<- exp(-DistSquare_hat_u_i_minus_hat_u_j_NNmatrix)
  diag(Exp_minus_DistSquare_hat_u_i_minus_hat_u_j_NNmatrix) <- 0 # exp(0)=1 in the diagonal
  Sum_Y <- sum(Y)
  
  PICL_Y_term <- -Sum_Y*log(sum(Exp_minus_DistSquare_hat_u_i_minus_hat_u_j_NNmatrix))-sum(Y*DistSquare_hat_u_i_minus_hat_u_j_NNmatrix)
  
  # Calculate the 2nd partial intergal + BIC approximation of the log(p(hat_U|hat_z)) term
  PICL_hat_U_term <- -0.5*K*log(N)
  max_w.r.t.tau_g_term <- function(tau_g,hat_n_g,d,omega2,DistSquare_Sum_hat_U_g,Sum_DistSquare_hat_U_g){ # hat_U_g := {hat_u_i: hat_z_i=g}
    return(0.5*d*(hat_n_g*log(tau_g)-log(tau_g*hat_n_g*omega2+1))+
             0.5*tau_g^2*omega2/(tau_g*hat_n_g*omega2+1)*DistSquare_Sum_hat_U_g-
             0.5*tau_g*Sum_DistSquare_hat_U_g)
  }
  tau_list <- c()
  for (g in 1:K){
    STOP <- FALSE; epsilon <- 1; tau_g_t <- 1; hat_n_g <- table_hat_z[[g]]
    hat_U_g <- matrix(hat_U[hat_z==g,],ncol=2); DistSquare_Sum_hat_U_g <- sum(colSums(hat_U_g)^2); Sum_DistSquare_hat_U_g <- sum(hat_U_g^2)
    max_term_t <- max_w.r.t.tau_g_term(tau_g=tau_g_t, hat_n_g=hat_n_g, d=d,omega2=omega2,
                                       DistSquare_Sum_hat_U_g=DistSquare_Sum_hat_U_g,Sum_DistSquare_hat_U_g=Sum_DistSquare_hat_U_g)
    while (STOP == FALSE){
      tau_g_tp <- tau_g_t*exp(epsilon*tau_g_t*(0.5*d*hat_n_g*(1/tau_g_t-omega2/(tau_g_t*hat_n_g*omega2+1))+
                                                 0.5*tau_g_t*omega2*(tau_g_t*hat_n_g*omega2+2)/((tau_g_t*hat_n_g*omega2+1)^2)*DistSquare_Sum_hat_U_g-
                                                 0.5*Sum_DistSquare_hat_U_g))
      max_term_tp <- max_w.r.t.tau_g_term(tau_g=tau_g_tp, hat_n_g=hat_n_g, d=d,omega2=omega2,
                                          DistSquare_Sum_hat_U_g=DistSquare_Sum_hat_U_g,Sum_DistSquare_hat_U_g=Sum_DistSquare_hat_U_g)
      if (max_term_tp>=max_term_t){
        if (max_term_tp-max_term_t<=tol){
          STOP <- TRUE
          PICL_hat_U_term <- PICL_hat_U_term + max_term_tp
          tau_list <- c(tau_list,tau_g_tp)
        }else{tau_g_t <- tau_g_tp;max_term_t <- max_term_tp}
      }else{epsilon <- epsilon/2}
    }
  }
  
  # Calculate the 3rd integrated log(p(hat_z|K)) term
  PICL_hat_z_term <- sum(lgamma(table_hat_z+delta))-lgamma(N+K*delta)+lgamma(K*delta)-K*lgamma(delta)
  
  # Calculate the constant term
  PICL_const_term <- (log(Sum_Y)-1)*Sum_Y-sum(log(factorial(Y)))-0.5*log(N*(N-1))-
    0.5*d*N*log(2*pi)
  
  # Return the final PICL value
  return(list(PICL=PICL_Y_term+PICL_hat_U_term+PICL_hat_z_term+PICL_const_term,max_tau=tau_list))
}

