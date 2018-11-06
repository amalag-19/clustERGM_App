#########################################################################################################
## Loading the required packages
require(lda)
library(Rcpp)
library(RcppArmadillo)
library(foreach)
library(doParallel)
library(ggplot2)
library(abind)
library(reshape)
library(combinat)
require(network)
library(Matrix)
library(igraph)
library(visNetwork)

#########################################################################################################
wrapper_ERGM_stat_undir_Dens<-function(sim.net,nclust,thres=10^(-6),theta_init,sim_indicator,theta_true=NA,K_true=NA,cluster_ids_true=NA){
  ## This is the combined final updated and working code for EM undirected case for all K
  ########################################################################################################
  ########################################################################################################
  ## Loading the required packages
  require(lda)
  library(quadprog)
  library(combinat)
  sourceCpp(file = "rcpparma_hello_world.cpp")

  #################################################
  ## Defining a function to update variational parameters gamma using quadratic program solver. Input:
  ## gamma.curr is a N*K matrix for current gamma estimates, pi.curr is a K*1 vector for current pi
  ## estimates, theta.curr is a T_grid*K array for current theta estimates. Given Network is assumed to be an N*N*T_data array
  gamma.update.wrapper<-function(gamma.curr,pi.curr,theta.curr,network,N,K){
    gamma.next<-matrix(NA_real_,N,K)
    constraint_matrix<-matrix(NA_real_,2+K,K)
    constraint_matrix[1,]<-rep(1,K)
    constraint_matrix[2,]<-rep(-1,K)
    constraint_matrix[3:(K+2),]<-diag(K)
    constraint_vector<-c(1,-1,rep(0,K))

    quad_lin_coeff<-gamma_update_HMM_stat_undir(gamma=gamma.curr, pi=pi.curr, theta=theta.curr, network=network, N=N, K=K)
    for (i in 1:N){
      gamma.next[i,]<-solve.QP((-2*diag(as.vector(quad_lin_coeff[i,,1]))),as.vector(quad_lin_coeff[i,,2]),t(constraint_matrix),constraint_vector)$solution
    }

    ## normalizing gamma_i. deal with case outside (0,1) later
    gamma.next<-t(apply(X = gamma.next,MARGIN = 1,FUN = function(x){
      x_norm<-x/sum(x)
      return(x_norm)
    }))
    return(gamma.next)
  }

  #################################################
  ## Defining a function to update K*1 vector of pi
  pi.update<-function(gamma.curr,N,K){
    pi.next<-as.vector(apply(X = gamma.curr,MARGIN = 2,FUN = mean))
    ## normalization of pi
    pi.next<-pi.next/sum(pi.next)
    return(pi.next)
  }

  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update<-function(theta.curr,pi,gamma,network,N,K){
    theta.next<-rep(NA_real_,K)
    gradient<-grad_HMM_stat_undir(theta=as.vector(theta.curr), gamma=gamma, network=network, N=N, K=K)
    hess<-hess_HMM_stat_undir(theta=as.vector(theta.curr), gamma=gamma, N=N, K=K)
    theta.next<-as.vector(theta.curr)-as.vector(solve(hess)%*%gradient)
    return(theta.next)
  }

  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator<-function(start,network,K,n_iter,thres){
    N<-dim(network)[1]

    ## initializing the arrays for parameters
    gamma<-array(NA_real_,dim=c(N,K,n_iter))
    pi<-matrix(NA_real_,K,n_iter)
    theta<-matrix(NA_real_,K,n_iter)
    gamma[,,1]<-start[[1]]
    pi[,1]<-start[[2]]
    theta[,1]<-start[[3]]
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while((error>thres)&(iter_index<120)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the N*K gamma matrix i.e. variational variational parameters
      gamma[,,iter_index]<-gamma.update.wrapper(gamma.curr=gamma[,,iter_index-1],pi.curr=pi[,iter_index-1], theta.curr=theta[,iter_index-1],network=network,N=N,K=K)
      ## Updating the pi vector
      pi[,iter_index]<-pi.update(gamma.curr=gamma[,,iter_index], N=N, K=K)
      ## Updating the theta matrix
      theta[,iter_index]<-theta.update(theta.curr=theta[,iter_index-1], pi=pi[,iter_index],gamma=gamma[,,iter_index], network=network, N=N, K=K)
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_HMM_stat_undir(gamma=gamma[,,iter_index], pi=pi[,iter_index], theta = theta[,iter_index], network=network, N=N, K=K)
      error<-abs(((ELBO_grid.prev-ELBO_grid.curr)/ELBO_grid.prev))
      print(error)
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(list(gamma,pi,theta))
  }

  ########################################################################################################
  ########################################################################################################
  ## Defining the functions for K=1

  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  ## Defining the functions for K=1

  #################################################
  ## Defining the update function of full matrix theta with dimensions T_grid*K
  theta.update_K1<-function(theta.curr,network,N){
    gradient<-grad_HMM_stat_undir_K1(theta=theta.curr, network=network, N=N)
    hess<-hess_HMM_stat_undir_K1(theta=theta.curr, N=N)
    theta.next<-theta.curr-(gradient/hess)
    return(theta.next)
  }

  #################################################
  ## Defining a function that takes inputs as: 1) initial values list of parameters, 2) network, 3) Number of clusters K and 4) number of iterations and outputs the whole list of all parameter iterates
  iterator_K1<-function(start,network,K,n_iter,thres){
    N<-dim(network)[1]
    ## initializing the arrays for parameters
    theta<-rep(NA_real_,n_iter)
    theta[1]<-start
    ## iterations
    ## Defining the iteration index
    iter_index<-2
    ## Initializing the error
    error<-Inf
    ## Initializing the current ELBO values over the whole grid
    ELBO_grid.curr<-10^10
    while((error>thres)&(iter_index<200)){
      ## Starting the stopwatch to calculate the per iteration time
      ptm<-proc.time()
      ## Updating the theta matrix
      theta[iter_index]<-theta.update_K1(theta.curr=theta[iter_index-1], network=network, N=N)
      ELBO_grid.prev<-ELBO_grid.curr
      ELBO_grid.curr<-ELBO_conv_HMM_stat_undir_K1(theta = theta[iter_index], network=network, N=N)
      error<-sqrt(sum((ELBO_grid.prev-ELBO_grid.curr)^2))/sqrt(sum(ELBO_grid.prev^2))
      print(iter_index)
      print(proc.time()-ptm)
      iter_index<-iter_index+1
    }
    return(theta)
  }

  ########################################################################################################
  ########################################################################################################
  ## Defining Model Selection functions based on converged paramters

  #################################################
  ## Defining a function to calculate the estimate of complete log likelihood
  comp_loglik<-function(gamma=NA,alpha=NA,theta,network,N,K,cluster_ids_est=NA){
    if(K!=1){
      ## 1st term
      t1<-0
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          cluster_id_i<-cluster_ids_est[i]
          cluster_id_j<-cluster_ids_est[j]
          exp_val<-exp(theta[cluster_id_i]+theta[cluster_id_j])
          t1<-t1+((network[i,j]*(theta[cluster_id_i]+theta[cluster_id_j]))-(log(1+exp_val)))
        }
      }
      ## 2nd term
      t2<-0
      for (i in 1:N){
        t2<-t2+log(alpha[cluster_ids_est[i]])
      }
      comp_val<-t1+t2
    }else if(K==1){
      comp_val<-0
      exp_val<-exp(2*theta)
      log_exp_val<-log(1+exp_val)
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          comp_val<-comp_val+((network[i,j]*(2*theta))-log_exp_val)
        }
      }
    }
    return(comp_val)
  }

  #################################################
  ## Defining a function to calculate the integrated classification likelihood
  ICL<-function(gamma=NA,alpha=NA,theta,network,N,K,cluster_ids_est=NA){
    if(K!=1){
      t1<-comp_loglik(gamma = gamma, alpha = alpha, theta = theta,network = network,N = N,K = K, cluster_ids_est = cluster_ids_est)
      t2<-K*log((N*(N-1))/2)
      ICL_val<-t1-t2
    }else if(K==1){
      t1<-comp_loglik(theta = theta,network = network,N = N,K = K)
      t2<-log((N*(N-1))/2)
      ICL_val<-t1-t2
    }
    return(ICL_val)
  }


  ########################################################################################################
  ## Defining the Rand Index and RASE functions for evaluating the performance of clustering and estimation    respectively
  ## Rand Index function
  RI<-function(cluster_ids_est, cluster_ids_true){
    n=length(cluster_ids_est)
    RI_val=0
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        RI_val=RI_val+as.numeric((cluster_ids_est[i]==cluster_ids_est[j])==(cluster_ids_true[i]==cluster_ids_true[j]))
      }
    }
    RI_mean=RI_val/(n*(n-1)/2)
    return(RI_mean)
  }

  #################################################
  ## RASE functions
  RASE_theta<-function(theta_est, theta_true){
    if(is.vector(theta_est)!=1){
      RASE_val<-sqrt(sum((theta_est-theta_true)^2)/dim(theta_est)[1])
    }else{RASE_val<-sqrt(sum((theta_est-theta_true)^2)/length(theta_est))}
    return(RASE_val)
  }

  #########################################################################################################
  ## Defining the parameters
  N<-dim(sim.net)[1] ## Number of nodes from the network
  ## Total time points in the network for which we have data in the form of adjacency matrix
  K<-nclust ## Defining the number of clusters
  ## Setting initial values using package lda which includes mixed membership stochastic block model (MMSB). Using first time point network to run MMSB. Using the mixed membership result of MMSB as our intial gamma. Next we use inital gamma to find initial pi (mixing proportion). Lastly, for network parameter theta we start with 0 matrix.
  set.seed((2))
  MMSB_result <- mmsb.collapsed.gibbs.sampler(beta.prior = list(1,1),K = K,network = sim.net,alpha = 1/2,num.iterations = 100)
  gamma.start <- with(MMSB_result, t(document_sums)/colSums(document_sums))
  ## There is some chance that some component of initial gamma will be exactly 0 which can cause problem in calculating log(gamma). Therfore adding very small amount (10^(-3)) to exact 0 value and rescaling to have sum to 1.
  for(i in 1:N){
    gamma.start[i, which(gamma.start[i,] == 0)] <- rep(10^(-3), length(which(gamma.start[i,] == 0)))
    gamma.start[i,] <- gamma.start[i,]/sum(gamma.start[i,])
  }

  #################################################
  ## Defining the starting values of the iterator and running the main algorithm
  start<-list()

  if(K==1){
    start<-theta_init
    param<-iterator_K1(start=start, network=sim.net, K=K, n_iter=1000, thres=thres)
  }else{
    start[[1]]<-gamma.start
    start[[2]]<-rep(1/K,K) ## alpha (initial distribution)
    start[[3]]<-theta_init ## theta
    #debug(iterator)
    param<-iterator(start=start, network=sim.net, K=K, n_iter=1000, thres=thres)
  }

  #################################################
  ## extracting the coverged parameter values and calculating BIC
  n_iter=1
  indicator_last<-0
  while(indicator_last==0){
    if(K==1){
      temp<-is.na(param[n_iter])
    }else{temp<-is.na(param[[1]][1,1,n_iter])}
    if(temp==TRUE){
      n_last<-n_iter-1
      indicator_last<-1
    }
    n_iter<-n_iter+1
  }
  param_converge<-list()
  if(K==1){
    param_converge<-param[n_last]
    ICL_val<-ICL(theta = param_converge,network = sim.net,N = N,K = K)
    if(sim_indicator==1){
      if(K==K_true){
        RASE_theta<-RASE_theta(theta_est = param_converge,theta_true = theta_true)
        output_list<-list(param_converge,ICL_val,RASE_theta)
      }else{
        output_list<-list(param_converge,ICL_val)
      }
    }else{output_list<-list(param_converge,ICL_val)}
  }else{
    param_converge[[1]]<-param[[1]][,,n_last]
    param_converge[[2]]<-param[[2]][,n_last]
    param_converge[[3]]<-param[[3]][,n_last]
    cluster_ids_est<-as.vector(apply(X = matrix(1:N),MARGIN = 1,FUN = function(x){
      cluster_id<-which.max(param_converge[[1]][x,])
      return(cluster_id)
    }))
    ICL_val<-ICL(gamma = param_converge[[1]], alpha=param_converge[[2]] ,theta = param_converge[[3]],network = sim.net,N = N,K = K, cluster_ids_est = cluster_ids_est)
    if(sim_indicator==1){
      RI_val<-RI(cluster_ids_est = cluster_ids_est,cluster_ids_true = cluster_ids_true)
      if(K==K_true){
        K_permute_mat<-do.call(rbind,permn(1:K))
        RASE_theta_vec<-rep(NA_real_,nrow(K_permute_mat))
        RASE_pi_vec<-rep(NA_real_,nrow(K_permute_mat))
        for (k in 1:nrow(K_permute_mat)){
          theta_est<-param_converge[[3]][K_permute_mat[k,]]
          RASE_theta_vec[k]<-RASE_theta(theta_est = theta_est,theta_true = theta_true)
        }
        permute_true_id_theta<-which.min(RASE_theta_vec)
        RASE_theta_val<-RASE_theta_vec[permute_true_id_theta]
        output_list<-list(param_converge,cluster_ids_est,ICL_val,RI_val,RASE_theta_val)
      }else{
        output_list<-list(param_converge,cluster_ids_est,ICL_val,RI_val)
      }
    }else{output_list<-list(param_converge,cluster_ids_est,ICL_val)}
  }
  return(output_list)
}

#########################################################################################################
net_plotter_igraph<-function(input_demo,adjacency_mat,net_result,nclust){
  igraph_obj<-graph.adjacency(adjacency_mat, mode="undirected")
  set.seed(123456)
  lo <- layout.fruchterman.reingold(igraph_obj)
  V(igraph_obj)$label.cex <- 0.6
  V(igraph_obj)$label.font <- 2
  V(igraph_obj)$Est = net_result[[2]]
  N<-dim(adjacency_mat)[1]
  for(i in 1:N){
    if(V(igraph_obj)$Est[i] == 1){
      V(igraph_obj)$color[i] = "skyblue"
    }
    if(V(igraph_obj)$Est[i] == 2){
      V(igraph_obj)$color[i] = "orange"
    }
    if(V(igraph_obj)$Est[i] == 3){
      V(igraph_obj)$color[i] = "lightgreen"
    }
    if(V(igraph_obj)$Est[i] == 4){
      V(igraph_obj)$color[i] = "lightpink"
    }
  }
  par(mar=c(1,1,1,1), oma=c(1,1,1,1))
  if(input_demo=="yes"){
    label_name <- c("Algeria", "Argentina", "Australia", "Austria", "Barbados", "Belgium", "Bolivia", "Brazil", "Canada", "Chile", "Colombia", "Costa Rica", "Cyprus", "Denmark", "Ecuador", "Egypt", "El Salvador", "Finland", "France", "Germany", "Greece", "Guatemala", "Honduras", "Iceland", "India", "Indonesia", "Ireland", "Israel", "Italy", "Jamaica", "Japan", "South Korea", "Malaysia", "Mauritius", "Mexico", "Morocco", "Nepal", "Netherlands", "New Zealand", "Norway", "Omen", "Panama", "Paraguay", "Peru", "Philippines", "Portugal", "Singapore", "Spain", "Sweden", "Switzerland", "Thailand", "Trinidad and Tobago", "Tunisia", "Turkey", "UK", "USA", "Uruguay", "Venezuela")
    p1<-plot(igraph_obj, vertex.frame.color=NA, vertex.size=7, edge.width=0.5, vertex.label=label_name,layout=lo, main="International Trade Network in 1990")
  }else if(input_demo=="no"){
    p1<-plot(igraph_obj, vertex.frame.color=NA, vertex.size=7, edge.width=0.5, layout=lo, main="Clustered Network")
  }
  return(p1)
}

net_plotter_visNetwork<-function(input_demo,adjacency_mat,net_result,nclust,node_labels=NA){
  igraph_obj<-graph.adjacency(adjacency_mat, mode="undirected")
  df_edges<-data.frame(get.edgelist(graph = igraph_obj))
  names(df_edges)[c(1,2)]<-c("from","to")
  color_vec<-c("skyblue","orange","lightpink","lightgreen")
  N<-dim(adjacency_mat)[1]
  if(input_demo=="yes"){
    label_name <- c("Algeria", "Argentina", "Australia", "Austria", "Barbados", "Belgium", "Bolivia", "Brazil", "Canada", "Chile", "Colombia", "Costa Rica", "Cyprus", "Denmark", "Ecuador", "Egypt", "El Salvador", "Finland", "France", "Germany", "Greece", "Guatemala", "Honduras", "Iceland", "India", "Indonesia", "Ireland", "Israel", "Italy", "Jamaica", "Japan", "South Korea", "Malaysia", "Mauritius", "Mexico", "Morocco", "Nepal", "Netherlands", "New Zealand", "Norway", "Omen", "Panama", "Paraguay", "Peru", "Philippines", "Portugal", "Singapore", "Spain", "Sweden", "Switzerland", "Thailand", "Trinidad and Tobago", "Tunisia", "Turkey", "UK", "USA", "Uruguay", "Venezuela")
    df_nodes <- data.frame(id = 1:N,label = 1:N,group = letters[net_result[[2]]],value = 2,shape = "circle",title = paste0("<p><b>",label_name,"</b></p>"))
    visNet_obj<-visNetwork(df_nodes,df_edges,main = "International Trade Network in 1990", height = "500px",width = "100%") %>% visOptions(highlightNearest = TRUE)%>%visPhysics(enabled = F)%>%visLayout(randomSeed = 1)
    for(k in 1:nclust){
      visNet_obj<-visNet_obj%>%visGroups(groupname = letters[k], color = color_vec[k])
    }
  }else if(input_demo=="no"){
    df_nodes <- data.frame(id = 1:N,label = 1:N,group = letters[net_result[[2]]],value = 1,shape = "circle",title = paste0("<p><b>",node_labels,"</b></p>"))
    visNet_obj<-visNetwork(df_nodes,df_edges,main = "Network Visualization", height = "500px",width = "100%") %>% visOptions(highlightNearest = TRUE)%>%visPhysics(enabled = F)%>%visLayout(randomSeed = 1)
    for(k in 1:nclust){
      visNet_obj<-visNet_obj%>%visGroups(groupname = letters[k], color = color_vec[k])
    }
  }
  return(visNet_obj)
}
#########################################################################################################
load(file = "Trade_data.Rdata")

#########################################################################################################
# df_adj<-read.csv(file = "/Users/Amal/Box Sync/PSU/Spring 2018/Main_Research/Network Models/Apps/clustERGM/example_2.csv")
# 
# net_result<-wrapper_ERGM_stat_undir_Dens(sim.net=adj_mat,nclust = 2,thres=10^(-6),theta_init=jitter(rep(0,2)),sim_indicator=0)
# # 
# debug(net_plotter_visNetwork)
# net_plot<-net_plotter_visNetwork(input_demo = "no",adjacency_mat = adj_mat,net_result = net_result,nclust = 2,node_labels = names(read.csv(file = "/Users/Amal/Box Sync/PSU/Spring 2018/Main_Research/Network Models/Apps/clustERGM/example_2.csv")))
# print(net_plot)

#########################################################################################################
