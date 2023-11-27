
data_process <- function(sc_exp, sc_label, spot_exp, spot_loc,
                         gene_det_in_min_cells_per = 0.01, expression_threshold = 0,
                         nUMI = 100, verbose = FALSE,depthscale = depthscale,
                         clean.only = TRUE){
  
  #gene by cell matrix
  if(ncol(sc_exp) != length(sc_label))
    stop("Require cell labels!")
  
  if(ncol(spot_exp) != nrow(spot_loc))
    stop("Require x , y coordinations")
  
  #### scRNA-seq data process
  sc_matrix = cleanCounts(sc_exp, gene_det_in_min_cells_per = gene_det_in_min_cells_per,
                          expression_threshold = expression_threshold, nUMI = nUMI,
                          verbose = verbose,depthscale = depthscale,
                          clean.only = clean.only)
  
  
  
  sc_matrix= as.matrix(sc_matrix)
  ind = match(colnames(sc_matrix), colnames(sc_exp))
  sc_label = sc_label[ind]
  
  
  #### SRT data process
  st_matrix = cleanCounts(spot_exp, gene_det_in_min_cells_per = gene_det_in_min_cells_per,
                          expression_threshold = expression_threshold, nUMI = nUMI,
                          verbose = verbose,depthscale = depthscale,
                          clean.only = clean.only)
  
  st_matrix= as.matrix(st_matrix)
  ind_sp = match(colnames(st_matrix), colnames(spot_exp))
  spot_loc = spot_loc[ind_sp, ]
  
  #### find common genes
  com_gene = intersect(rownames(sc_matrix),rownames(st_matrix))
  sc_exp = sc_matrix[com_gene,]
  st_exp = st_matrix[com_gene,]
  
  ### rechecking nUMI
  index_sc <- colSums(sc_exp) >= nUMI
  sc_exp_filter <- sc_exp[,index_sc]
  sc_label_filter <- sc_label[index_sc]
  
  index_st <- colSums(st_exp) >= nUMI
  st_exp_filter = st_exp[,index_st]
  spot_loc_filter <- spot_loc[index_st,]
  
  database <- list(sc_exp = sc_exp_filter, sc_label = sc_label_filter,
                   spot_exp = st_exp_filter, spot_loc = spot_loc_filter)
  return(database)
}


cleanCounts <- function (counts, gene_det_in_min_cells_per = 0.01,
                         expression_threshold = 0 ,
                         nUMI = 100, verbose = FALSE, depthscale = 1,
                         clean.only = FALSE) {
  
  #counts: gene by cells matrix
  
  if (clean.only){
    
    n = ncol(counts)
    
    ##### select of the genes
    filter_index_genes = Matrix::rowSums(counts > expression_threshold) >
      gene_det_in_min_cells_per*n
    
    #### filter the cell
    filter_index_cells = Matrix::colSums(counts[filter_index_genes,] >
                                           expression_threshold) > nUMI
    
    x = counts[ filter_index_genes,filter_index_cells]
    
    if (verbose) {
      message("Resulting matrix has ", ncol(counts), " cells and ", nrow(counts), " genes")
    }
    
    return(x) 
    
  } else {
    
    n = ncol(counts)
    
    ##### select of the genes
    filter_index_genes = Matrix::rowSums(counts > expression_threshold) >
      gene_det_in_min_cells_per*n
    
    #### filter the cell
    filter_index_cells = Matrix::colSums(counts[filter_index_genes,] >
                                           expression_threshold) >nUMI
    
    x = counts[ filter_index_genes,filter_index_cells]
    
    
    sf <- colSums(x)/median(colSums(x))
    
    return(log(sweep(x, 2, sf, '/')*depthscale+1))
    
  }
  
}

### Winsorize expression values to prevent outliers
winsorize <- function (x, qt=.05, both = FALSE) {
  
  if(length(qt) != 1 || qt < 0 ||
     qt > 0.5) {
    stop("bad value for quantile threashold")
  }
  
  lim <- quantile(x, probs=c(qt, 1-qt))
  
  if(both){
    
    x[ x < lim[1] ] <- lim[1]
    x[ x > lim[2] ] <- lim[2]
    
  }else{
    
    x[ x > lim[2] ] <- lim[2] 
    
  }
  return(x)
}


### construct cell-type-specific gene expression data
create_group_exp <- function(sc_exp,sc_label) {
  
  #sc_exp,  single cell gene expression datasets
  #sc_label,  cell annotation of the single cells of the reference
  ##group cells
  
  cell_type = sort(unique(sc_label))
  group = list()
  for(i in 1:length(cell_type)){
    temp_use <- which(sc_label == cell_type[i])
    names(temp_use) <- NULL
    group[[i]] <- temp_use
  }
  sc_group_exp = sapply(group,function(x) Matrix::rowMeans(as.matrix(sc_exp[,x])))
  
  sc_group_exp = as.matrix(sc_group_exp)
  colnames(sc_group_exp) = cell_type
  return(sc_group_exp)
}

#### the similarity matrix matrix based on the spot location information
dis_weight = function(spot_loc = spot_loc ,spot_exp = spot_exp, k = 6,  
                      quantile_prob_bandwidth = 1/3, method = "Hex", 
                      coord_type = "grid"){
  
  if(method == "Hex"){
    
    anndata <- reticulate::import("anndata")
    np <- reticulate::import("numpy")
    sq <- reticulate::import("squidpy")
    
    obsm.meta <- list(spatial = as.matrix(as.data.frame(spot_loc)))
    
    anndata_ST <- anndata$AnnData(X = t(spot_exp) , obsm = obsm.meta)
    
    sq$gr$spatial_neighbors(adata = anndata_ST,
                            spatial_key = 'spatial',
                            coord_type = coord_type,
                            n_neighs = as.integer(k))
    
    mat_adj <- as.matrix(anndata_ST$obsp['spatial_connectivities'])
    
    rownames( mat_adj) <- rownames(spot_loc)
    colnames( mat_adj) <- rownames(spot_loc)
    mat_adj = mat_adj+t(mat_adj)
    diag( mat_adj) = 0
    mat_adj[mat_adj!=0] = 1
    
  }else{
    
    print("Please chose Hex")
  }
  
  #add weight to each edges
  dis_euc  = dist(spot_loc , method = "euclidean")^2
  dis_euc = as.matrix(dis_euc)
  
  weight_network <-  mat_adj * dis_euc
  
  bandwidth_selecting <- apply(weight_network, 2, non_zeros_quantile, prob = quantile_prob_bandwidth)
  
  similarity_network_Gaussian <- exp(-sweep(weight_network, 2, bandwidth_selecting, "/")) * mat_adj
  
  similarity_network_Gaussian_sys <- 0.5 * (similarity_network_Gaussian + t(similarity_network_Gaussian))
  
  
  return(list(dis_weight = similarity_network_Gaussian_sys,
              weight_adj = mat_adj))
}


non_zeros_quantile <- function(x, prob){
  if(sum(x == 0) == length(x)){
    values <- 1
  }else{
    x_non_zero <- x[x > 0]
    values <- quantile(x_non_zero, probs = prob)
  }
  return(values)
}


## reshape the input reference mu and  cell type proportion matrix
reshapemat = function(ref_exp = ref_exp, beta.type = beta.type, cutoff= cutoff){
  
  # ref_exp is a p by K signature matrix, where K  is the number of cell type from reference sc data
  # beta.type is a n by K cell type proportions matrix on captured spots
  
  ind = beta.type > cutoff
  beta.type = beta.type*ind
  beta.type = sweep(beta.type, 1, rowSums( beta.type), '/')
  
  n = nrow(beta.type)
  K = ncol(beta.type)
  
  p = nrow(ref_exp)
  
  #cell type names 
  if(is.null(colnames(beta.type))){
    
    celltype = paste0("celltype", 1:K)
    colnames(beta.type) =  celltype
    
  }else{
    
    celltype = colnames(beta.type)
    
  }
  
  colnames(ref_exp) = celltype
  
  A = matrix(list(),K,1)
  A.adj = B = A
  names(A.adj)= names(A) = names(B) =  celltype
  
  
  #for each cell type 
  for(i in celltype){
    
    A[[i]]  = matrix( rep(beta.type[,i],p), nrow = p, ncol =  n, byrow = TRUE)
    
    B[[i]]  = matrix( rep(ref_exp[,i],n), nrow = p, byrow = FALSE) 
    
    A.adj[[i]]= (A[[i]] > 0)*1
    
  }
  
  #sm <- Matrix(m, sparse = T)
  out = list(A = A, B = B,  A.adj =  A.adj,  beta.type =  beta.type)
  
  return(out)
}


hadamard.sum = function(srt_exp = srt_exp, F_list = F_list,  A_list = A_list, C_mat = C_mat){
  
  K = length(F_list)
  n = ncol(srt_exp)
  
  Qk_list = matrix( list(), K, 1)
  
  for (k in 1:K){
    
    F_list_new = F_list[-k]
    A_list_new = A_list[-k]
    
    hadsum_k = Reduce("+",  Map('*',F_list_new,A_list_new))
    
    Qk = srt_exp- hadsum_k * C_mat
    
    Qk_list[[k]] = Qk
    
  }
  
  return(Qk_list)
}



############ADMM to optimize the problem
ctexp.admm = function(srt_exp = srt_exp, A_list = A_list, A_adj = A_adj, B_list = B_list,
                      L = L, lambda1 = lambda1, lambda2 = lambda2, epsilon = 1e-5,  
                      maxiter = 100, rho = 1,  rho.incr = 1.05, 
                      rho.max = 1e10){
  
  
  timestart<-Sys.time()
  p = nrow(srt_exp)
  n = ncol(L)
  K = length(A_list)
  
  F_list = B_list
  V_list = B_list
  
  P_list = matrix(list(), K, 1)
  for (k in 1:K){
    P_list[[k]] = matrix(0, p, n)
  }
  
  #check the length of tuning parameters of lambdas
  if( length(lambda1)==1 ){  lambda1 = rep( lambda1,K)
  
  }else{ lambda1 = lambda1 }
  
  if( length(lambda2)==1 ){  lambda2 = rep( lambda2,K)
  
  }else{ lambda2 = lambda2 }
  
  
  #The correct factors
  r = matrix( rep(1,p), nrow = p, 1)  
  
  ## the main algorithm
  obj.fun = c()
  
  for (i in 1:maxiter){
    
    ## update S
    recon_est_temp = Reduce("+",  Map('*',F_list, A_list))
    
    S_update_temp = ( matrix(r,p,1) %*% matrix(1,1,n)) * recon_est_temp 
    
    S = sapply(1:n,function (i){ coef(nnls::nnls(as.matrix(S_update_temp[,i]),
                                                 srt_exp[,i]))})
    
    C_mat = (matrix(r,p,1) %*% matrix(S,1,n))
    
    Qk_list = hadamard.sum (srt_exp = srt_exp, F_list = F_list, A_list = A_list, C_mat = C_mat)
    
    #update for each cell type
    F_list.old = F_list
    V_list.old = V_list
    
    reErrorf <- rep(0,K)
    reErrorv <- rep(0,K)
    
    for (k in 1:K){
      
      #update Fk
      F_list[[k]] = update.Fk(Qk = Qk_list[[k]], Ak = A_list[[k]], Bk = B_list[[k]], Vk = V_list[[k]], Pk =  P_list[[k]], 
                              Ak_adj = A_adj[[k]], C_mat = C_mat, rho = rho, lambda2 = lambda2[k])
      
      #update Vk
      V_list[[k]] = update.Vk(Fk = F_list[[k]],  Pk = P_list[[k]], rho = rho, L = L,  
                              Ak_adj = A_adj[[k]],lambda1 = lambda1[k] )
      
      #update Pk
      P_list[[k]] =  P_list[[k]] + rho*(F_list[[k]] - V_list[[k]])
      
      
      reErrorf[k] <- norm(F_list[[k]]-F_list.old[[k]])/(norm(F_list[[k]])+1E-10)
      reErrorv[k] <- norm(V_list[[k]]-V_list.old[[k]])/(norm(V_list[[k]])+1E-10)
    }
    
    
    ## stop criterion
    if(max(c(reErrorf,reErrorv)) < epsilon){
      break
    }
    
    rho = min(rho*rho.incr,rho.max )
    
  }
  
  
  timeend<-Sys.time()
  runningtime <- as.numeric(timeend - timestart) 
  
  recon_est_temp =  Reduce("+",  Map('*',V_list, A_list))
  
  srt_exp_est = ( matrix(r,p,1) %*% matrix(S,1,n)) * recon_est_temp
  
  
  out = list(V.hat = V_list, F.hat = F_list, runningtime =  runningtime,
             srt_exp_est = srt_exp_est)
  
  return(out)
}




### update primal parameters
update.Fk = function(Qk = Qk, Ak = Ak, Bk = Bk, Vk = Vk, Pk = Pk, Ak_adj = Ak_adj,
                     C_mat = C_mat, rho = rho , lambda2 = lambda2  ){
  
  dimn = dim(Ak)
  
  n = ncol( C_mat)
  
  J = matrix(1,dimn[1],dimn[2])
  
  #the Numerators
  nume = rho * Vk - Pk + Qk*Ak*C_mat + 2*lambda2*Bk
  
  #the Denominators 
  deno = Ak *Ak *C_mat*C_mat+  (2*lambda2+rho) * J 
  
  #denoinv =  1/ deno
  updatefk =   nume * (1/ deno)
  
  
  ##########
  ##the non negative constrains 
  
  updatefk = updatefk*Ak_adj
  
  updatefk = updatefk * (updatefk >  1e-5)
  
  
  return(updatefk)
  
}

### update primal parameters
update.Vk = function(Fk = Fk, Pk = Pk, rho = rho, L = L, Ak_adj = Ak_adj, lambda1 =lambda1 ){
  
  nspot = ncol(L)
  
  iden = diag(nspot)
  
  nume =  rho*Fk + Pk
  
  deno = rho*iden + 2* lambda1 *L 
  
  denoinv = MASS::ginv(deno)  
  
  updatevk =   nume %*% denoinv
  
  ##the non negative constrains 
  updatevk = updatevk * Ak_adj
  updatevk = updatevk * (updatevk > 1e-5)
  
  return(updatevk)
  
}

#tuning test selection for the model
tuning.sel.ADMM = function(srt_exp = srt_exp, ref_exp = ref_exp, beta.type = beta.type,
                           L = L,  lambda1 = lambda1, lambda2 = lambda2, cutoff = 0.05, 
                           epsilon = 1e-5, maxiter = 100,  rho = 1,  rho.incr = 1.05, rho.max = 1e10){
  
  p = nrow(srt_exp)
  n = ncol(L)
  K = ncol(beta.type)
  
  #reshape the input data
  reshape.mat = reshapemat(ref_exp = ref_exp, beta.type = beta.type, cutoff = cutoff)
  A_list = reshape.mat$A
  B_list = reshape.mat$B
  A_adj = reshape.mat$A.adj
  beta = reshape.mat$beta.type 
  rm(reshape.mat)
  
  
  ### do not perform batch effect correction
  r = rep(1,p) 
  
  #lambda1 selection
  if(is.null(lambda1)){
    
    cat("We will adpote a value for lambda 1 in our algorithm...", "\n")
    
    #calculate the lambda 1
    recon_est_temp =  Reduce("+",  Map('*', B_list, A_list))
    
    S = sapply(1:n,function (i){ coef(nnls::nnls(as.matrix(recon_est_temp[,i]),
                                                 srt_exp[,i]))})
    ##the final loss
    srt_exp_est = ( matrix(r,p,1) %*% matrix( S,1,n)) * recon_est_temp
    
    obj.loss =  0.5*(norm(srt_exp - srt_exp_est , type = "F"))^2
    
    rm(srt_exp_est)
    rm( recon_est_temp)
    
    
    reg.loss = rep(0,K)
    
    for (k in 1:K){reg.loss[k] <-  sum(diag(t(B_list[[k]]) %*% B_list[[k]] %*% L))}
    
    lambda1 = obj.loss/sum(as.numeric(reg.loss))
    
    
    
  }
  
  ###lambda2 selection
  
  if(is.null(lambda2)){
    
    cat("tuning for lambda 2 in our algorithm...", "\n")
    
    lambda2  = sd(ref_exp)*2
  }
  
  cat("Select value of lambda2", lambda2, "\n")
  
  
  cat("Run the main algorithm...", "\n")
  
  model.final =  ctexp.admm(srt_exp = srt_exp, A_list = A_list, A_adj = A_adj, B_list = B_list,
                            L = L, lambda1 = lambda1, lambda2 = lambda2,
                            epsilon = epsilon, maxiter = maxiter, rho = rho,  rho.incr = rho.incr, rho.max = rho.max)
  
  out = list(V.hat = model.final$V.hat, lambda1 =  lambda1, lambda2 = lambda2, 
             beta = beta, srt_exp_est = model.final$srt_exp_est)
  
  return(out)
  
}

# Main algorithm of the model
#' @param sc_exp scRNA-seq matrix, genes * cells. The format should be raw-counts. The matrix need include gene names and cell names.
#' @param sc_label cell type information. The cell need be divided into multiple categories.
#' @param spot_exp stRNA-seq matrix, genes * spots. The format should be raw-counts. The matrix need include gene names and spot names.
#' @param spot_loc coordinate matrix, spots * coordinates. The matrix need include spot names and coordinate name (x, y).
#' @param beta cell type proportion matrix, spots * cell types. 
#' @param gene_det_in_min_cells_per a floor variable. minimum percent # of genes that need to be detected in a cell.
#' @param expression_threshold a floor variable. Threshold to consider a gene expressed.
#' @param nUMI a floor variable. 	minimum # of read count that need to be detected in a cell or spot.
#' @param verbose a logical variable that defines whether to print the processing flow of data process.
#' @param clean.only a logical variable that defines whether to normalize data with log transform.
#' @param python_env the path of python environment. 
#' @param truncate a logical variable. that defines whether to truncate the gene expression data for both scRNA-seq and SRT.
#' @param qt Winsorize expression values to prevent outliers. Values below this quantile and above 1-this quantile will be set to the quantile value.
#' @param knei number of neighbor spots for construct spatial neighboring graph, details refer to Squidpy.
#' @param methodL the method used to construct spatial neighboring graph, details refer to Squidpy.
#' @param coord_type  grid coordinates, details refer to Squidpy.
#' @param quantile_prob_bandwidth selection of bandwidth of the spatial kernel of each spot.
#' @param lambda1 tuning parameter to balance the graph regularization term. If the tuning parameter is set to NULL, then, we will adopt the value in our algorithm.
#' @param lambda2 tuning parameter to balance the prior regularization term. If the tuning parameter is set to NULL, then, we will adopt the value in our algorithm.
#' @param cutoff  a cutoff value of cell type proportion. 
#' @param rho  the penalty parameter in the ADMM algorithm.
#' @param rho.incr the increase step parameter for varying penalty parameter rho.
#' @param rho.max the maximum value of rho.
#' @param maxiter a positive integer represents the maximum number of updating algorithm. Default setting is 100.
#' @param epsilon a parameter represents the stop criterion.

STged = function(sc_exp = sc_exp, sc_label = sc_label, 
                 spot_exp = spot_exp, spot_loc = spot_loc,   beta = beta,
                 gene_det_in_min_cells_per = 0.01, expression_threshold = 0,
                 nUMI = 100, verbose = FALSE, clean.only = FALSE,python_env = python_env,
                 truncate = TRUE, qt = 0.01,
                 knei = 6,  methodL =  "Hex",coord_type = "grid", quantile_prob_bandwidth = 1/3,
                 lambda1 = NULL, lambda2 = NULL, cutoff = 0.05,  rho = 1,  rho.incr = 1.05,
                 rho.max = 1e10, maxiter = 100,  epsilon = 1e-5){
  
  #Set the environment of Python
  if(is.null(python_env)){
    cat("Please give a Python environment", "\n")
  }
  
  # clear both spatial and scRNA-seq data
  cat("Clear data", "\n")
  datax = data_process(sc_exp = sc_exp, sc_label = sc_label, 
                       spot_exp = spot_exp, spot_loc = spot_loc, depthscale = 1,
                       gene_det_in_min_cells_per = gene_det_in_min_cells_per, 
                       expression_threshold = expression_threshold,
                       nUMI =  nUMI, verbose = verbose, clean.only = clean.only )
  
  #construct spatial correlation
  cat("Construct spatial correlation", "\n")
  L.mat = dis_weight(spot_loc = datax$spot_loc, spot_exp = datax$spot_exp, k = knei, 
                     quantile_prob_bandwidth = quantile_prob_bandwidth, method = methodL, 
                     coord_type = coord_type)
  
  #Set the environment of Python
  if(truncate){
    datax$sc_exp  =  winsorize(x =  datax$sc_exp, qt = qt)
    datax$spot_exp  =  winsorize(x =  datax$spot_exp, qt = qt)
  }
  
  
  #construct reference gene matrix
  cat("Construct reference gene matrix", "\n")
  ref_exp = create_group_exp(sc_exp = datax$sc_exp,sc_label = datax$sc_label)
  
  #the corresponding cell type proportion
  beta = beta[colnames(datax$spot_exp),]
  
  #run the main model
  cat("Run the STged", "\n")
  
  start_time <- Sys.time()
  model.est = tuning.sel.ADMM(srt_exp = datax$spot_exp, ref_exp = ref_exp, 
                              beta.type = beta, L = L.mat$dis_weight, lambda1 =  NULL, lambda2 =  NULL, 
                              cutoff = 0.05, epsilon = 1e-3, 
                              maxiter = 1000, rho = 1,  rho.incr = 1.1, rho.max = 1e10)
  end_time <- Sys.time()
  
  cat("Run time of STged", end_time - start_time,"\n")
  
  
  return(model.est)
  
}

