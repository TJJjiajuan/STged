
#' Data Processing for Spatial and Single-cell Transcriptomics.
#' Processes and cleans spatial transcriptomics (SRT) and single-cell RNA-seq (scRNA-seq)
#' data by filtering genes and cells/spots based on specified criteria. It aligns the
#' gene expression data between scRNA-seq and SRT datasets by retaining only the common genes.
#' @param sc_exp A matrix of scRNA-seq data (genes x cells) containing expression counts.
#' @param sc_label A vector of labels corresponding to the cells in the sc_exp matrix.
#' @param spot_exp A matrix of spatial transcriptomics data (genes x spots) containing expression counts.
#' @param spot_loc A matrix of coordinates for each spot in the spot_exp matrix.
#' @param gene_det_in_min_cells_per Minimum percentage of cells in which a gene must be detected.
#' @param expression_threshold The minimum expression level a gene must exhibit to be considered expressed.
#' @param nUMI Minimum number of UMIs (unique molecular identifiers) required for a cell/spot to be included.
#' @param verbose Logical; if TRUE, detailed processing information will be printed.
#' @param depthscale A scaling factor used for normalization (default is 100,000).
#' @param clean.only Logical; if TRUE, only performs data cleaning without additional normalization.
#' @return A list containing filtered and processed scRNA-seq and SRT data:
#'   - `sc_exp`: Filtered scRNA-seq expression matrix.
#'   - `sc_label`: Corresponding cell labels for the filtered scRNA-seq data.
#'   - `spot_exp`: Filtered SRT expression matrix.
#'
#' @examples
#' data(Fishplus)
#' datax = data_process(sc_exp = sc_exp, sc_label = sc_label, spot_exp = spot_exp, spot_loc = spot_loc)
#' @export
data_process <- function(sc_exp, sc_label, spot_exp, spot_loc,
                         gene_det_in_min_cells_per = 0.01, expression_threshold = 0,
                         nUMI = 100, verbose = FALSE, depthscale = 1e6, clean.only = TRUE) {
  # Ensure correct dimensions of inputs
  if(ncol(sc_exp) != length(sc_label)) {
    stop("Mismatch between number of cells in sc_exp and length of sc_label.")
  }
  if(ncol(spot_exp) != nrow(spot_loc)) {
    stop("Mismatch between number of spots in spot_exp and number of rows in spot_loc.")
  }

  # Process scRNA-seq data
  sc_matrix <- cleanCounts(sc_exp, gene_det_in_min_cells_per = gene_det_in_min_cells_per,
                           expression_threshold = expression_threshold, nUMI = nUMI,
                           verbose = verbose, depthscale = depthscale, clean.only = clean.only)
  sc_matrix <- as.matrix(sc_matrix)
  ind <- match(colnames(sc_matrix), colnames(sc_exp))
  sc_label <- sc_label[ind]

  # Process SRT data
  st_matrix <- cleanCounts(spot_exp, gene_det_in_min_cells_per = gene_det_in_min_cells_per,
                           expression_threshold = expression_threshold, nUMI = nUMI,
                           verbose = verbose, depthscale = depthscale, clean.only = clean.only)
  st_matrix <- as.matrix(st_matrix)
  ind_sp <- match(colnames(st_matrix), colnames(spot_exp))
  spot_loc <- spot_loc[ind_sp, ]

  # Find common genes between datasets
  com_gene <- intersect(rownames(sc_matrix), rownames(st_matrix))
  sc_exp <- sc_matrix[com_gene,]
  st_exp <- st_matrix[com_gene,]

  # Re-check and filter by nUMI counts
  index_sc <- colSums(sc_exp) >= nUMI
  sc_exp_filter <- sc_exp[, index_sc]
  sc_label_filter <- sc_label[index_sc]

  index_st <- colSums(st_exp) >= nUMI
  st_exp_filter <- st_exp[, index_st]
  spot_loc_filter <- spot_loc[index_st,]

  # Compile and return the database of cleaned data
  database <- list(sc_exp = sc_exp_filter, sc_label = sc_label_filter,
                   spot_exp = st_exp_filter, spot_loc = spot_loc_filter)
  return(database)
}


cleanCounts <- function(counts, gene_det_in_min_cells_per = 0.01, expression_threshold = 0,
                        nUMI = 100, verbose = FALSE, depthscale = 1, clean.only = FALSE) {

  if (!inherits(counts, "matrix")) stop("The 'counts' must be a matrix.")
  n = ncol(counts)
  if (verbose) message("Starting with ", n, " cells.")

  # Selecting genes detected in sufficient number of cells
  sufficient_cells = rowSums(counts > expression_threshold) > gene_det_in_min_cells_per * n
  counts = counts[sufficient_cells, ]

  # Filtering cells based on total UMI counts after gene filtering
  sufficient_umi = colSums(counts) > nUMI
  counts = counts[, sufficient_umi]

  if (clean.only) {
    if (verbose) message("Matrix filtered without scaling. Dimension: ", dim(counts))
    return(counts)
  } else {
    sf = colSums(counts) / median(colSums(counts))
    scaled_counts = log(sweep(counts, 2, sf, FUN = "/") * depthscale + 1)
    if (verbose) message("Matrix scaled and log-transformed. Dimension: ", dim(scaled_counts))
    return(scaled_counts)
  }
}



### Winsorize expression values to prevent outliers
winsorize <- function(x, qt = 0.05, both = FALSE) {
  if (length(qt) != 1 || qt < 0 || qt > 0.5) stop("Bad value for quantile threshold")

  lim <- quantile(x, probs = c(qt, 1 - qt))
  if (both) {
    x[x < lim[1]] <- lim[1]
    x[x > lim[2]] <- lim[2]
  } else {
    x[x > lim[2]] <- lim[2]
  }
  return(x)
}


# Function to create cell-type-specific gene expression data
#' Create Grouped Expression Matrix by Cell Type
#'
#' This function computes the mean expression for each cell type from scRNA-seq data.
#' It groups cells by their specified type and calculates the mean of all cells within
#' each type to create a summarized expression profile per cell type.
#' @param sc_exp A matrix of single-cell RNA-seq data (genes x cells), containing expression counts.
#' @param sc_label A vector of labels corresponding to the cells in the sc_exp matrix.
#' @return A matrix of mean expression values for each cell type (genes x cell types).
#' # Assuming `sc_data` is a matrix of gene expression data and `cell_types` is a vector of cell type labels
#' @examples
#' data(Fishplus)
#' refmu = create_group_exp(sc_exp,sc_label)
#' @export
#'
create_group_exp <- function(sc_exp, sc_label) {
  # Identify unique cell types and sort them
  cell_types <- sort(unique(sc_label))
  # Initialize a list to hold indices for each cell type
  group_indices <- vector("list", length(cell_types))

  # Populate the list with indices of cells for each cell type
  for (i in seq_along(cell_types)) {
    group_indices[[i]] <- which(sc_label == cell_types[i])
  }

  # Compute mean expression across cells for each cell type
  group_expression <- sapply(group_indices, function(indices) {
    rowMeans(sc_exp[, indices], na.rm = TRUE)
  }, simplify = "array")

  # Set column names as the cell type names
  colnames(group_expression) <- cell_types

  return(as.matrix(group_expression))
}


# Function to calculate the weighted distance matrix based on spot location
#' Calculate Weighted Adjacency Matrix for Spatial Data
#'
#' This function computes the weighted adjacency matrix for spatial data based on the
#' connectivity and distances between spots. It supports only the 'Hex' method for
#' defining spatial neighbors in a hexagonal grid system.
#'
#' @param spot_loc A matrix of coordinates for each spot, typically named (x, y).
#' @param spot_exp A matrix of spatial transcriptomics data (genes x spots) containing expression counts.
#' @param k Integer; number of neighboring spots to consider when constructing the spatial neighboring graph.
#' @param quantile_prob_bandwidth The bandwidth quantile for weighting the adjacency matrix.
#' @param method Method used for constructing the spatial neighboring graph. Currently supports only 'Hex'.
#' @param coord_type Type of coordinates system used, typically 'grid' for standard spatial transcriptomics data.
#' @return A list containing:
#'   - `dis_weight`: The Gaussian-weighted adjacency matrix based on spatial distances and connectivity.
#'   - `weight_adj`: The unweighted adjacency matrix derived from spatial connectivity alone.
#'
#' @examples
#' # Assuming `coordinates` is a matrix of spatial coordinates and `expression` is the corresponding expression data
#' data(Fishplus)
#' weights <- dis_weight(spot_loc = spot_loc, spot_exp = spot_exp, k = 6, method = "Hex", coord_type = "grid")
#' @export
#'
dis_weight <- function(spot_loc, spot_exp, k = 6, quantile_prob_bandwidth = 1/3, method = "Hex", coord_type = "grid") {
  if (method != "Hex") {
    stop("Unsupported method: Please choose 'Hex'")
  }

  anndata <- reticulate::import("anndata")
  np <- reticulate::import("numpy")
  sq <- reticulate::import("squidpy")

  # Convert location data into a suitable format for Anndata
  obsm.meta <- list(spatial = as.matrix(as.data.frame(spot_loc)))

  # Create an AnnData object using numpy's transpose function for expression data
  anndata_ST <- anndata$AnnData(X = t(spot_exp), obsm = obsm.meta)

  # Use Squidpy to compute spatial neighbors
  sq$gr$spatial_neighbors(adata = anndata_ST, spatial_key = 'spatial', coord_type = coord_type, n_neighs = as.integer(k))

  # Extract and symmetrize the adjacency matrix
  mat_adj <- as.matrix(anndata_ST$obsp[['spatial_connectivities']])


  colnames( mat_adj)  <- rownames( mat_adj) <- rownames(spot_loc)
  mat_adj <- mat_adj + t(mat_adj)
  diag(mat_adj) <- 0

  # Compute squared Euclidean distances between spots
  dis_euc <- as.matrix(dist(spot_loc, method = "euclidean")^2)
  # Apply distances to the adjacency matrix
  weight_network <- mat_adj * dis_euc

  # Compute weights using a quantile of the non-zero distances
  bandwidth_selecting <- apply(weight_network, 2, non_zeros_quantile, prob = quantile_prob_bandwidth)
  similarity_network_Gaussian <- exp(-sweep(weight_network, 2, bandwidth_selecting, "/")) * mat_adj

  # Ensure symmetry in the final weighted matrix
  similarity_network_Gaussian_sys <- 0.5 * (similarity_network_Gaussian + t(similarity_network_Gaussian))

  return(list(dis_weight = similarity_network_Gaussian_sys, weight_adj = mat_adj))
}

# Helper function to find quantiles of non-zero values
non_zeros_quantile <- function(x, prob) {
  x_non_zero <- x[x > 0]
  if(length(x_non_zero) == 0) {
    return(1)
  } else {
    return(quantile(x_non_zero, probs = prob))
  }
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



############ADMM to optimize the problem
MUR = function(srt_exp = srt_exp, A_list = A_list, A_adj = A_adj, B_list = B_list,
               W = W , D = D,
               lambda1 = lambda1, lambda2 = lambda2, epsilon = epsilon,
               maxiter = maxiter){


  timestart<-Sys.time()

  K = length(A_list)

  F_list = B_list
  V_list = B_list

  L = D-W

  #check the length of tuning parameters of lambdas
  if( length(lambda1) == 1 ){ lambda1 = rep( lambda1,  K)

  }else{ lambda1 = lambda1 }

  if( length(lambda2)==1 ){  lambda2 = rep( lambda2, K)

  }else{ lambda2 = lambda2 }


  ### the loss
  srt_exp_est =  Reduce("+",  Map('*', F_list, A_list))
  obj.loss =  (norm(srt_exp - srt_exp_est , type = "F"))^2

  graph.reg.loss =  sapply(seq_along(F_list), function(x) {sum(diag(t(F_list[[x]]*A_adj[[x]]) %*% (F_list[[x]]*A_adj[[x]] ) %*% L))})

  loss.all = obj.loss + sum(graph.reg.loss)

  ## the main algorithm
  obj.fun = c()

  for (i in 1:maxiter){

    #update for each cell type
    F_list.old = F_list
    V_list.old = V_list

    loss.all.old =  loss.all

    Rk_list =  lapply(seq_along(F_list), function(x) F_list[[x]]+ A_list[[x]])

    # update for the F

    for (k in 1:K){

      F_list[[k]] = update.Fk(Fk_old = F_list.old[[k]],srt_exp = srt_exp, Rk = Rk_list[[k]],
                              Ak = A_list[[k]], Bk = B_list[[k]],  Wk = W, Dk = D,
                              lambda1 = lambda1[k], lambda2 = lambda2[k])

    }



    srt_exp_est =  Reduce("+",  Map('*',F_list, A_list))

    obj.loss =  (norm(srt_exp - srt_exp_est , type = "F"))^2


    graph.reg.loss =  sapply(seq_along(F_list), function(x) {sum(diag(t(F_list[[x]]) %*% (F_list[[x]]) %*% L))})

    prior.reg.loss =  sapply(seq_along(F_list), function(x) {
      (norm(F_list[[k]] - B_list[[k]], type = "F"))^2
    })

    loss.all = obj.loss + sum(graph.reg.loss) + sum(prior.reg.loss)


    obj.fun = c(obj.fun , loss.all)

    ## stop criterion

    if(i>5 & ((abs(loss.all - loss.all.old) / abs(loss.all  +loss.all.old))  < epsilon) ){
      break
    }


  }

  timeend<-Sys.time()

  runningtime <- as.numeric(timeend - timestart)

  F.hat =  lapply(seq_along(F_list), function(x) {F_list[[x]]*A_adj[[x]] })


  srt_exp_est =  Reduce("+",  Map('*', F_list, A_list))

  out = list(F_list = F_list, F.hat = F.hat,  srt_exp_est = srt_exp_est,
             runningtime =  runningtime, obj.loss =  obj.fun)

  return(out)
}



### update primal parameters
update.Fk = function(Fk_old = Fk_old, srt_exp = srt_exp, Rk = Rk,
                     Ak = Ak, Bk = Bk,
                     Wk = Wk, Dk = Dk,
                     lambda1 = lambda1 ,lambda2 = lambda2){

  nume  =  srt_exp * Ak + lambda1 * Fk_old %*% Wk + lambda2 * Bk

  #the Denominators
  deno = Rk * Ak + lambda1 * Fk_old %*%Dk + lambda2*Fk_old

  updatefk =  Fk_old * nume * (1/ deno)

  updatefk[!is.finite(updatefk)] = 0

  return(updatefk)

}



#' Multiplicative update rule (MUR) for Spatial Transcriptomics Gene Expression Deconvolution (MUR.STged)
#' This function performs the gene expression deconvolution using a multilevel perceptron approach for
#' spatial transcriptomics data.
#'
#' @param srt_exp A matrix of spatial transcriptomics RNA-seq gene expression data (genes x spots), containing raw counts.
#' @param ref_exp A matrix of reference single-cell RNA-seq gene expression data (genes x cell types).
#' @param beta.type A matrix describing the initial cell type proportions at each spot (spots x cell types).
#' @param W A spatial weight matrix (spots x spots) describing the spatial correlation between spots.
#' @param lambda1 Regularization parameter for the graph regularization term. Automatically determined if set to NULL.
#' @param lambda2 Regularization parameter for the prior regularization term. Automatically determined if set to NULL.
#' @param cutoff Cutoff value for cell type proportion to filter insignificant cell type contributions. Default is 0.05.
#' @param epsilon Convergence threshold for stopping the algorithm. Default is 1e-5.
#' @param maxiter Maximum number of iterations allowed in the algorithm. Default is 100.
#'
#' @return A list containing the following elements:
#'   \item{V.hat}{Estimated gene expression matrix (genes x spots).}
#'   \item{F_list}{List of estimated cell type-specific gene expression matrices.}
#'   \item{lambda1}{Selected value of lambda1.}
#'   \item{lambda2}{Selected value of lambda2.}
#'   \item{beta}{Final estimated cell type proportions matrix (spots x cell types).}
#'   \item{obj.loss}{Objective loss value of the final model.}
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' result <- MUR.STged(srt_exp = spatial_exp, ref_exp = reference_exp, beta.type = beta,
#'                     W = spatial_weights, lambda1 = NULL, lambda2 = NULL, cutoff = 0.05,
#'                     epsilon = 1e-5, maxiter = 100)
#' }
#' @export

MUR.STged = function(srt_exp = srt_exp, ref_exp = ref_exp, beta.type = beta.type,
                     W = W,  lambda1 = lambda1, lambda2 = lambda2, cutoff = 0.05,
                     epsilon = 1e-5, maxiter = 100){


  p = nrow(srt_exp)
  n = ncol(W)
  K = ncol(beta.type)

  #reshape the input data
  reshape.mat = reshapemat(ref_exp = ref_exp, beta.type = beta.type, cutoff = cutoff)
  A_list = reshape.mat$A
  B_list = reshape.mat$B
  A_adj = reshape.mat$A.adj
  beta = reshape.mat$beta.type
  rm(reshape.mat)


  D  = diag(apply(W,1,sum))

  L = D - W

  #lambda1 selection
  if(is.null(lambda1)){

    cat("We will adpote a value for lambda 1 in our algorithm...", "\n")

    #calculate the lambda 1
    srt_exp_est = Reduce("+",  Map('*', B_list, A_list))

    obj.loss =  (norm(srt_exp - srt_exp_est , type = "F"))^2

    rm(srt_exp_est)


    reg.loss = rep(0,K)

    for (k in 1:K){reg.loss[k] <-  sum(diag(t(B_list[[k]]*A_adj[[k]]) %*% (B_list[[k]]*A_adj[[k]] )%*% L))}

    lambda1 = 0.5*obj.loss/sum(as.numeric(reg.loss))

  }
  cat("Select value of lambda1", lambda1, "\n")

  ###lambda2 selection

  if(is.null(lambda2)){

    cat("tuning for lambda 2 in our algorithm...", "\n")

    lambda2  = sd(ref_exp)
  }

  cat("Select value of lambda2", lambda2, "\n")


  cat("Run the main algorithm...", "\n")

  model.final = MUR(srt_exp = srt_exp, A_list = A_list, A_adj = A_adj, B_list = B_list,
                     W = W , D = D,
                     lambda1 = lambda1, lambda2 = lambda2, epsilon = epsilon,
                     maxiter = maxiter)


  out = list(V.hat = model.final$F.hat, F_list = model.final$F_list , ambda1 =  lambda1, lambda2 = lambda2,
             beta = beta, obj.loss =  model.final $obj.loss)

  return(out)

}


#' Spatial Transcriptomics Gene Expression Deconvolution (STged)
#'
#' @param sc_exp A matrix of single-cell RNA-seq gene expression data (genes x cells), containing raw counts.
#' @param sc_label A vector or factor of cell type labels corresponding to each cell in the sc_exp data.
#' @param spot_exp A matrix of spatial transcriptomics RNA-seq gene expression data (genes x spots), containing raw counts.
#' @param spot_loc A matrix specifying the coordinates of each spot (spots x coordinates), typically named (x, y).
#' @param beta A matrix describing the cell type proportions at each spot (spots x cell types).
#' @param gene_det_in_min_cells_per Minimum percentage of genes that must be detected in a cell, used as a filtering criterion.
#' @param expression_threshold Minimum expression level to consider a gene as expressed, used for data cleaning.
#' @param nUMI Minimum number of read counts required to consider a cell or spot, used for data cleaning.
#' @param verbose Logical; if TRUE, prints details of the data processing steps.
#' @param clean.only Logical; if TRUE, applies only data cleaning steps without normalization.
#' @param depthscale A scaling factor used for normalization (default is 100,000).
#' @param python_env Path to the Python environment, used for calling Python functions.
#' @param truncate Logical; if TRUE, truncates gene expression data to reduce the impact of extreme values.
#' @param qt Quantile used for Winsorizing expression values to control the effects of outliers.
#' @param knei Integer; number of neighboring spots to consider when constructing the spatial neighboring graph.
#' @param methodL Method used for constructing the spatial neighboring graph, details provided by the Squidpy package.
#' @param coord_type Type of coordinates system used, e.g., 'grid' as described by Squidpy.
#' @param quantile_prob_bandwidth Bandwidth for the spatial kernel used in constructing the neighboring graph.
#' @param lambda1 Tuning parameter for balancing the graph regularization term. Automatically determined if set to NULL.
#' @param lambda2 Tuning parameter for balancing the prior regularization term. Automatically determined if set to NULL.
#' @param cutoff Cutoff value for cell type proportion to filter insignificant cell type contributions.
#' @param rho Penalty parameter in the ADMM algorithm, with adaptive adjustments.
#' @param rho.incr Increment step for varying the penalty parameter rho.
#' @param rho.max Maximum allowable value for rho.
#' @param maxiter Maximum number of iterations allowed in the algorithm.
#' @param epsilon Convergence threshold for stopping the algorithm.
#' @return A list containing model estimation results, including estimated parameters and diagnostics.
##'@usage model_results <- STged(sc_exp = sc_exp, sc_label = sc_label,
#'                       spot_exp = spot_exp, spot_loc = spot_loc, beta = beta)
#' @export
STged <- function(sc_exp, sc_label, spot_exp, spot_loc, beta,
                  gene_det_in_min_cells_per = 0.01, expression_threshold = 0,
                  nUMI = 100, verbose = FALSE, clean.only = FALSE, depthscale = 1e6,python_env,
                  truncate = TRUE, qt = 0.0001, knei = 6, methodL = "Hex",
                  coord_type = "grid", quantile_prob_bandwidth = 1/3,
                  lambda1 = NULL, lambda2 = NULL, cutoff = 0.05,
                  maxiter = 100, epsilon = 1e-5) {

  if(is.null(python_env)) {
    stop("Python environment path is not specified. Please provide a valid path.")
  }

  if(verbose) cat("Data cleaning and preprocessing...\n")
  datax <- data_process(sc_exp = sc_exp, sc_label = sc_label, spot_exp = spot_exp, spot_loc = spot_loc,
                        gene_det_in_min_cells_per = gene_det_in_min_cells_per, expression_threshold = expression_threshold,
                        nUMI = nUMI ,verbose = verbose,
                        depthscale = depthscale,clean.only = clean.only)
  if(truncate) {
    datax$sc_exp <- winsorize(datax$sc_exp, qt)
    datax$spot_exp <- winsorize(datax$spot_exp, qt)
  }

  if(verbose) cat("Constructing spatial correlation matrix...\n")
  L.mat <- dis_weight(spot_loc = datax$spot_loc, spot_exp = datax$spot_exp,
                      k = knei, quantile_prob_bandwidth = quantile_prob_bandwidth,
                      method = methodL, coord_type = coord_type)


  if(verbose) cat("Constructing reference gene expression matrix...\n")
  ref_exp <- create_group_exp(sc_exp = datax$sc_exp, sc_label = sc_label)

  beta <- beta[colnames(datax$spot_exp), ]

  if(verbose) cat("Running the STged model...\n")
  start_time <- Sys.time()
  model.est <- MUR.STged(srt_exp  = datax$spot_exp, ref_exp = ref_exp, beta.type = beta,
                         w = L.mat$dis_weight,
                         lambda1 =lambda1, lambda2= lambda2, cutoff = cutoff,
                         epsilon =epsilon, maxiter = maxiter)
  end_time <- Sys.time()


  if(verbose) cat("Total run time: ", end_time - start_time, " seconds.\n")

  return(model.est)
}


