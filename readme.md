# STged
R package supporting the paper **"STged: Gene expression deconvolution for spatial transcriptomic data"**. 

To simultaneously take both spatial correlation structure and gene expression similarity into account, we propose a graph- and prior-guided Gene Expression Deconvolution method on SRT data (STged).
STged is based on a non-negative least-squares regression framework, assuming that gene expression levels of each spot are expressed as a weighted linear combination of cell type-specific 
gene expression and corresponding cell type proportions. The unique features of STged are its ability to account for the spatial correlation structure in cell-type expressions across tissue locations 
through a spatial neighborhood graph and to account for gene expression similarity through integrating prior of cell-type specific gene expression prior information from scRNA-seq data.  As a result, 
STged can exploit both spatial correlation structure and gene expression similarity to achieve accurate and robust gene expression deconvolution of SRT data.

### Overview of STged
![alt
text](https://github.com/TJJjiajuan/STged/blob/main/docs/STged_mian.PNG?raw=true)

## Install STged
In the STged, there are several Python modules and R packages are needed to install first.
-   **Install Python dependencies**
``` buildoutcfg
pip install squidpy
pip install numpy
pip install anndata
```

-   **Install R dependencies**
``` buildoutcfg
install.packages("Matrix", repos="http://R-Forge.R-project.org")
install.packages("nnls")
install.packages("MASS")
install.packages('reticulate')
```

-   **Install STged**
``` buildoutcfg
devtools::install_github("TJJjiajuan/STged")
```

## Run the example
``` buildoutcfg
#load the input data set
data("PDAC")

##### path on the Windows platform on our computer
python_env <- 'C:/Users/visitor01/.conda/envs/stged/python.exe'

# We run STged as a toy examples
model.est = STged(sc_exp = sc_exp, sc_label = sc_label, 
                  spot_exp = spot_exp, spot_loc = spot_loc,  beta = beta,
                  gene_det_in_min_cells_per = 0.01, expression_threshold = 0,
                  nUMI = 100, verbose = FALSE, clean.only = FALSE, python_env = python_env,
                  truncate = TRUE, qt = 0.01,
                  knei = 6,  methodL =  "Hex",coord_type = "grid", quantile_prob_bandwidth = 1/3,
                  lambda1 = NULL, lambda2 = NULL, cutoff = 0.05,  rho = 1,  rho.incr = 1.05,
                  rho.max = 1e10, maxiter = 20,  epsilon = 1e-5)

```

## Tutorials
We also give a small example of how to run STged. Here we use the PDAC data ser as an example. There are two ways to run STged, (1) run STged step by step, (2) run the STged by a main function.
- [PDAC with `STged`](https://github.com/TJJjiajuan/STged/blob/main/docs/Demo_STged_PDAC.html)
  
Please do not hesitate to contact Dr. Tu at tujiajuan@163.com
to seek any clarifications regarding any content or operation of the
archive.
