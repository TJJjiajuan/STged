# STged
R package supporting the paper **"STged: Gene expression deconvolution for spatial transcriptomic data"**. 

STged integrates spatial correlation patterns of gene expression and intra-cell type expression similarity to achieve precise and robust deconvolution results. Implemented within a non-negative least-squares regression framework, STged models gene expression levels at each spot as a weighted linear combination of cell type-specific gene expression, with the weights corresponding to the respective cell type proportions. By incorporating a spatial neighborhood graph prior, STged captures spatial correlation structures in cell type expressions across spots. Moreover, it integrates cell type-specific gene expression information prior from scRNA-seq data to enhance accuracy.
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
data(Fishplus)

##### path on the Windows platform on our computer
python_env <- 'C:/Users/visitor01/.conda/envs/stged/python.exe'

# We run STged as a toy examples
model.est = STged(sc_exp, sc_label, spot_exp, spot_loc, beta,
                  gene_det_in_min_cells_per = 0.01, expression_threshold = 0,
                  nUMI = 100, verbose = FALSE, clean.only = FALSE, depthscale = 1e6,  python_env,
                  truncate = TRUE, qt = 0.0001, knei = 6,  methodL = "Hex",
                  coord_type = "grid", quantile_prob_bandwidth = 1/3,
                  lambda1 = NULL, lambda2 = NULL, cutoff = 0.05,
                  maxiter = 100, epsilon = 1e-5) 

```

## Tutorials
We also give a small example of how to run STged. Here we use the PDAC data ser as an example. There are two ways to run STged, (1) run STged step by step, (2) run the STged by a main function.
- [PDAC with `STged`](https://github.com/TJJjiajuan/STged/blob/main/docs/Demo_STged_FISH.Rmd)
  
Please do not hesitate to contact Dr. Tu at tujiajuan@163.com
to seek any clarifications regarding any content or operation of the
archive.
