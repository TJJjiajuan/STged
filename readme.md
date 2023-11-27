# STged
To simultaneously take both spatial correlation structure and gene expression similarity into account, we propose a graph- and prior-guided Gene Expression Deconvolution method on SRT data (STged).
STged is based on a non-negative least-squares regression framework, assuming that gene expression levels of each spot are expressed as a weighted linear combination of cell type-specific 
gene expression and corresponding cell type proportions. The unique features of STged are its ability to account for the spatial correlation structure in cell-type expressions across tissue locations 
through a spatial neighborhood graph and to account for gene expression similarity through integrating prior of cell-type specific gene expression prior information from scRNA-seq data.  As a result, 
STged can exploit both spatial correlation structure and gene expression similarity to achieve accurate and robust gene expression deconvolution of SRT data.
![alt
text](https://github.com/TJJjiajuan/STged/new/main/docs/Figure1.png?raw=true)
