![me](https://github.com/Bigardcode/Tripartite_Netwok_tutorial/assets/84800557/f487c73d-7737-49df-9507-da20058a611b)






## Tripartite ceRNA Netwok tutorial


- Step  A. Set directory                             
   Get the current working directory, make sure that it is writable, otherwise, change to a new directory

        getwd()   #shows the directory where R is currently looking for files and saving files to
        setwd("D:/Biard_All_course/ceRNA_Course/ceRNA_Net_Material_Datasets&Codes")   #setwd(working Directory)
        Return file and folder names to console  dir() 
         
- Step  B. CeRNA Network package installation

  if (!require("BiocManager"))  BiocManager::install("BiocManager")

BiocManager::install(c("knitr"))
BiocManager::install(c("SummarizedExperiment"))
BiocManager::install(c("edgeR"))
BiocManager::install(c("reshape2"))
BiocManager::install(c("DESeq2"))
BiocManager::install(c("limma"))
BiocManager::install(c("igraph"))
BiocManager::install(c("dplyr"))
BiocManager::install(c("edgeR"))
BiocManager::install(c("tidyverse"))
BiocManager::install(c("ggalluvial"))
BiocManager::install(c("tidyr"))
BiocManager::install(c("ggsankey"))
BiocManager::install(c("ggplot2"))
BiocManager::install(c("gprofiler2")) 
BiocManager::install(c("tidyr"))
BiocManager::install(c("ggsankey"))
BiocManager::install(c("ggplot2"))
BiocManager::install(c("EnhancedVolcano")) 
BiocManager::install(c("BiocParallel"))
BiocManager::install(c("readr")) 
BiocManager::install(c("data.table"))
BiocManager::install(c("BiocParallel"))
BiocManager::install(c("GEOquery"))
BiocManager::install(c("BiocGeneric"))
BiocManager::install(c("remotes"))


