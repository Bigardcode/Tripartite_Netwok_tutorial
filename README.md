![me](https://github.com/Bigardcode/Tripartite_Netwok_tutorial/assets/84800557/f487c73d-7737-49df-9507-da20058a611b)






### Tripartite ceRNA Netwok tutorial


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


- Step C; Get DEGs from CRC  mRNA  datasets                                                      

         Step 1. Retrieval GEO datasets(mRNA datasets) from GEO database   
         Step 2. load datasets into a table
         #Read the data table
         GSE138202 = as.matrix(read.table("GSE138202_mRNA.txt"))

        View(GSE138202)
        dim(GSE138202)
        head(GSE138202, 2)
        tail(GSE138202)
        class(GSE138202)
        mode(GSE138202)


- Step D; Investigating the normalization  of the datasets

          #Check the Normalization
          barplot(GSE138202[1:100,])
          boxplot((GSE138202[1:100,]))

          #log transformation
          logGSE138202 <- log2(GSE138202 + 1)

          barplot(logGSE138202[1:100,])
          boxplot((logGSE138202[1:100,]))

- Step E. Creating design formula for DESeq2                         

        library(DESeq2)
        vignette("DESeq2")

        group = factor(c( rep("Normal", 8), rep("Tumor", 8)))
        head(group)
        class(group)

       ##Create a coldata frame and instantiate the DESeqDataSet. 
       #See ?DESeqDataSetFromMatrix

       colData <- data.frame(group=group, type="paired-end")
       colData
       head(colData)

       #Construct DESEQDataSet Object
       #With the count matrix, cts, and the sample information, coldata, we can construct a DESeqDataSet:
       cds <- DESeqDataSetFromMatrix(GSE138202, colData, design = ~group)

  - Step F. #Normalization with DESeq2
    
        ##set control condition as reference
        cds$group <- relevel(cds$group, ref = "Normal")
        cds <- DESeq(cds)
        resultsNames(cds)
        #results table
         res <- results(cds)

         #let's look at the results table
          head(results(cds, tidy=TRUE)) 
         res = as.data.frame(res)

          #Summary of differential gene expression
           summary(res)
           dim(res)
           class(res)
 - Step 6. Differential expression analysis of mRNA (DEmRNA)          


        #Get DEGs
        dif <- results(cds, pAdjustMethod = "BH", alpha = 0.001)
        dif$padj <- p.adjust(dif$pvalue, method="BH")
        dif <- dif[order(dif$padj),]
        dim(dif)
        mRNAdif = as.data.frame(dif)
        dim(mRNAdif)

        mRNAdown = subset(dif, log2FoldChange < - 2 & padj< 0.01)
        mRNAup = subset(dif,   log2FoldChange > 2 & padj< 0.01)
        nrow(mRNAup)
        nrow(mRNAdown)

        class(mRNAdown)

        mRNAup = as.data.frame(mRNAup)
        mRNAdown= as.data.frame(mRNAdown)

        mRNAUp = rownames(mRNAup)
        mRNADown = rownames(mRNAdown)

- Step 7. volcano plot
  
        res = as.data.frame(res)
        dim(res)

        res1 = res[1:2000, 1:6]
        dim(res1)


        tiff('mRNA.tiff', units="in", width=5.5, height=4.1, res=700, compression = 'lzw')
         EnhancedVolcano(res1, lab = rownames(res1), 
          x = 'log2FoldChange', y = 'padj', titleLabSize = 10, 
          labSize = 1.8,  title = '', subtitle = "", legendPosition = 'bottom', col = c("grey50", "darkred",  "orange", "lawngreen")) 
          dev.off()

- Step 8. Save the DEGs
        export a dataframe or matrix to a csv file

        write.csv(mRNAup , "mRNAnup.csv" , row.names = T)
        write.csv(mRNAdown , "mRNAdown.csv" , row.names = T)
       write.csv(mRNAdif , "mRNAdif.csv" , row.names = T)


## Get DEGs from CRC MiRNA (DEMiRNA) datasets#                    
            #Step 1. set directory
            setwd("D:/BigardCode/ceRNA_course/CRC_ceRNA")
            options(stringsAsFactors=F)
            #setwd()
            getwd() # shows the directory where R is currently looking for files and saving files to
            dir() # You can change the working directory


     #Step 2. Retrieval GEO datasets from GEO database                     
     GSE130084_miRNANA
     #Step 3. #loading gene expression matrix into a table                   
     #MiRNA
     GSE130084 = as.matrix(read.table("GSE130084_miRNA.txt")
     view(GSE130084)
     dim(GSE130084)
     head(GSE130084)
     tail(GSE130084)
     class(GSE130084)
     mode(GSE130084)

## Step 4. #Investigating the normalization  of the data             

      barplot(GSE130084)
      boxplot(GSE130084[1:30,])

   #Log transformation
    logGSE130084 <- log2(GSE130084 + 1)
    barplot(logGSE130084)
    boxplot(logGSE130084)

#Step 5. #Creating the group of sample               
        library(DESeq2)
        group = factor(c( rep("Tumor", 2), rep("Normal", 2)))
        head(group)

##Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix

colData <- data.frame(group=group, type="paired-end")
colData
head(colData)

#Construct DESEQDataSet Object
#With the count matrix, cts, and the sample information, coldata, we can construct a DESeqDataSet:
cds <- DESeqDataSetFromMatrix(GSE130084, colData, design = ~group)

# set control condition as reference
        cds$group <- relevel(cds$group, ref = "Normal")

 
 #Step 5. #Normalization with DESeq2                 
         cds <- DESeq(cds)

##  see all comparisons (here there is only one
           resultsNames(cds)
           res <- results(cds)
           res = as.data.frame(res)
           summary(res)
           dim(res)
           mcols(res, use.names=TRUE)

##  Step 6. Differential expression analysis of mRNA (DEmRNA)           


         summary(results(cds, alpha=0.05))
         dif <- results(cds, pAdjustMethod = "BH", alpha = 0.05)
         dif$padj <- p.adjust(dif$pvalue, method="BH")
         di <- dif[order(dif$padj),]
         dim(dif)
         miRNAdif =  as.data.frame(dif)
         dim(miRNAdif)

        MiRNAdown = subset(miRNAdif, log2FoldChange < - 0.5 & pvalue< 0.05)
        MiRNAup = subset(miRNAdif, log2FoldChange > 0.5 & pvalue< 0.05)
        nrow(MiRNAdown)
        nrow(MiRNAup)


        MiRNAup = as.data.frame(MiRNAup)
        MiRNAdown= as.data.frame(MiRNAdown)

        MiRNAUp = rownames(MiRNAup)
        MiRNADown = rownames(MiRNAdown)


## Step 7. Save the DEGs                                    

##export a dataframe or matrix to a csv file
write.csv(MiRNAup , "MiRNAup.csv" , row.names = T)
write.csv(MiRNAdown , "MiRNAdown.csv", row.names = T)
write.csv(miRNAdif , "miRNAdif.csv", row.names = T)

## Get DEGs from lncRNA CRC  datasets#                    

- Step 1. Set directory                             
setwd("D:/BigardCode/ceRNA_course/CRC_ceRNA")
options(stringsAsFactors=F)
#setwd()
getwd() # shows the directory where R is currently looking for files and saving files to
dir()


## You can change the working directory

- Step 2. Retrieval LncRNA datasets from GEO database                    

## Import the data and look at the first six rows

GSE1048361 <- read.csv(file = 'GSE104836_LncRNA.csv')
GSE104836 = as.matrix(GSE1048361[,-1])
rownames(GSE104836) <- GSE1048361[,1]
GSE104836 = as.data.frame(GSE104836)
dim(GSE104836)
class(GSE104836)
mode(GSE104836)
head(GSE104836)
tail(GSE104836)


GSE104836 = as.numeric(GSE104836)
na.omit(GSE104836)


## Step 3. #Investigating the normalization  of the data                  

barplot(GSE104836[1:100,])
boxplot(GSE104836[1:100,])

#Log transformation
logGSE104836 <- log2(GSE104836 + 1)

barplot(logGSE104836[1:100,])
boxplot(logGSE104836[1:100,])

## Step 4. #Creating the group of samples              

#Creating the group of samples
library(ggplot2)
library(DESeq2)
vignette("DESeq2")
group = factor(c(rep("Tumor", 10), rep("Normal", 10)))
head(group)


## Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix

colData <- data.frame(group=group, type="paired-end")
colData
head(colData)
#Construct DESEQDataSet Object
#With the count matrix, cts, and the sample information, coldata, we can construct a DESeqDataSet:
cds <- DESeqDataSetFromMatrix(GSE104836, colData, design = ~group)

## converting counts to integer mode
cds

## Step 5. #Normalization with Deseq2                      
# set control condition as reference
cds$group <- relevel(cds$group, ref = "Normal")
?DESeqDataSetFromMatrix

cds <- DESeq(cds)
resultsNames(cds)
res <- results(cds)
res = as.data.frame(res)
summary(res)
dim(res)

## Using DEseq2 built in method

normalized_counts <- counts(cds, normalized=T)
cntadolog <- log2(1+counts(cds, normalized=T))

head(normalized_counts)

boxplot(normalized_counts)
boxplot(cntadolog)

## Step 6. Differential expression analysis of mRNA (DEmRNA)          

summary(results(cds, alpha=0.05))

dif <- results(cds, pAdjustMethod = "BH", alpha = 0.05)
dif$padj <- p.adjust(dif$pvalue, method="BH")
dif <- dif[order(dif$padj),]
dim(dif)
lncdif1 =  as.data.frame(dif)
dim(lncdif1)

lncdown = subset(lncdif1, log2FoldChange < - 0.5 & padj< 0.05)
lncup = subset(lncdif1, log2FoldChange > 0.5 & padj< 0.05)
nrow(lncup)
nrow(lncdown)


lncup = as.data.frame(lncdown)
lncdown= as.data.frame(lncdown)

lncUp = rownames(lncup)
lncDown = rownames(lncdown)


## Step 7. Save the DEGs                              

##export a dataframe or matrix to a csv file

write.csv(lncup,"lncup.csv", row.names=T)
write.csv(lncdown,"lncdown.csv", row.names=T)
write.csv(lncdif,"lncdif.csv", row.names=T)

## Exploring interaction between lncRNA-miRNA  & miRNA-mRNA               

- Step 1. loading mRNAlist,mirlist,lnclist                 
setwd("D:/BigardCode/ceRNA_course/CRC_ceRNA")

#mRNAlist
genelist = read.table("genelist.txt")

dim(genelist)
head(genelist)
genelist = genelist[,1]
class(geneList)
head(genelist)

#mirlist
mirlist = read.table("mirlist.txt")
head(mirlist)
dim(mirlist)
mirlist = mirlist[,1]
class(mirlist)
head(mirlist)

#lnclist
lnclist = read.table("lnclist.txt")
head(lnclist)
dim(lnclist)
lnclist = lnclist[,1]
class(lnclist)
head(lnclist)


## Step 2. Retrieval mir-mRNA through  Mirwalk  

setwd("D:/BigardCode/ceRNA_course/CRC_ceRNA")
#load Mirwalk------------------------------------------------------------------------
miRWalk <- read.csv("miRWalk.csv", head = TRUE, sep =",", stringsAsFactors=F)
dim(miRWalk)
head(miRWalk)
miRWalk = as.data.frame(miRWalk)
MirWalk_subset <- subset(miRWalk, subset = genesymbol %in% genelist)
dim(MirWalk_subset)
head(MirWalk_subset)
write.csv(MirWalk_subset,"MirWalk_subset.csv", row.names=T)


## Step 3. Retrieval miRNA-lncRNA through  miRTarBase                           

#load miRTarBase---------------------------------------------------------------------
setwd("D:/BigardCode/ceRNA_course/CRC_ceRNA")
miRTarBase <- read.csv("miRTarBase.csv", head = TRUE, sep =",", stringsAsFactors=F)
dim(miRTarBase)
head(miRTarBase)
miRTarBase = as.data.frame(miRTarBase)
mirtarbase_subset <- subset(miRTarBase, subset = tarName %in% mirlist)
dim(mirtarbase_subset)
head(mirtarbase_subset)

mirtarbase_lncRNA_mirRNA_subset <- subset(mirtarbase_subset, subset = ncName %in% lnclist)
dim(mirtarbase_lncRNA_mirRNA_subset)
head(mirtarbase_lncRNA_mirRNA_subset)

mirtarbase_lncRNA_mRNA_subset <- subset(mirtarbase_subset, subset = tarName %in% genelist)
dim(mirtarbase_lncRNA_mRNA_subset)
head(mirtarbase_lncRNA_mRNA_subset)


write.csv(mirtarbase_subset,"mirtarbase_subset.csv", row.names=T)
write.csv(mirtarbase_lncRNA_mirRNA_subset,"mirtarbase_lncRNA_mirRNA_subset.csv", row.names=T)
write.csv(mirtarbase_lncRNA_mRNA_subset,"mirtarbase_lncRNA_mRNA_subset.csv", row.names=T)


## Step 4. Merge LncRNA-MiRNA & MiRNA                           

#Merge LncRNA-MirRNA & MiRNA----------------------------------------------------------

setwd("D:/BigardCode/ceRNA_course/CRC_ceRNA")

# subset (Mirwalk)
MiRNA_MRNA = MirWalk_subset[,c("mirnaid","genesymbol")]  # returns a data.frame
head(MiRNA_MRNA)
dim(MiRNA_MRNA)
class(MiRNA_MRNA)

# Remove duplicate rows
MiRNA_MRNA <- MiRNA_MRNA[!duplicated(MiRNA_MRNA), ]
head(MiRNA_MRNA)
dim(MiRNA_MRNA)

# subset (mirtarbase)
LncRNA_MiRNA = mirtarbase_lncRNA_mirRNA_subset[,c("ncName","tarName")]  # returns a data.frame
head(LncRNA_MiRNA)
dim(LncRNA_MiRNA)
class(LncRNA_MiRNA)

colnames(LncRNA_MiRNA)[2]  <- "mirnaid"    # change column name for x column

## Step 4. Merge two data frames(LncRNA_MirRNA, LncRNA_MirRNA) by ID                  

LncRNA_MiRNA_MRNA <- merge(MiRNA_MRNA,LncRNA_MiRNA, by="mirnaid")
head(LncRNA_MiRNA_MRNA)
dim(LncRNA_MiRNA_MRNA)
class(LncRNA_MiRNA_MRNA)

LncRNA_MiRNA_MRNA = LncRNA_MiRNA_MRNA %>% relocate(ncName, .before=mirnaid)

#Save the results---------------------------------------------------
write.csv(LncRNA_MiRNA_MRNA,"LncRNA_MiRNA_MRNA.csv", row.names=T)
write.csv(MiRNA_MRNA,"MiRNA_MRNA.csv", row.names=T)
write.csv(LncRNA_MiRNA,"LncRNA_MiRNA.csv", row.names=T)

## Step 5. investagate the ceRNA network interaction                           

view(LncRNA_MiRNA_MRNA)


   
