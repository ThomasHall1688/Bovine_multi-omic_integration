########################################################################################
##########################           DESeq2            #################################
########################################################################################

#LANGUAGE : R
#This is the peliminary pipeline for normalisation. The input data will vary, 
#so writing the table will also vary in the following R code. This will mainly be 
#in R

#Libraries needed
library(stringr)
library(DESeq2)
library(tidyverse)
library(AnnotationFuncs)  
library(Biobase)          
library(plyr)             
library(dplyr)            
library(edgeR)            
library(extrafont)        
library(ggplot2)          
library(GO.db)            
library(goseq)            
library(grDevices)        
library(grid)             
library(gridExtra)         
library(limma)            
library(magrittr)         
library(MASS)             
library(org.Bt.eg.db)        
library(svglite)          
library(tools)            
library(VennDiagram)     
library(pheatmap) 
library(RColorBrewer)
library(BiocParallel)
library(biomaRt)

#####Load in data#####
countData <- read.csv(".csv",sep=",", row.names=1)
colData <- read.csv("Alv_mac_colData.csv", row.names=1)

#####clean rownames#####
rownames(countData) %<>%
    str_replace("BGD.*,", "") %>%
    str_replace(",miRBase:.*", "") %>% 
    str_replace("GeneID:", "")
	
##########note##########
#for the Alv_mac data I had to redo the count data construction as for some reason read count and gene name were misaligned. Additionally, FeatureCounts
#will return the command used to generate the file in the first line and least the directory of the sample name. 
#This is the code I used to delete first line and cat all count data colums for each CN, MB and TB count directory (this should be done from within the directory)
#havent found a way to remove the directory listing as the file name. 
awk 'FNR==1{f++}{a[f,FNR]=$7}END{for(x=1;x<=FNR;x++){for(y=1;y<ARGC;y++)printf("%s ",a[y,x]);print ""}}' *.txt > Alv_mac_new_counts.csv 
sed -i '1d' Alv_mac_new_counts.csv
#This can be used to cat all of the data into one file (raw counts for instance)


########## Construct DESeqDataSet ###########


dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ TT)
							  
#Assign factor level with “control” or “untreated”as “reference level”. In this case the reference level was set to the negative samples.
#log2 fold changes are compared against reference level logs 
#However, for the TT design format, I do not use an rlevel and instead specify groups in the contrast section. 
#dds$condition <- relevel(dds$condition, "untreated")

dds <- DESeq(dds)
dds

#result summary of the normalisation. Also order the results based on criteria as well as selcting a subset of those genes for visualization 
res <- results(dds)
res
summary(res)

#####multiple contrasts for timepoint and type (TT) specific DE results 

MB2_res <- results(dds, contrast=c("TT","MB2","CN2")) #res <- results(dds, contrast=c("condition","treated","untreated"))
MB6_res <- results(dds, contrast=c("TT","MB6","CN6"))
MB24_res <- results(dds, contrast=c("TT","MB24","CN24"))
MB48_res <- results(dds, contrast=c("TT","MB48","CN48"))
TB2_res <- results(dds, contrast=c("TT","TB2","CN2"))
TB6_res <- results(dds, contrast=c("TT","TB6","CN6"))
TB24_res <- results(dds, contrast=c("TT","TB24","CN24"))
TB48_res <- results(dds, contrast=c("TT","TB48","CN48"))

# To extract the normalised reads as a .CSV
# After you run the dds <- DESeq(dds)
# Calculate and save the annotated normalised counts

normalized.counts <- as.data.frame(counts( dds, normalized=TRUE))
write.csv(normalized.counts, file="Alv_mac_normalized_counts_ARD.csv")

#write results to .csv
write.csv(as.data.frame(MB2_res), file="MB_2hpi_DE_genes.csv")
write.csv(as.data.frame(MB6_res), file="MB_6hpi_DE_genes.csv")
write.csv(as.data.frame(MB24_res), file="MB_24hpi_DE_genes.csv")
write.csv(as.data.frame(MB48_res), file="MB_48hpi_DE_genes.csv")

write.csv(as.data.frame(TB2_res), file="TB_2hpi_DE_genes.csv")
write.csv(as.data.frame(TB6_res), file="TB_6hpi_DE_genes.csv")
write.csv(as.data.frame(TB24_res), file="TB_24hpi_DE_genes.csv")
write.csv(as.data.frame(TB48_res), file="TB_48hpi_DE_genes.csv")
############ Data transformation ##########

rld <- rlog(dds, blind=FALSE)
vst <- vst(dds, blind=FALSE)

MB2_rld <- rld[ , rld$TT %in% c("MB2", "CN2") ]
MB2_vst <- vst[ , vst$TT %in% c("MB2", "CN2") ]
MB6_rld <- rld[ , rld$TT %in% c("MB6", "CN6") ]
MB6_vst <- vst[ , vst$TT %in% c("MB6", "CN6") ]
MB24_rld <- rld[ , rld$TT %in% c("MB24", "CN24") ]
MB24_vst <- vst[ , vst$TT %in% c("MB24", "CN24") ]
MB48_rld <- rld[ , rld$TT %in% c("MB48", "CN48") ]
MB48_vst <- vst[ , vst$TT %in% c("MB48", "CN48") ]

TB2_rld <- rld[ , rld$TT %in% c("TB2", "CN2") ]
TB2_vst <- vst[ , vst$TT %in% c("TB2", "CN2") ]
TB6_rld <- rld[ , rld$TT %in% c("TB6", "CN6") ]
TB6_vst <- vst[ , vst$TT %in% c("TB6", "CN6") ]
TB24_rld <- rld[ , rld$TT %in% c("TB24", "CN24") ]
TB24_vst <- vst[ , vst$TT %in% c("TB24", "CN24") ]
TB48_rld <- rld[ , rld$TT %in% c("TB48", "CN48") ]
TB48_vst <- vst[ , vst$TT %in% c("TB48", "CN48") ]

########################################################################################
##################################  Data analysis  #####################################
########################################################################################

# res as a data frame for dow

MB2 <- results(dds, contrast=c("TT","MB2","CN2"), tidy = TRUE)
MB6 <- results(dds, contrast=c("TT","MB6","CN6"), tidy = TRUE)
MB24 <- results(dds, contrast=c("TT","MB24","CN24"), tidy = TRUE)
MB48 <- results(dds, contrast=c("TT","MB48","CN48"), tidy = TRUE)
TB2 <- results(dds, contrast=c("TT","TB2","CN2"), tidy = TRUE)
TB6 <- results(dds, contrast=c("TT","TB6","CN6"), tidy = TRUE)
TB24 <- results(dds, contrast=c("TT","TB24","CN24"), tidy = TRUE)
TB48 <- results(dds, contrast=c("TT","TB48","CN48"), tidy = TRUE)

#Table for DE summary
#MB24
Table_MB24_FDR_0.05_l2fc1 <- subset(MB24_0.05_Full, abs(log2FoldChange) > 1)
Table_MB24_FDR_0.01 <- subset(MB24_0.05_Full, padj <= 0.01)
Table_MB24_FDR_0.01_l2fc2 <- subset(Table_MB24_FDR_0.01, abs(log2FoldChange) > 2)
Table_MB24_FDR_0.000001 <- subset(MB24_0.05_Full, padj <= 0.000001)
Table_MB24_FDR_0.000001_l2fc2 <- subset(Table_MB24_FDR_0.000001, abs(log2FoldChange) > 2)

Table_MB24_FDR_0.05 <- MB24_0.05_Full
Table_MB24_FDR_0.05_0p <- subset(MB24_0.05_Full, log2FoldChange > 0)
Table_MB24_FDR_0.05_0n <- subset(MB24_0.05_Full, log2FoldChange < 0)
Table_MB24_FDR_0.05_1p <- subset(MB24_0.05_Full, log2FoldChange > 1)
Table_MB24_FDR_0.05_1n <- subset(MB24_0.05_Full, log2FoldChange < -1)
Table_MB24_FDR_0.01_0p <- subset(Table_MB24_FDR_0.01, log2FoldChange > 0)
Table_MB24_FDR_0.01_0n <- subset(Table_MB24_FDR_0.01, log2FoldChange < 0)
Table_MB24_FDR_0.01_1p <- subset(Table_MB24_FDR_0.01, log2FoldChange > 1)
Table_MB24_FDR_0.01_1n <- subset(Table_MB24_FDR_0.01, log2FoldChange < -1)

#TB24
TB24_0.05_Full <- TB24_FDR_0.05 
Table_TB24_FDR_0.01 <- subset(TB24_0.05_Full, padj <= 0.01)
Table_TB24_FDR_0.05 <- TB24_0.05_Full
Table_TB24_FDR_0.05_0p <- subset(TB24_0.05_Full, log2FoldChange > 0)
Table_TB24_FDR_0.05_0n <- subset(TB24_0.05_Full, log2FoldChange < 0)
Table_TB24_FDR_0.05_1p <- subset(TB24_0.05_Full, log2FoldChange > 1)
Table_TB24_FDR_0.05_1n <- subset(TB24_0.05_Full, log2FoldChange < -1)
Table_TB24_FDR_0.01_0p <- subset(Table_TB24_FDR_0.01, log2FoldChange > 0)
Table_TB24_FDR_0.01_0n <- subset(Table_TB24_FDR_0.01, log2FoldChange < 0)
Table_TB24_FDR_0.01_1p <- subset(Table_TB24_FDR_0.01, log2FoldChange > 1)
Table_TB24_FDR_0.01_1n <- subset(Table_TB24_FDR_0.01, log2FoldChange < -1)


#MB6
Table_MB6_FDR_0.05_l2fc1 <- subset(MB6_FDR_0.05, abs(log2FoldChange) > 1)
Table_MB6_FDR_0.01 <- subset(MB6_FDR_0.05, padj <= 0.01)
Table_MB6_FDR_0.01_l2fc2 <- subset(Table_MB6_FDR_0.01, abs(log2FoldChange) > 2)
Table_MB6_FDR_0.000001 <- subset(MB6_FDR_0.05, padj <= 0.000001)
Table_MB6_FDR_0.000001_l2fc2 <- subset(Table_MB6_FDR_0.000001, abs(log2FoldChange) > 2)

Table_MB6_FDR_0.05 <- MB6_FDR_0.05
Table_MB6_FDR_0.05_0p <- subset(MB6_FDR_0.05, log2FoldChange > 0)
Table_MB6_FDR_0.05_0n <- subset(MB6_FDR_0.05, log2FoldChange < 0)
Table_MB6_FDR_0.05_1p <- subset(MB6_FDR_0.05, log2FoldChange > 1)
Table_MB6_FDR_0.05_1n <- subset(MB6_FDR_0.05, log2FoldChange < -1)
Table_MB6_FDR_0.01_0p <- subset(Table_MB6_FDR_0.01, log2FoldChange > 0)
Table_MB6_FDR_0.01_0n <- subset(Table_MB6_FDR_0.01, log2FoldChange < 0)
Table_MB6_FDR_0.01_1p <- subset(Table_MB6_FDR_0.01, log2FoldChange > 1)
Table_MB6_FDR_0.01_1n <- subset(Table_MB6_FDR_0.01, log2FoldChange < -1)


#MB48
Table_MB48_FDR_0.05_l2fc1 <- subset(MB48_0.05_Full, abs(log2FoldChange) > 1)
Table_MB48_FDR_0.01 <- subset(MB48_0.05_Full, padj <= 0.01)
Table_MB48_FDR_0.01_l2fc2 <- subset(Table_MB48_FDR_0.01, abs(log2FoldChange) > 2)
Table_MB48_FDR_0.000001 <- subset(MB48_0.05_Full, padj <= 0.000001)
Table_MB48_FDR_0.000001_l2fc2 <- subset(Table_MB48_FDR_0.000001, abs(log2FoldChange) > 2)

Table_MB48_FDR_0.05 <- MB48_0.05_Full
Table_MB48_FDR_0.05_0p <- subset(MB48_0.05_Full, log2FoldChange > 0)
Table_MB48_FDR_0.05_0n <- subset(MB48_0.05_Full, log2FoldChange < 0)
Table_MB48_FDR_0.05_1p <- subset(MB48_0.05_Full, log2FoldChange > 1)
Table_MB48_FDR_0.05_1n <- subset(MB48_0.05_Full, log2FoldChange < -1)
Table_MB48_FDR_0.01_0p <- subset(Table_MB48_FDR_0.01, log2FoldChange > 0)
Table_MB48_FDR_0.01_0n <- subset(Table_MB48_FDR_0.01, log2FoldChange < 0)
Table_MB48_FDR_0.01_1p <- subset(Table_MB48_FDR_0.01, log2FoldChange > 1)
Table_MB48_FDR_0.01_1n <- subset(Table_MB48_FDR_0.01, log2FoldChange < -1)

#MB2
Table_MB2_FDR_0.05_l2fc1 <- subset(MB2_FDR_0.05, abs(log2FoldChange) > 1)
Table_MB2_FDR_0.01 <- subset(MB2_FDR_0.05, padj <= 0.01)
Table_MB2_FDR_0.01_l2fc2 <- subset(Table_MB2_FDR_0.01, abs(log2FoldChange) > 2)
Table_MB2_FDR_0.000001 <- subset(MB2_FDR_0.05, padj <= 0.000001)
Table_MB2_FDR_0.000001_l2fc2 <- subset(Table_MB2_FDR_0.000001, abs(log2FoldChange) > 2)

Table_MB2_FDR_0.05 <- MB2_FDR_0.05
Table_MB2_FDR_0.05_0p <- subset(MB2_FDR_0.05, log2FoldChange > 0)
Table_MB2_FDR_0.05_0n <- subset(MB2_FDR_0.05, log2FoldChange < 0)
Table_MB2_FDR_0.05_1p <- subset(MB2_FDR_0.05, log2FoldChange > 1)
Table_MB2_FDR_0.05_1n <- subset(MB2_FDR_0.05, log2FoldChange < -1)
Table_MB2_FDR_0.01_0p <- subset(Table_MB2_FDR_0.01, log2FoldChange > 0)
Table_MB2_FDR_0.01_0n <- subset(Table_MB2_FDR_0.01, log2FoldChange < 0)
Table_MB2_FDR_0.01_1p <- subset(Table_MB2_FDR_0.01, log2FoldChange > 1)
Table_MB2_FDR_0.01_1n <- subset(Table_MB2_FDR_0.01, log2FoldChange < -1)

#Change the id's of genes (if required) via 
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("org.Bt.eg.db", version = "3.8")
library(org.Bt.eg.db)

rownames <- as.data.frame(mapIds(org.Bt.eg.db, row.names(MB2_res), 'SYMBOL', 'ENTREZID'))
rownames(MB2)<-MB2[,1]
row.names(MB2_res) <- make.names(rownames$`mapIds(org.Bt.eg.db, row.names(MB2_res), "SYMBOL", "ENTREZID")`, unique = TRUE)

rownames <- as.data.frame(mapIds(org.Bt.eg.db, row.names(MB6_res), 'SYMBOL', 'ENTREZID'))
rownames(MB6)<-MB6[,1]
row.names(MB6) <- make.names(rownames$`mapIds(org.Bt.eg.db, row.names(MB6_res), "SYMBOL", "ENTREZID")`, unique = TRUE)

rownames <- as.data.frame(mapIds(org.Bt.eg.db, row.names(MB24_res), 'SYMBOL', 'ENTREZID'))
rownames(MB24)<-MB24[,1]
row.names(MB24) <- make.names(rownames$`mapIds(org.Bt.eg.db, row.names(MB24_res), "SYMBOL", "ENTREZID")`, unique = TRUE)

rownames <- as.data.frame(mapIds(org.Bt.eg.db, row.names(MB48_res), 'SYMBOL', 'ENTREZID'))
rownames(MB48)<-MB48[,1]
row.names(MB48) <- make.names(rownames$`mapIds(org.Bt.eg.db, row.names(MB48_res), "SYMBOL", "ENTREZID")`, unique = TRUE)

rownames <- as.data.frame(mapIds(org.Bt.eg.db, row.names(TB2_res), 'SYMBOL', 'ENTREZID'))
rownames(TB2)<-TB2[,1]
row.names(TB2) <- make.names(rownames$`mapIds(org.Bt.eg.db, row.names(TB2_res), "SYMBOL", "ENTREZID")`, unique = TRUE)

rownames <- as.data.frame(mapIds(org.Bt.eg.db, row.names(TB6_res), 'SYMBOL', 'ENTREZID'))
rownames(TB6)<-TB6[,1]
row.names(TB6) <- make.names(rownames$`mapIds(org.Bt.eg.db, row.names(TB6_res), "SYMBOL", "ENTREZID")`, unique = TRUE)

rownames <- as.data.frame(mapIds(org.Bt.eg.db, row.names(TB24_res), 'SYMBOL', 'ENTREZID'))
rownames(TB24)<-TB24[,1]
row.names(TB24) <- make.names(rownames$`mapIds(org.Bt.eg.db, row.names(TB24_res), "SYMBOL", "ENTREZID")`, unique = TRUE)

rownames <- as.data.frame(mapIds(org.Bt.eg.db, row.names(TB48_res), 'SYMBOL', 'ENTREZID'))
rownames(TB48)<-TB48[,1]
row.names(TB48) <- make.names(rownames$`mapIds(org.Bt.eg.db, row.names(TB48_res), "SYMBOL", "ENTREZID")`, unique = TRUE)

#write results to .csv
write.csv(as.data.frame(MB2), file="MB_2hpi_DE_genes.csv")
write.csv(as.data.frame(MB6), file="MB_6hpi_DE_genes.csv")
write.csv(as.data.frame(MB24), file="MB_24hpi_DE_genes.csv")
write.csv(as.data.frame(MB48), file="MB_48hpi_DE_genes.csv")

write.csv(as.data.frame(TB2), file="TB_2hpi_DE_genes.csv")
write.csv(as.data.frame(TB6), file="TB_6hpi_DE_genes.csv")
write.csv(as.data.frame(TB24), file="TB_24hpi_DE_genes.csv")
write.csv(as.data.frame(TB48), file="TB_48hpi_DE_genes.csv")

#results restricted to FDR < 0.05

MB2_FDR_0.05<-as.data.frame(subset(MB2_res, padj < 0.05))
MB6_FDR_0.05<-as.data.frame(subset(MB6_res, padj < 0.05))
MB24_FDR_0.05<-as.data.frame(subset(MB24_res, padj < 0.05))
MB48_FDR_0.05<-as.data.frame(subset(MB48_res, padj < 0.05))
TB2_FDR_0.05<-as.data.frame(subset(TB2_res, padj < 0.05))
TB6_FDR_0.05<-as.data.frame(subset(TB6_res, padj < 0.05))
TB24_FDR_0.05<-as.data.frame(subset(TB24_res, padj < 0.05))
TB48_FDR_0.05<-as.data.frame(subset(TB48_res, padj < 0.05))

#gene names restricted to FDR < 0.05 

MB2_FDR_0.05_genes<-as.data.frame(rownames(subset(MB2_res, padj < 0.05)))
MB6_FDR_0.05_genes<-as.data.frame(rownames(subset(MB6_res, padj < 0.05)))
MB24_FDR_0.05_genes<-as.data.frame(rownames(subset(MB24_res, padj < 0.05)))
MB48_FDR_0.05_genes<-as.data.frame(rownames(subset(MB48_res, padj < 0.05)))
TB2_FDR_0.05_genes<-as.data.frame(rownames(subset(TB2_res, padj < 0.05)))
TB6_FDR_0.05_genes<-as.data.frame(rownames(subset(TB6_res, padj < 0.05)))
TB24_FDR_0.05_genes<-as.data.frame(rownames(subset(TB24_res, padj < 0.05)))
TB48_FDR_0.05_genes<-as.data.frame(rownames(subset(TB48_res, padj < 0.05)))

#LF2C restricted to FDR < 0.05

MB2_FDR_0.05_genes<-as.data.frame(rownames(subset(MB2_res, padj < 0.05)))
MB6_FDR_0.05_genes<-as.data.frame(rownames(subset(MB6_res, padj < 0.05)))
MB24_FDR_0.05_genes<-as.data.frame(rownames(subset(MB24_res, padj < 0.05)))
MB48_FDR_0.05_genes<-as.data.frame(rownames(subset(MB48_res, padj < 0.05)))
TB2_FDR_0.05_genes<-as.data.frame(rownames(subset(TB2_res, padj < 0.05)))
TB6_FDR_0.05_genes<-as.data.frame(rownames(subset(TB6_res, padj < 0.05)))
TB24_FDR_0.05_genes<-as.data.frame(rownames(subset(TB24_res, padj < 0.05)))
TB48_FDR_0.05_genes <- as.data.frame(rownames(subset(TB48_res, padj < 0.05)))

#Full DE sets with no restriction
MB2_DE<-as.data.frame(subset(MB2_res, padj > 0))
MB6_DE<-as.data.frame(subset(MB6_res, padj > 0))
MB24_DE<-as.data.frame(subset(MB24_res, padj > 0))
MB48_DE<-as.data.frame(subset(MB48_res, padj > 0))

MB2_DE$entrezgene <- row.names(MB2_DE)
MB6_DE$entrezgene <- row.names(MB6_DE)
MB24_DE$entrezgene <- row.names(MB24_DE)
MB48_DE$entrezgene <- row.names(MB48_DE)

########################################################################################
#############################         Annotation        ################################
########################################################################################

#Gene coordinates, as well as other information is needed for donwstream analysis and integration. 
#This pipeline is convulted, mainly because the annotation on most databases is still out of date. It will be refined in the future i.e. only the biomaRt step.


#first check available ID's 
library(org.Bt.eg.db)
columns(org.Bt.eg.db)

#Then append ID's based on row names 
#Convert one list of genes to another id type. First '' is what you want, second '' is what you have i.e the first line converst Entrez to Symbol

MB2_FDR_0.05$SYMBOL <- (mapIds(org.Bt.eg.db, row.names(MB2_FDR_0.05), 'SYMBOL','ENTREZID'))
MB2_FDR_0.05$ENSEMBL <- (mapIds(org.Bt.eg.db, row.names(MB2_FDR_0.05), 'ENSEMBL','ENTREZID'))
MB2_FDR_0.05$ENSEMBLTRANSCRIPT <- (mapIds(org.Bt.eg.db, row.names(MB2_FDR_0.05), 'ENSEMBLTRANS','ENTREZID'))
MB2_FDR_0.05$GENENAME <- (mapIds(org.Bt.eg.db, row.names(MB2_FDR_0.05), 'GENENAME','ENTREZID'))
MB2_FDR_0.05$entrezgene <- (mapIds(org.Bt.eg.db, row.names(MB2_FDR_0.05), 'ENTREZID','ENTREZID')) #this is for combining later on
MB2_FDR_0.05$HS_entrezgene <- (mapIds(org.Hs.eg.db, MB2_FDR_0.05$SYMBOL, 'ENTREZID', 'SYMBOL')) #for human orthologs. This list will be incomplete.

MB24_FDR_0.05$SYMBOL <- (mapIds(org.Bt.eg.db, row.names(MB24_FDR_0.05), 'SYMBOL','ENTREZID'))
MB24_FDR_0.05$ENSEMBL <- (mapIds(org.Bt.eg.db, row.names(MB24_FDR_0.05), 'ENSEMBL','ENTREZID'))
MB24_FDR_0.05$ENSEMBLTRANSCRIPT <- (mapIds(org.Bt.eg.db, row.names(MB24_FDR_0.05), 'ENSEMBLTRANS','ENTREZID'))
MB24_FDR_0.05$GENENAME <- (mapIds(org.Bt.eg.db, row.names(MB24_FDR_0.05), 'GENENAME','ENTREZID'))
MB24_FDR_0.05$entrezgene <- (mapIds(org.Bt.eg.db, row.names(MB24_FDR_0.05), 'ENTREZID','ENTREZID')) 
MB24_FDR_0.05$HS_entrezgene <- (mapIds(org.Hs.eg.db, MB24_FDR_0.05$SYMBOL, 'ENTREZID', 'SYMBOL'))

MB6_FDR_0.05$SYMBOL <- (mapIds(org.Bt.eg.db, row.names(MB6_FDR_0.05), 'SYMBOL','ENTREZID'))
MB6_FDR_0.05$ENSEMBL <- (mapIds(org.Bt.eg.db, row.names(MB6_FDR_0.05), 'ENSEMBL','ENTREZID'))
MB6_FDR_0.05$ENSEMBLTRANSCRIPT <- (mapIds(org.Bt.eg.db, row.names(MB6_FDR_0.05), 'ENSEMBLTRANS','ENTREZID'))
MB6_FDR_0.05$GENENAME <- (mapIds(org.Bt.eg.db, row.names(MB6_FDR_0.05), 'GENENAME','ENTREZID'))
MB6_FDR_0.05$entrezgene <- (mapIds(org.Bt.eg.db, row.names(MB6_FDR_0.05), 'ENTREZID','ENTREZID')) 
MB6_FDR_0.05$HS_entrezgene <- (mapIds(org.Hs.eg.db, MB6_FDR_0.05$SYMBOL, 'ENTREZID', 'SYMBOL'))

MB48_FDR_0.05$SYMBOL <- (mapIds(org.Bt.eg.db, row.names(MB48_FDR_0.05), 'SYMBOL','ENTREZID'))
MB48_FDR_0.05$ENSEMBL <- (mapIds(org.Bt.eg.db, row.names(MB48_FDR_0.05), 'ENSEMBL','ENTREZID'))
MB48_FDR_0.05$ENSEMBLTRANSCRIPT <- (mapIds(org.Bt.eg.db, row.names(MB48_FDR_0.05), 'ENSEMBLTRANS','ENTREZID'))
MB48_FDR_0.05$GENENAME <- (mapIds(org.Bt.eg.db, row.names(MB48_FDR_0.05), 'GENENAME','ENTREZID'))
MB48_FDR_0.05$entrezgene <- (mapIds(org.Bt.eg.db, row.names(MB48_FDR_0.05), 'ENTREZID','ENTREZID')) 
MB48_FDR_0.05$HS_entrezgene <- (mapIds(org.Hs.eg.db, MB48_FDR_0.05$SYMBOL, 'ENTREZID', 'SYMBOL'))

TB2_FDR_0.05$SYMBOL <- (mapIds(org.Bt.eg.db, row.names(TB2_FDR_0.05), 'SYMBOL','ENTREZID'))
TB2_FDR_0.05$ENSEMBL <- (mapIds(org.Bt.eg.db, row.names(TB2_FDR_0.05), 'ENSEMBL','ENTREZID'))
TB2_FDR_0.05$ENSEMBLTRANSCRIPT <- (mapIds(org.Bt.eg.db, row.names(TB2_FDR_0.05), 'ENSEMBLTRANS','ENTREZID'))
TB2_FDR_0.05$GENENAME <- (mapIds(org.Bt.eg.db, row.names(TB2_FDR_0.05), 'GENENAME','ENTREZID'))
TB2_FDR_0.05$entrezgene <- (mapIds(org.Bt.eg.db, row.names(TB2_FDR_0.05), 'ENTREZID','ENTREZID'))
TB2_FDR_0.05$HS_entrezgene <- (mapIds(org.Hs.eg.db, TB2_FDR_0.05$SYMBOL, 'ENTREZID', 'SYMBOL'))

TB24_FDR_0.05$SYMBOL <- (mapIds(org.Bt.eg.db, row.names(TB24_FDR_0.05), 'SYMBOL','ENTREZID'))
TB24_FDR_0.05$ENSEMBL <- (mapIds(org.Bt.eg.db, row.names(TB24_FDR_0.05), 'ENSEMBL','ENTREZID'))
TB24_FDR_0.05$ENSEMBLTRANSCRIPT <- (mapIds(org.Bt.eg.db, row.names(TB24_FDR_0.05), 'ENSEMBLTRANS','ENTREZID'))
TB24_FDR_0.05$GENENAME <- (mapIds(org.Bt.eg.db, row.names(TB24_FDR_0.05), 'GENENAME','ENTREZID'))
TB24_FDR_0.05$entrezgene <- (mapIds(org.Bt.eg.db, row.names(TB24_FDR_0.05), 'ENTREZID','ENTREZID')) 
TB24_FDR_0.05$HS_entrezgene <- (mapIds(org.Hs.eg.db, TB24_FDR_0.05$SYMBOL, 'ENTREZID', 'SYMBOL'))

TB6_FDR_0.05$SYMBOL <- (mapIds(org.Bt.eg.db, row.names(TB6_FDR_0.05), 'SYMBOL','ENTREZID'))
TB6_FDR_0.05$ENSEMBL <- (mapIds(org.Bt.eg.db, row.names(TB6_FDR_0.05), 'ENSEMBL','ENTREZID'))
TB6_FDR_0.05$ENSEMBLTRANSCRIPT <- (mapIds(org.Bt.eg.db, row.names(TB6_FDR_0.05), 'ENSEMBLTRANS','ENTREZID'))
TB6_FDR_0.05$GENENAME <- (mapIds(org.Bt.eg.db, row.names(TB6_FDR_0.05), 'GENENAME','ENTREZID'))
TB6_FDR_0.05$entrezgene <- (mapIds(org.Bt.eg.db, row.names(TB6_FDR_0.05), 'ENTREZID','ENTREZID')) 
TB6_FDR_0.05$HS_entrezgene <- (mapIds(org.Hs.eg.db, TB6_FDR_0.05$SYMBOL, 'ENTREZID', 'SYMBOL'))

TB48_FDR_0.05$SYMBOL <- (mapIds(org.Bt.eg.db, row.names(TB48_FDR_0.05), 'SYMBOL','ENTREZID'))
TB48_FDR_0.05$ENSEMBL <- (mapIds(org.Bt.eg.db, row.names(TB48_FDR_0.05), 'ENSEMBL','ENTREZID'))
TB48_FDR_0.05$ENSEMBLTRANSCRIPT <- (mapIds(org.Bt.eg.db, row.names(TB48_FDR_0.05), 'ENSEMBLTRANS','ENTREZID'))
TB48_FDR_0.05$GENENAME <- (mapIds(org.Bt.eg.db, row.names(TB48_FDR_0.05), 'GENENAME','ENTREZID'))
TB48_FDR_0.05$entrezgene <- (mapIds(org.Bt.eg.db, row.names(TB48_FDR_0.05), 'ENTREZID','ENTREZID')) 
TB48_FDR_0.05$HS_entrezgene <- (mapIds(org.Hs.eg.db, TB48_FDR_0.05$SYMBOL, 'ENTREZID', 'SYMBOL'))

#TXdb for gene coordinates (This currently has UMD sites, so is out of date. This is why we use stupid biomaRt, which has the current gene coordinates)
cols <- c( "TXCHROM", "TXSTART", "TXEND", "TXSTRAND")

MB24_FDR_0.05_coordinates <- select(TxDb.Btaurus.UCSC.bosTau8.refGene, keys=row.names(MB24_FDR_0.05), columns=cols, keytype="GENEID")
MB48_FDR_0.05_coordinates <- select(TxDb.Btaurus.UCSC.bosTau8.refGene, keys=row.names(MB48_FDR_0.05), columns=cols, keytype="GENEID")
MB24_FDR_0.05_coordinates <- select(TxDb.Btaurus.UCSC.bosTau8.refGene, keys=row.names(MB24_FDR_0.05), columns=cols, keytype="GENEID")
MB48_FDR_0.05_coordinates <- select(TxDb.Btaurus.UCSC.bosTau8.refGene, keys=row.names(MB48_FDR_0.05), columns=cols, keytype="GENEID")

MB24_FDR_0.05_coordinates <- select(TxDb.Btaurus.UCSC.bosTau8.refGene, keys=row.names(MB24_FDR_0.05), columns=cols, keytype="GENEID")
MB48_FDR_0.05_coordinates <- select(TxDb.Btaurus.UCSC.bosTau8.refGene, keys=row.names(MB48_FDR_0.05), columns=cols, keytype="GENEID")
MB24_FDR_0.05_coordinates <- select(TxDb.Btaurus.UCSC.bosTau8.refGene, keys=row.names(MB24_FDR_0.05), columns=cols, keytype="GENEID")
MB48_FDR_0.05_coordinates <- select(TxDb.Btaurus.UCSC.bosTau8.refGene, keys=row.names(MB48_FDR_0.05), columns=cols, keytype="GENEID")

#Get annotation data for the FDR < 0.05

install.packages("biomaRt")
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "btaurus_gene_ensembl")

MB2_annotation_0.05 <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'external_gene_name', 'entrezgene', 'refseq_mrna', 'description', 'chromosome_name',
                                             'start_position', 'end_position', 'strand' ),  
                              filter = 'entrezgene',
                              values = MB2_FDR_0.05_genes$`rownames(subset(MB2_res, padj < 0.05))`, 
                              mart = ensembl)

MB24_annotation_0.05 <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'external_gene_name', 'entrezgene', 'refseq_mrna', 'description', 'chromosome_name',
                                             'start_position', 'end_position', 'strand' ),  
                              filter = 'entrezgene',
                              values = MB24_FDR_0.05_genes$`rownames(subset(MB24_res, padj < 0.05))`, 
                              mart = ensembl)							  
							 
MB6_annotation_0.05 <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'external_gene_name', 'entrezgene', 'refseq_mrna', 'description', 'chromosome_name',
                                             'start_position', 'end_position', 'strand' ),  
                              filter = 'entrezgene',
                              values = MB6_FDR_0.05_genes$`rownames(subset(MB6_res, padj < 0.05))`, 
                              mart = ensembl)
							  
MB48_annotation_0.05 <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'external_gene_name', 'entrezgene', 'refseq_mrna', 'description', 'chromosome_name',
                                             'start_position', 'end_position', 'strand' ),  
                              filter = 'entrezgene',
                              values = MB48_FDR_0.05_genes$`rownames(subset(MB48_res, padj < 0.05))`, 
                              mart = ensembl)

TB2_annotation_0.05 <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'external_gene_name', 'entrezgene', 'refseq_mrna', 'description', 'chromosome_name',
                                             'start_position', 'end_position', 'strand' ),  
                              filter = 'entrezgene',
                              values = TB2_FDR_0.05_genes$`rownames(subset(TB2_res, padj < 0.05))`, 
                              mart = ensembl)

TB24_annotation_0.05 <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'external_gene_name', 'entrezgene', 'refseq_mrna', 'description', 'chromosome_name',
                                             'start_position', 'end_position', 'strand' ),  
                              filter = 'entrezgene',
                              values = TB24_FDR_0.05_genes$`rownames(subset(TB24_res, padj < 0.05))`, 
                              mart = ensembl)							  
							 
TB6_annotation_0.05 <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'external_gene_name', 'entrezgene', 'refseq_mrna', 'description', 'chromosome_name',
                                             'start_position', 'end_position', 'strand' ),  
                              filter = 'entrezgene',
                              values = TB6_FDR_0.05_genes$`rownames(subset(TB6_res, padj < 0.05))`, 
                              mart = ensembl)
							  
TB48_annotation_0.05 <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'external_gene_name', 'entrezgene', 'refseq_mrna', 'description', 'chromosome_name',
                                             'start_position', 'end_position', 'strand' ),  
                              filter = 'entrezgene',
                              values = TB48_FDR_0.05_genes$`rownames(subset(TB48_res, padj < 0.05))`, 
                              mart = ensembl)


#Combine Annotation files. First remove duplicated rows based on entrez id's. The coordinates for each duplicted entry are the same 
MB2_annotation_0.05 <- MB2_annotation_0.05[!duplicated(MB2_annotation_0.05$entrezgene), ]	
MB24_annotation_0.05 <- MB24_annotation_0.05[!duplicated(MB24_annotation_0.05$entrezgene), ]	
MB6_annotation_0.05 <- MB6_annotation_0.05[!duplicated(MB6_annotation_0.05$entrezgene), ]	
MB48_annotation_0.05 <- MB48_annotation_0.05[!duplicated(MB48_annotation_0.05$entrezgene), ]	

TB2_annotation_0.05 <- TB2_annotation_0.05[!duplicated(TB2_annotation_0.05$entrezgene), ]	
TB24_annotation_0.05 <- TB24_annotation_0.05[!duplicated(TB24_annotation_0.05$entrezgene), ]	
TB6_annotation_0.05 <- TB6_annotation_0.05[!duplicated(TB6_annotation_0.05$entrezgene), ]	
TB48_annotation_0.05 <- TB48_annotation_0.05[!duplicated(TB48_annotation_0.05$entrezgene), ]	

#Then combine the two files based on Entrez id's 						  
MB2_0.05_Full <- merge(MB2_FDR_0.05, MB2_annotation_0.05, by='entrezgene', all=TRUE)				
MB24_0.05_Full <- merge(MB24_FDR_0.05, MB24_annotation_0.05, by='entrezgene', all=TRUE)					
MB6_0.05_Full <- merge(MB6_FDR_0.05, MB6_annotation_0.05, by='entrezgene', all=TRUE)				
MB48_0.05_Full <- merge(MB48_FDR_0.05, MB48_annotation_0.05, by='entrezgene', all=TRUE)

TB2_0.05_Full <- merge(TB2_FDR_0.05, TB2_annotation_0.05, by='entrezgene', all=TRUE)				
TB24_0.05_Full <- merge(TB24_FDR_0.05, TB24_annotation_0.05, by='entrezgene', all=TRUE)					
TB6_0.05_Full <- merge(TB6_FDR_0.05, TB6_annotation_0.05, by='entrezgene', all=TRUE)				
TB48_0.05_Full <- merge(TB48_FDR_0.05, TB48_annotation_0.05, by='entrezgene', all=TRUE)

#If desired, you can subset the final list by logFC as well, and restrict FDR cutoff. To reduce the DE gene list to below 500, 
#I reduced MB24 & 48 to logfc <= -2 / >= 2 and MB48 to padj <= 0.00001

MB24_0_05_LFC_Full_DE <- subset(MB24_0_05_Full_DE, log2FoldChange >= 2 | log2FoldChange <= -2)
MB48_0_05_LFC_Full_DE <- subset(MB48_0_05_Full_DE, log2FoldChange >= 2 | log2FoldChange <= -2)
 
MB48_0_05_LFC_Full_DE <- subset(MB48_0_05_LFC_Full_DE, padj <= 0.000001)
MB24_0_05_LFC_Full_DE <- subset(MB24_0_05_LFC_Full_DE, padj <= 0.01)

write.csv(MB24_0_05_LFC_Full_DE, file = "MB24_0_05_LFC_Full_DE.csv")
write.csv(MB48_0_05_LFC_Full_DE, file = "MB48_0_05_LFC_Full_DE.csv")

#write the results for GWAS integration 
write.csv(MB2_0.05_Full , file="MB2_0.05_Full.csv")
write.csv(MB24_0.05_Full , file="MB24_0.05_Full.csv")
write.csv(MB48_0.05_Full , file="MB48_0.05_Full.csv")
write.csv(MB6_0.05_Full , file="MB6_0.05_Full.csv")

write.csv(TB2_0.05_Full , file="TB2_0.05_Full.csv")
write.csv(TB24_0.05_Full , file="TB24_0.05_Full.csv")
write.csv(TB48_0.05_Full , file="TB48_0.05_Full.csv")
write.csv(TB6_0.05_Full , file="TB6_0.05_Full.csv")			


##Full DE list for MB###

MB2_DE_annotation<- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'external_gene_name', 'entrezgene', 'refseq_mrna', 'description', 'chromosome_name',
                                           'start_position', 'end_position', 'strand' ),,  
                             filter = 'entrezgene',
                             values = MB2_DE$entrezgene, 
                             mart = ensembl)
							 
MB24_DE_annotation <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'external_gene_name', 'entrezgene', 'refseq_mrna', 'description', 'chromosome_name',
                                             'start_position', 'end_position', 'strand' ),  
                              filter = 'entrezgene',
                              values = MB6_DE$entrezgene, 
                              mart = ensembl)							  
							 
MB6_DE_annotation <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'external_gene_name', 'entrezgene', 'refseq_mrna', 'description', 'chromosome_name',
                                             'start_position', 'end_position', 'strand' ),  
                              filter = 'entrezgene',
                              values = MB6_DE$entrezgene, 
                              mart = ensembl)
							  
MB48_DE_annotation <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'external_gene_name', 'entrezgene', 'refseq_mrna', 'description', 'chromosome_name',
                                             'start_position', 'end_position', 'strand' ),  
                              filter = 'entrezgene',
                              values = MB48_DE$entrezgene, 
                              mart = ensembl)

MB2_DE_annotation <- MB2_DE_annotation[!duplicated(MB2_DE_annotation$entrezgene), ]	
MB24_DE_annotation <- MB24_DE_annotation[!duplicated(MB24_DE_annotation$entrezgene), ]	
MB6_DE_annotation <- MB6_DE_annotation[!duplicated(MB6_DE_annotation$entrezgene), ]	
MB48_DE_annotation <- MB48_DE_annotation[!duplicated(MB48_DE_annotation$entrezgene), ]	

MB2_DE_full <- merge(MB2_DE, MB2_DE_annotation, by='entrezgene', all=TRUE)				
MB24_DE_full <- merge(MB24_DE, MB24_DE_annotation, by='entrezgene', all=TRUE)					
MB6_DE_full <- merge(MB6_DE, MB6_DE_annotation, by='entrezgene', all=TRUE)				
MB48_DE_full<- merge(MB48_DE, MB48_DE_annotation, by='entrezgene', all=TRUE)

write.csv(MB2_DE_full , file="MB2_DE_fulll.csv")
write.csv(MB24_DE_full , file="MB24_DE_full.csv")
write.csv(MB6_DE_full , file="MB6_DE_full.csv")
write.csv(MB48_DE_full , file="MB48_DE_full.csv")

#####Conversion to human genes#####
# Due to the outdated annotation on ENSEMBL, it takes a few ways to convert bovine genes to human genes.  
# Basic function to convert human to mouse gene names


BovineGenes <- c(x)

convertBovineGeneList <- function(x){
require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
cow = useMart("ensembl", dataset = "btaurus_gene_ensembl")
genesV2 = getLDS(attributes = c("entrezgene"), filters = "entrezgene", values = x , mart = cow, attributesL = c("entrezgene"), martL = human, uniqueRows=T)
humanx <- unique(genesV2[, 2])
 
# Print the first 6 genes found to the screen
}

#####Annotation for the DE network files#####

MB24_0_05_Full_DE$HS_entrezgene <- (mapIds(org.Hs.eg.db, MB24_0_05_Full_DE$SYMBOL, 'ENTREZID', 'SYMBOL'))
MB48_0_05_Full_DE$HS_entrezgene <- (mapIds(org.Hs.eg.db, MB48_0_05_Full_DE$SYMBOL, 'ENTREZID', 'SYMBOL'))

write.csv(MB24_0_05_Full_DE , file="MB24_0.05_Full.csv")
write.csv(MB48_0_05_Full_DE , file="MB48_0.05_Full.csv")

#And coversion and population of DEN modules to Bovine ID's. Needed for gene coordiantes
MB24_Top5_Full_DEN$entrezgene <- (mapIds(org.Bt.eg.db, MB24_Top5_Full_DEN$SYMBOL, 'ENTREZID', 'SYMBOL'))
MB48_Top5_Full_DEN$entrezgene <- (mapIds(org.Bt.eg.db, MB48_Top5_Full_DEN$SYMBOL, 'ENTREZID', 'SYMBOL'))

MB24_Top5_Full_DEN$ENSEMBL <- (mapIds(org.Bt.eg.db, MB24_Top5_Full_DEN$SYMBOL, 'ENSEMBL', 'SYMBOL'))
MB48_Top5_Full_DEN$ENSEMBL <- (mapIds(org.Bt.eg.db, MB48_Top5_Full_DEN$SYMBOL, 'ENSEMBL', 'SYMBOL'))

MB24_Top5_Full_DEN$ENSEMBLTRANSCRIPT <- (mapIds(org.Bt.eg.db, MB24_Top5_Full_DEN$SYMBOL, 'ENSEMBLTRANS', 'SYMBOL'))
MB48_Top5_Full_DEN$ENSEMBLTRANSCRIPT <- (mapIds(org.Bt.eg.db, MB48_Top5_Full_DEN$SYMBOL, 'ENSEMBLTRANS', 'SYMBOL'))

MB24_Top5_Full_DEN$GENENAME <- (mapIds(org.Bt.eg.db, MB24_Top5_Full_DEN$SYMBOL, 'GENENAME', 'SYMBOL'))
MB48_Top5_Full_DEN$GENENAME <- (mapIds(org.Bt.eg.db, MB48_Top5_Full_DEN$SYMBOL, 'GENENAME', 'SYMBOL'))

MB24_Top5_annotation <- (getBM(attributes = c('entrezgene', 'description', 'chromosome_name', 'start_position', 'end_position', 'strand' ),  
                                   filter = 'entrezgene',
                                   values = MB24_Top5_Full_DEN$entrezgene, 
                                   mart = ensembl))

MB24_Top5_Full_DEN <- merge(MB24_Top5_Full_DEN, MB24_Top5_annotation, by = "entrezgene", sort = TRUE)

MB48_Top5_annotation <- (getBM(attributes = c('entrezgene', 'description', 'chromosome_name', 'start_position', 'end_position', 'strand' ),  
                                   filter = 'entrezgene',
                                   values = MB48_Top5_Full_DEN$entrezgene, 
                                   mart = ensembl))

MB48_Top5_Full_DEN <- merge(MB48_Top5_Full_DEN, MB48_Top5_annotation, by = "entrezgene", sort = TRUE)


write.csv(MB24_Top5_Full_DEN , file="MB24_Top5_Full_DEN.csv")
write.csv(MB48_Top5_Full_DEN , file="MB48_Top5_Full_DEN.csv")
