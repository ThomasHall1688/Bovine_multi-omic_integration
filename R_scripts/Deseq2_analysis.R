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

MB2_res <- results(dds, contrast=c("TT","MB2","CN2")) 
MB6_res <- results(dds, contrast=c("TT","MB6","CN6"))
MB24_res <- results(dds, contrast=c("TT","MB24","CN24"))
MB48_res <- results(dds, contrast=c("TT","MB48","CN48"))
TB2_res <- results(dds, contrast=c("TT","TB2","CN2"))
TB6_res <- results(dds, contrast=c("TT","TB6","CN6"))
TB24_res <- results(dds, contrast=c("TT","TB24","CN24"))
TB48_res <- results(dds, contrast=c("TT","TB48","CN48"))

MB2_LFC2_res <- results(dds, contrast=c("TT","MB2","CN2"), lfcThreshold=2) 
MB6_LFC2_res <- results(dds, contrast=c("TT","MB6","CN6"), lfcThreshold=2)
MB24_LFC2_res <- results(dds, contrast=c("TT","MB24","CN24"), lfcThreshold=2)
MB48_LFC2_res <- results(dds, contrast=c("TT","MB48","CN48"), lfcThreshold=2)
TB2_LFC2_res <- results(dds, contrast=c("TT","TB2","CN2"), lfcThreshold=2)
TB6_LFC2_res <- results(dds, contrast=c("TT","TB6","CN6"), lfcThreshold=2)
TB24_LFC2_res <- results(dds, contrast=c("TT","TB24","CN24"), lfcThreshold=2)
TB48_LFC2_res <- results(dds, contrast=c("TT","TB48","CN48"), lfcThreshold=2)

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

write.csv(as.data.frame(MB2_res), file="MB_2hpi_DE_LFC2_genes.csv")
write.csv(as.data.frame(MB6_res), file="MB_6hpi_DE_LFC2_genes.csv")
write.csv(as.data.frame(MB24_res), file="MB_24hpi_DE_LFC2_genes.csv")
write.csv(as.data.frame(MB48_res), file="MB_48hpi_DE_LFC2_genes.csv")

write.csv(as.data.frame(TB2_res), file="TB_2hpi_DE_LFC2_genes.csv")
write.csv(as.data.frame(TB6_res), file="TB_6hpi_DE_LFC2_genes.csv")
write.csv(as.data.frame(TB24_res), file="TB_24hpi_DE_LFC2_genes.csv")
write.csv(as.data.frame(TB48_res), file="TB_48hpi_DE_LFC2_genes.csv")
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

