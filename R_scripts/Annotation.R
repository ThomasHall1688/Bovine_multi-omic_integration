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
