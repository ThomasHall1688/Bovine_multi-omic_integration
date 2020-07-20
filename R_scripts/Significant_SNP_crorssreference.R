######################################################################################
############################ Crossreferencing hits ###################################
######################################################################################
#If there are 100's of SNP hits, its too much to do to find the SNP metadata by hand, so we combine more datasets. 
#First we need the GWAS files 

ARS_SNPs = read.delim(file = "/home/workspace/thall/Analysis/Networks/GWAS_integration/R_analysis/New_GWAS/GWAS_remapping_data/raw_GWAS/ARS1.2PlusY_BQSR.vcf", +
						  header = F, sep = "\t", fill = TRUE)
ARS_SNPs <- tail(ARS_SNPs, -4)
colnames(ARS_SNPs) <- as.character(unlist(ARS_SNPs[1,]))
ARS_SNPs = ARS_SNPs[-1, ]
colnames(ARS_SNPs)[1] <- "CHROM"
ARS_SNPs <- ARS_SNPs[,-c(6,7,8)]
colnames(ARS_SNPs) <- c("CHROM", "POS", "ID", "REF", "ALT")
						 
significant_SNPs_HOFR <- read.csv(file = "Significant_SNPs_HOFR.csv")
significant_SNPs_LM <-  read.csv(file = "Significant_SNPs_LM.csv")
significant_SNPs_CH <-  read.csv(file = "Significant_SNPs_CH.csv")

#normalize column names
colnames(significant_SNPs_HOFR) <-  c("GWAS_set", "File_name", "Analysis_type", "ID", "P_value", "Original_qvalue", "New_qvalue", "Probability_%")
colnames(significant_SNPs_CH) <-  c("GWAS_set", "File_name", "Analysis_type", "ID", "P_value", "Original_qvalue", "New_qvalue", "Probability_%")
colnames(significant_SNPs_LM) <-  c("GWAS_set", "File_name", "Analysis_type", "ID", "P_value", "Original_qvalue", "New_qvalue", "Probability_%")

#Remove duplications, leaving the lowest qvalues
Test <- significant_SNPs_HOFR[order(significant_SNPs_HOFR$ID, -abs(significant_SNPs_HOFR$New_qvalue) ), ] #sort by id and reverse of abs(value)
Test <- Test[ !duplicated(Test$ID), ]

significant_SNPs_HOFR <- significant_SNPs_HOFR[order(significant_SNPs_HOFR$ID, -abs(significant_SNPs_HOFR$New_qvalue) ), ] #sort by id and reverse of abs(value)
significant_SNPs_HOFR <- significant_SNPs_HOFR[ !duplicated(significant_SNPs_HOFR$ID), ]

significant_SNPs_CH <- significant_SNPs_CH[order(significant_SNPs_CH$ID, -abs(significant_SNPs_CH$New_qvalue) ), ] #sort by id and reverse of abs(value)
significant_SNPs_CH <- significant_SNPs_CH[ !duplicated(significant_SNPs_CH$ID), ]

significant_SNPs_LM <- significant_SNPs_LM[order(significant_SNPs_LM$ID, -abs(significant_SNPs_LM$New_qvalue) ), ] #sort by id and reverse of abs(value)
significant_SNPs_LM <- significant_SNPs_LM[ !duplicated(significant_SNPs_LM$ID), ]

#Merge based on ID to get all positions and chromosomes
significant_SNPs_HOFR  <- merge(significant_SNPs_HOFR , ARS_SNPs, by='ID')
significant_SNPs_CH <- merge(significant_SNPs_CH, ARS_SNPs, by='ID')
significant_SNPs_LM <- merge(significant_SNPs_LM, ARS_SNPs, by='ID')

write.csv(significant_SNPs_HOFR, file = "significant_SNPs_HOFR_annotated.csv", row.names = FALSE)
write.csv(significant_SNPs_CH, file = "significant_SNPs_CH_annotated.csv", row.names = FALSE)
write.csv(significant_SNPs_LM, file = "significant_SNPs_LM_annotated.csv", row.names = FALSE)


#In order to merge expression values with the annotated SNP file (combined)

significant_SNPs_HOFR  <- merge(significant_SNPs_HOFR , ARS_SNPs, by='ID')
Significant_SNPs_All$entrezgene <- (mapIds(org.Bt.eg.db, Significant_SNPs_All$SYMBOL, 'ENTREZID', 'SYMBOL'))
Significant_SNPs_All <- merge(Significant_SNPs_All, MB_48hpi_DE_genes_entrez, by='entrezgene', all.x = TRUE)
write.csv(Significant_SNPs_All, file = "Significant_SNPs.csv", row.names = FALSE)
