############### GO analysis on modules #################
#Get Gene ontology terms for each of the sub modules 
#I have added an extra step to change the gene id's for the modules data frame. This is because GOstats require ENSEMBL ID's.
#HGNC will convert, but the author reccommends to do this outside the function of moduleGO, so we will use org.Bt.eg.db
#Universe is the background set. This gene list most likely should be the gene set post-filtering, but prior to differential correlation analysis.

#description 	Takes input vectors of gene symbols, labels of corresponding modules, and a universe gene set and leverages the GOstats package to perform GO enrichment analysis.
#param 			genes A character vector specifying gene symbols, present as rows in the inputMat, corresponding to each module label in the labels argument.
#param 			labels A character vector specifying module label names, one for each gene symbol in the genes argument, with overlap allowed (i.e., each gene can be in more than one module).
#param 			universe Character vector of gene symbols which should be used as the background in the hypergeomtric test. If using this in the context of a DGCA experiment, this gene list most likely should be the gene set post-filtering, but prior to differential correlation analysis.
#param 			HGNC_clean Logical indicating whether the input gene symbols should be switched to clean HGNC symbols using the checkGeneSymbols function from the R package HGNChelper. Only applies if HGNC symbols are inputted.
#param 			pval_GO_cutoff Cutoff for the unadjusted p-values of gene ontology terms in the enrichment tests that should be displayed in the resulting table.
#param 			HGNC_switch Logical indicating whether or not the input gene symbols need to be switched from HGNC to Ensembl, the latter of which is required for GOstats enrichment test. Note that this is done by selecting the first Ensembl symbol that maps to a particular HGNC symbol, 
#	   			which is not always unique. If you need more precision on the conversion, you should do this outside of the function and insert the Ensembl list to the function.
#param 			gene_ontology A string specifying the branch of GO that should be used for enrichment analysis. One of "BP" (Biological Process), "MF" (Molecular Function), 
#	   			"CC" (Cellular Component), or "all". If "all" is chosen, then this function finds the enrichment for all of the terms and combines them into one table. Default = "all"
#param 			annotation The library indicating the GO annotation database from which the Go terms should be mapped to gene symbols. Default = "org.Hs.eg.db", which is the table for Homo sapiens. Other common choices include "org.Mm.eg.db", "org.Rn.eg.db". The corresponding annotation library needs to be installed.
#param 			calculateVariance Optionally, find the variance of the odds ratio for each enrichment test. In particular, this finds the standard error of the log odds ratio, which converges to a normal distribution much more quickly than the non-log OR.
#param 			conditional Logical specifying whether the GO analysis should be done conditionally to take into account the hierarchical structure of the GO database in making sense of the gene set enrichments.
#return 		A list of lists of df's, one corresponding to each module, containing GO enrichment information for each module in each of the GO categories selected.


library(GOstats, quietly = TRUE)
library(HGNChelper, quietly = TRUE)
library(org.Bt.eg.db, quietly = TRUE)
columns(org.Bt.eg.db)

#First create a copy of modules_df for editing. THis is because the GO enrichment rarely works when you include moodules that are smaller than 5-6 genes
#Seeing as we are not interested in <20, its best to remove this. I have a line that restricts MEGENA to modules >20, but this is not working, so for now it must be done manually. 

Alv_mac_MF_MB_CN_2hpi_MEGENA_modules_GO <- Alv_mac_MF_MB_CN_2hpi_MEGENA_modules_df
Alv_mac_MF_MB_CN_24hpi_MEGENA_modules_GO <- Alv_mac_MF_MB_CN_24hpi_MEGENA_modules_df
Alv_mac_MF_MB_CN_6hpi_MEGENA_modules_GO <- Alv_mac_MF_MB_CN_6hpi_MEGENA_modules_df
Alv_mac_MF_MB_CN_48hpi_MEGENA_modules_GO <- Alv_mac_MF_MB_CN_48hpi_MEGENA_modules_df
Alv_mac_MF_TB_CN_2hpi_MEGENA_modules_GO <- Alv_mac_MF_TB_CN_2hpi_MEGENA_modules_df
Alv_mac_MF_TB_CN_24hpi_MEGENA_modules_GO <- Alv_mac_MF_TB_CN_24hpi_MEGENA_modules_df
Alv_mac_MF_TB_CN_6hpi_MEGENA_modules_GO <- Alv_mac_MF_TB_CN_6hpi_MEGENA_modules_df
Alv_mac_MF_TB_CN_48hpi_MEGENA_modules_GO <- Alv_mac_MF_TB_CN_48hpi_MEGENA_modules_df

#Next, we want to extract the modules and their genes that are greater than 20. This list is in Alv_mac_MF_MB_CN_24hpi_MEGENA_module_table$id

modules_2hpi <- as.data.frame(Alv_mac_MF_MB_CN_2hpi_MEGENA_module_table$id)
colnames(modules_2hpi) <- c("modules")
Alv_mac_MF_MB_CN_2hpi_MEGENA_modules_GO <- Alv_mac_MF_MB_CN_2hpi_MEGENA_modules_GO[Alv_mac_MF_MB_CN_2hpi_MEGENA_modules_GO$Modules %in% modules_2hpi$modules,]

modules_24hpi <- as.data.frame(Alv_mac_MF_MB_CN_24hpi_MEGENA_module_table$id)
colnames(modules_24hpi) <- c("modules")
Alv_mac_MF_MB_CN_24hpi_MEGENA_modules_GO <- Alv_mac_MF_MB_CN_24hpi_MEGENA_modules_GO[Alv_mac_MF_MB_CN_24hpi_MEGENA_modules_GO$Modules %in% modules_24hpi$modules,]

modules_6hpi <- as.data.frame(Alv_mac_MF_MB_CN_6hpi_MEGENA_module_table$id)
colnames(modules_6hpi) <- c("modules")
Alv_mac_MF_MB_CN_6hpi_MEGENA_modules_GO <- Alv_mac_MF_MB_CN_6hpi_MEGENA_modules_GO[Alv_mac_MF_MB_CN_6hpi_MEGENA_modules_GO$Modules %in% modules_6hpi$modules,]

modules_48hpi <- as.data.frame(Alv_mac_MF_MB_CN_48hpi_MEGENA_module_table$id)
colnames(modules_48hpi) <- c("modules")
Alv_mac_MF_MB_CN_48hpi_MEGENA_modules_GO <- Alv_mac_MF_MB_CN_48hpi_MEGENA_modules_GO[Alv_mac_MF_MB_CN_48hpi_MEGENA_modules_GO$Modules %in% modules_48hpi$modules,]

modules_2hpi_TB <- as.data.frame(Alv_mac_MF_TB_CN_2hpi_MEGENA_module_table$id)
colnames(modules_2hpi_TB) <- c("modules")
Alv_mac_MF_TB_CN_2hpi_MEGENA_modules_GO <- Alv_mac_MF_TB_CN_2hpi_MEGENA_modules_GO[Alv_mac_MF_TB_CN_2hpi_MEGENA_modules_GO$Modules %in% modules_2hpi_TB$modules,]

modules_24hpi_TB <- as.data.frame(Alv_mac_MF_TB_CN_24hpi_MEGENA_module_table$id)
colnames(modules_24hpi_TB) <- c("modules")
Alv_mac_MF_TB_CN_24hpi_MEGENA_modules_GO <- Alv_mac_MF_TB_CN_24hpi_MEGENA_modules_GO[Alv_mac_MF_TB_CN_24hpi_MEGENA_modules_GO$Modules %in% modules_24hpi_TB$modules,]

modules_6hpi_TB <- as.data.frame(Alv_mac_MF_TB_CN_6hpi_MEGENA_module_table$id)
colnames(modules_6hpi_TB) <- c("modules")
Alv_mac_MF_TB_CN_6hpi_MEGENA_modules_GO <- Alv_mac_MF_TB_CN_6hpi_MEGENA_modules_GO[Alv_mac_MF_TB_CN_6hpi_MEGENA_modules_GO$Modules %in% modules_6hpi$modules,]

modules_48hpi_TB <- as.data.frame(Alv_mac_MF_TB_CN_48hpi_MEGENA_module_table$id)
colnames(modules_48hpi_TB) <- c("modules")
Alv_mac_MF_TB_CN_48hpi_MEGENA_modules_GO <- Alv_mac_MF_TB_CN_48hpi_MEGENA_modules_GO[Alv_mac_MF_TB_CN_48hpi_MEGENA_modules_GO$Modules %in% modules_48hpi_TB$modules,]

#Then use the refined module list for GO enrichment 

moduleGO_res_MB2hpi = moduleGO(genes = Alv_mac_MF_MB_CN_2hpi_MEGENA_modules_GO$Genes, labels = Alv_mac_MF_MB_CN_2hpi_MEGENA_modules_GO$Modules,
universe = rownames(Alv_mac_normalized_counts_MF_MB_CN_2hpi_filtered), pval_GO_cutoff = 1, annotation = "org.Bt.eg.db")

moduleGO_res_MB24hpi = moduleGO(genes = Alv_mac_MF_MB_CN_24hpi_MEGENA_modules_GO$Genes, labels = Alv_mac_MF_MB_CN_24hpi_MEGENA_modules_GO$Modules,
universe = rownames(Alv_mac_normalized_counts_MF_MB_CN_24hpi_filtered), pval_GO_cutoff = 1, annotation = "org.Bt.eg.db")

moduleGO_res_MB6hpi = moduleGO(genes = Alv_mac_MF_MB_CN_6hpi_MEGENA_modules_GO$Genes, labels = Alv_mac_MF_MB_CN_6hpi_MEGENA_modules_GO$Modules,
universe = rownames(Alv_mac_normalized_counts_MF_MB_CN_6hpi_filtered), pval_GO_cutoff = 1, annotation = "org.Bt.eg.db")

moduleGO_res_MB48hpi = moduleGO(genes = Alv_mac_MF_MB_CN_48hpi_MEGENA_modules_GO$Genes, labels = Alv_mac_MF_MB_CN_48hpi_MEGENA_modules_GO$Modules,
universe = rownames(Alv_mac_normalized_counts_MF_MB_CN_48hpi_filtered), pval_GO_cutoff = 1, annotation = "org.Bt.eg.db")

moduleGO_res_TB2hpi = moduleGO(genes = Alv_mac_MF_TB_CN_2hpi_MEGENA_modules_GO$Genes, labels = Alv_mac_MF_TB_CN_2hpi_MEGENA_modules_GO$Modules,
universe = rownames(Alv_mac_normalized_counts_MF_TB_CN_2hpi_filtered), pval_GO_cutoff = 1, annotation = "org.Bt.eg.db")

moduleGO_res_TB24hpi = moduleGO(genes = Alv_mac_MF_TB_CN_24hpi_MEGENA_modules_GO$Genes, labels = Alv_mac_MF_TB_CN_24hpi_MEGENA_modules_GO$Modules,
universe = rownames(Alv_mac_normalized_counts_MF_TB_CN_24hpi_filtered), pval_GO_cutoff = 1, annotation = "org.Bt.eg.db")

moduleGO_res_TB6hpi = moduleGO(genes = Alv_mac_MF_TB_CN_6hpi_MEGENA_modules_GO$Genes, labels = Alv_mac_MF_TB_CN_6hpi_MEGENA_modules_GO$Modules,
universe = rownames(Alv_mac_normalized_counts_MF_TB_CN_6hpi_filtered), pval_GO_cutoff = 1, annotation = "org.Bt.eg.db")

moduleGO_res_TB48hpi = moduleGO(genes = Alv_mac_MF_TB_CN_48hpi_MEGENA_modules_GO$Genes, labels = Alv_mac_MF_TB_CN_48hpi_MEGENA_modules_GO$Modules,
universe = rownames(Alv_mac_normalized_counts_MF_TB_CN_48hpi_filtered), pval_GO_cutoff = 1, annotation = "org.Bt.eg.db")


#Then collate and plot the results 

moduleGO_df_MB2hpi = extractModuleGO(moduleGO_res_MB2hpi)
MB2hpi_module_GO_plot <- plotModuleGO(moduleGO_df_MB2hpi, nTerms = 4, text_size = 8, coord_flip = TRUE)
MB2hpi_module_GO_plot

moduleGO_df_MB24hpi = extractModuleGO(moduleGO_res_MB24hpi)
MB24hpi_module_GO_plot <- plotModuleGO(moduleGO_df_MB24hpi, nTerms = 3, text_size = 10,  adjust = TRUE, coord_flip = FALSE)
MB24hpi_module_GO_plot

moduleGO_df_MB48hpi = extractModuleGO(moduleGO_res_MB48hpi)
MB48hpi_module_GO_plot <- plotModuleGO(moduleGO_df_MB48hpi, nTerms = 3, text_size = 10,  adjust = TRUE, coord_flip = FALSE)
MB48hpi_module_GO_plot

moduleGO_df_MB6hpi = extractModuleGO(moduleGO_res_MB6hpi)
MB6hpi_module_GO_plot <- plotModuleGO(moduleGO_df_MB6hpi, nTerms = 4, text_size = 8, coord_flip = TRUE)
MB6hpi_module_GO_plot
