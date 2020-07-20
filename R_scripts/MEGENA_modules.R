########################################################################################
#################################      MEGENA     ######################################
########################################################################################

#The idea was to use the DGCA code hosted on github to integrate the DGCA results into MEGENA.
#I found however that this always caused errors and there were some hidden objects that 
#I actually need. So below is is the DGCA MEGENA (ddMEGENA) function broken down up until summary.output.
#From there its the same as the MEGENA github code.  

#ddcor_res_perm_500x <- read_csv("/home/workspace/thall/PPDbRNAseqTimeCourse/Count_summarisation/sense/Fresh/correlation_analysis/DGCA_MEGENA/ddcor_res_W1_perm_500x.csv")	#If you run the pipeline from the start, this is not needed. 

#Limits the dataset to a gene p-value threshold
#ddcor_res_sig = ddcor_res_perm_5000x[ddcor_res_perm_5000x$pValDiff < pval_gene_thresh, ]

Alv_mac_normalized_counts_MF_MB_CN_2hpi_ddcor_res_perm_5 = Alv_mac_normalized_counts_MF_MB_CN_2hpi_ddcor_res_perm_5[Alv_mac_normalized_counts_MF_MB_CN_2hpi_ddcor_res_perm_5$empPVals < pval_gene_thresh, ]
Alv_mac_normalized_counts_MF_MB_CN_24hpi_ddcor_res_perm_5 = Alv_mac_normalized_counts_MF_MB_CN_24hpi_ddcor_res_perm_5[Alv_mac_normalized_counts_MF_MB_CN_24hpi_ddcor_res_perm_5$empPVals < pval_gene_thresh, ]
Alv_mac_normalized_counts_MF_MB_CN_6hpi_ddcor_res_perm_5 = Alv_mac_normalized_counts_MF_MB_CN_6hpi_ddcor_res_perm_5[Alv_mac_normalized_counts_MF_MB_CN_6hpi_ddcor_res_perm_5$empPVals < pval_gene_thresh, ]
Alv_mac_normalized_counts_MF_MB_CN_48hpi_ddcor_res_perm_5 = Alv_mac_normalized_counts_MF_MB_CN_48hpi_ddcor_res_perm_5[Alv_mac_normalized_counts_MF_MB_CN_48hpi_ddcor_res_perm_5$empPVals < pval_gene_thresh_48, ]
Alv_mac_normalized_counts_MF_TB_CN_2hpi_ddcor_res_perm_5 = Alv_mac_normalized_counts_MF_TB_CN_2hpi_ddcor_res_perm_5[Alv_mac_normalized_counts_MF_TB_CN_2hpi_ddcor_res_perm_5$empPVals < pval_gene_thresh, ]
Alv_mac_normalized_counts_MF_TB_CN_24hpi_ddcor_res_perm_5 = Alv_mac_normalized_counts_MF_TB_CN_24hpi_ddcor_res_perm_5[Alv_mac_normalized_counts_MF_TB_CN_24hpi_ddcor_res_perm_5$empPVals < pval_gene_thresh, ]
Alv_mac_normalized_counts_MF_TB_CN_6hpi_ddcor_res_perm_5 = Alv_mac_normalized_counts_MF_TB_CN_6hpi_ddcor_res_perm_5[Alv_mac_normalized_counts_MF_TB_CN_6hpi_ddcor_res_perm_5$empPVals < pval_gene_thresh, ]
Alv_mac_normalized_counts_MF_TB_CN_48hpi_ddcor_res_perm_5 = Alv_mac_normalized_counts_MF_TB_CN_48hpi_ddcor_res_perm_5[Alv_mac_normalized_counts_MF_TB_CN_48hpi_ddcor_res_perm_5$empPVals < pval_gene_thresh, ]
#Turns the dataset into 
#ddcor_res_megena = ddcor_res_sig[ , colnames(ddcor_res_sig) %in% c("Gene1", "Gene2", "zScoreDiff")]

Alv_mac_normalized_counts_MF_MB_CN_2hpi_ddcor_res_Megena = Alv_mac_normalized_counts_MF_MB_CN_2hpi_ddcor_res_perm_5[ , colnames(Alv_mac_normalized_counts_MF_MB_CN_2hpi_ddcor_res_perm_5) %in% c("Gene1", "Gene2", "zScoreDiff")]
Alv_mac_normalized_counts_MF_MB_CN_24hpi_ddcor_res_Megena = Alv_mac_normalized_counts_MF_MB_CN_24hpi_ddcor_res_perm_5[ , colnames(Alv_mac_normalized_counts_MF_MB_CN_24hpi_ddcor_res_perm_5) %in% c("Gene1", "Gene2", "zScoreDiff")]
Alv_mac_normalized_counts_MF_MB_CN_6hpi_ddcor_res_Megena = Alv_mac_normalized_counts_MF_MB_CN_6hpi_ddcor_res_perm_5[ , colnames(Alv_mac_normalized_counts_MF_MB_CN_6hpi_ddcor_res_perm_5) %in% c("Gene1", "Gene2", "zScoreDiff")]
Alv_mac_normalized_counts_MF_MB_CN_48hpi_ddcor_res_Megena = Alv_mac_normalized_counts_MF_MB_CN_48hpi_ddcor_res_perm_5[ , colnames(Alv_mac_normalized_counts_MF_MB_CN_48hpi_ddcor_res_perm_5) %in% c("Gene1", "Gene2", "zScoreDiff")]
Alv_mac_normalized_counts_MF_TB_CN_2hpi_ddcor_res_Megena = Alv_mac_normalized_counts_MF_TB_CN_2hpi_ddcor_res_perm_5[ , colnames(Alv_mac_normalized_counts_MF_TB_CN_2hpi_ddcor_res_perm_5) %in% c("Gene1", "Gene2", "zScoreDiff")] 
Alv_mac_normalized_counts_MF_TB_CN_24hpi_ddcor_res_Megena = Alv_mac_normalized_counts_MF_TB_CN_24hpi_ddcor_res_perm_5[ , colnames(Alv_mac_normalized_counts_MF_TB_CN_24hpi_ddcor_res_perm_5) %in% c("Gene1", "Gene2", "zScoreDiff")]
Alv_mac_normalized_counts_MF_TB_CN_6hpi_ddcor_res_Megena = Alv_mac_normalized_counts_MF_TB_CN_6hpi_ddcor_res_perm_5[ , colnames(Alv_mac_normalized_counts_MF_TB_CN_6hpi_ddcor_res_perm_5) %in% c("Gene1", "Gene2", "zScoreDiff")] 
Alv_mac_normalized_counts_MF_TB_CN_48hpi_ddcor_res_Megena = Alv_mac_normalized_counts_MF_TB_CN_48hpi_ddcor_res_perm_5[ , colnames(Alv_mac_normalized_counts_MF_TB_CN_48hpi_ddcor_res_perm_5) %in% c("Gene1", "Gene2", "zScoreDiff")]

#ddcor_res_megena$zScoreDiff = abs(ddcor_res_megena$zScoreDiff)

Alv_mac_normalized_counts_MF_MB_CN_2hpi_ddcor_res_Megena$zScoreDiff = abs(Alv_mac_normalized_counts_MF_MB_CN_2hpi_ddcor_res_Megena$zScoreDiff)
Alv_mac_normalized_counts_MF_MB_CN_24hpi_ddcor_res_Megena$zScoreDiff = abs(Alv_mac_normalized_counts_MF_MB_CN_24hpi_ddcor_res_Megena$zScoreDiff)
Alv_mac_normalized_counts_MF_MB_CN_6hpi_ddcor_res_Megena$zScoreDiff = abs(Alv_mac_normalized_counts_MF_MB_CN_6hpi_ddcor_res_Megena$zScoreDiff)
Alv_mac_normalized_counts_MF_MB_CN_48hpi_ddcor_res_Megena$zScoreDiff = abs(Alv_mac_normalized_counts_MF_MB_CN_48hpi_ddcor_res_Megena$zScoreDiff)
Alv_mac_normalized_counts_MF_TB_CN_2hpi_ddcor_res_Megena$zScoreDiff = abs(Alv_mac_normalized_counts_MF_TB_CN_2hpi_ddcor_res_Megena$zScoreDiff)
Alv_mac_normalized_counts_MF_TB_CN_24hpi_ddcor_res_Megena$zScoreDiff = abs(Alv_mac_normalized_counts_MF_TB_CN_24hpi_ddcor_res_Megena$zScoreDiff)
Alv_mac_normalized_counts_MF_TB_CN_6hpi_ddcor_res_Megena$zScoreDiff = abs(Alv_mac_normalized_counts_MF_TB_CN_6hpi_ddcor_res_Megena$zScoreDiff)
Alv_mac_normalized_counts_MF_TB_CN_48hpi_ddcor_res_Megena$zScoreDiff = abs(Alv_mac_normalized_counts_MF_TB_CN_48hpi_ddcor_res_Megena$zScoreDiff)

#In this step, Planar Filtered Network (PFN) is calculated by taking significant correlation pairs, ijw. 
#In the case of utilizing a different similarity measure, one can independently format the results into 3-column 
#data frame with column names c("row","col","weight"), and make sure the weight column ranges within 0 to 1. 
#Using this as an input to calculate.PFN() will work just as fine.
#pfn_res = MEGENA::calculate.PFN(ddcor_res_megena, doPar = parallelize, num.cores = nCores)

Alv_mac_normalized_counts_MF_MB_CN_2hpi_MEGENA_pfn_res = MEGENA::calculate.PFN(Alv_mac_normalized_counts_MF_MB_CN_2hpi_ddcor_res_Megena, doPar = parallelize, num.cores = nCores)
Alv_mac_normalized_counts_MF_MB_CN_24hpi_MEGENA_pfn_res = MEGENA::calculate.PFN(Alv_mac_normalized_counts_MF_MB_CN_24hpi_ddcor_res_Megena, doPar = parallelize, num.cores = nCores)
Alv_mac_normalized_counts_MF_MB_CN_6hpi_MEGENA_pfn_res = MEGENA::calculate.PFN(Alv_mac_normalized_counts_MF_MB_CN_6hpi_ddcor_res_Megena, doPar = parallelize, num.cores = nCores)
Alv_mac_normalized_counts_MF_MB_CN_48hpi_MEGENA_pfn_res = MEGENA::calculate.PFN(Alv_mac_normalized_counts_MF_MB_CN_48hpi_ddcor_res_Megena, doPar = parallelize, num.cores = nCores)
Alv_mac_normalized_counts_MF_TB_CN_2hpi_MEGENA_pfn_res = MEGENA::calculate.PFN(Alv_mac_normalized_counts_MF_TB_CN_2hpi_ddcor_res_Megena, doPar = parallelize, num.cores = nCores)
Alv_mac_normalized_counts_MF_TB_CN_24hpi_MEGENA_pfn_res = MEGENA::calculate.PFN(Alv_mac_normalized_counts_MF_TB_CN_24hpi_ddcor_res_Megena, doPar = parallelize, num.cores = nCores)
Alv_mac_normalized_counts_MF_TB_CN_6hpi_MEGENA_pfn_res = MEGENA::calculate.PFN(Alv_mac_normalized_counts_MF_TB_CN_6hpi_ddcor_res_Megena, doPar = parallelize, num.cores = nCores)
Alv_mac_normalized_counts_MF_TB_CN_48hpi_MEGENA_pfn_res = MEGENA::calculate.PFN(Alv_mac_normalized_counts_MF_TB_CN_48hpi_ddcor_res_Megena, doPar = parallelize, num.cores = nCores)

#normalize to follow the convention that weights are in [0, 1]
#pfn_res$weight = (pfn_res$weight/max(pfn_res$weight)) * 0.999999999

Alv_mac_normalized_counts_MF_MB_CN_2hpi_MEGENA_pfn_res$weight = (Alv_mac_normalized_counts_MF_MB_CN_2hpi_MEGENA_pfn_res$weight/max(Alv_mac_normalized_counts_MF_MB_CN_2hpi_MEGENA_pfn_res$weight)) * 0.999999999
Alv_mac_normalized_counts_MF_MB_CN_24hpi_MEGENA_pfn_res$weight = (Alv_mac_normalized_counts_MF_MB_CN_24hpi_MEGENA_pfn_res$weight/max(Alv_mac_normalized_counts_MF_MB_CN_24hpi_MEGENA_pfn_res$weight)) * 0.999999999
Alv_mac_normalized_counts_MF_MB_CN_6hpi_MEGENA_pfn_res$weight = (Alv_mac_normalized_counts_MF_MB_CN_6hpi_MEGENA_pfn_res$weight/max(Alv_mac_normalized_counts_MF_MB_CN_6hpi_MEGENA_pfn_res$weight)) * 0.999999999
Alv_mac_normalized_counts_MF_MB_CN_48hpi_MEGENA_pfn_res$weight = (Alv_mac_normalized_counts_MF_MB_CN_48hpi_MEGENA_pfn_res$weight/max(Alv_mac_normalized_counts_MF_MB_CN_48hpi_MEGENA_pfn_res$weight)) * 0.999999999
Alv_mac_normalized_counts_MF_TB_CN_2hpi_MEGENA_pfn_res$weight = (Alv_mac_normalized_counts_MF_TB_CN_2hpi_MEGENA_pfn_res$weight/max(Alv_mac_normalized_counts_MF_TB_CN_2hpi_MEGENA_pfn_res$weight)) * 0.999999999
Alv_mac_normalized_counts_MF_TB_CN_24hpi_MEGENA_pfn_res$weight = (Alv_mac_normalized_counts_MF_TB_CN_24hpi_MEGENA_pfn_res$weight/max(Alv_mac_normalized_counts_MF_TB_CN_24hpi_MEGENA_pfn_res$weight)) * 0.999999999
Alv_mac_normalized_counts_MF_TB_CN_6hpi_MEGENA_pfn_res$weight = (Alv_mac_normalized_counts_MF_TB_CN_6hpi_MEGENA_pfn_res$weight/max(Alv_mac_normalized_counts_MF_TB_CN_6hpi_MEGENA_pfn_res$weight)) * 0.999999999
Alv_mac_normalized_counts_MF_TB_CN_48hpi_MEGENA_pfn_res$weight = (Alv_mac_normalized_counts_MF_TB_CN_48hpi_MEGENA_pfn_res$weight/max(Alv_mac_normalized_counts_MF_TB_CN_48hpi_MEGENA_pfn_res$weight)) * 0.999999999

#Create an igraph object of the the pfn
#g = igraph::graph.data.frame(pfn_res, directed = FALSE)

Alv_mac_MF_MB_CN_2hpi_MEGENA_igraph = igraph::graph.data.frame(Alv_mac_normalized_counts_MF_MB_CN_2hpi_MEGENA_pfn_res, directed = FALSE) 
Alv_mac_MF_MB_CN_24hpi_MEGENA_igraph = igraph::graph.data.frame(Alv_mac_normalized_counts_MF_MB_CN_24hpi_MEGENA_pfn_res, directed = FALSE)
Alv_mac_MF_MB_CN_6hpi_MEGENA_igraph = igraph::graph.data.frame(Alv_mac_normalized_counts_MF_MB_CN_6hpi_MEGENA_pfn_res, directed = FALSE) 
Alv_mac_MF_MB_CN_48hpi_MEGENA_igraph = igraph::graph.data.frame(Alv_mac_normalized_counts_MF_MB_CN_48hpi_MEGENA_pfn_res, directed = FALSE)
Alv_mac_MF_TB_CN_2hpi_MEGENA_igraph = igraph::graph.data.frame(Alv_mac_normalized_counts_MF_TB_CN_2hpi_MEGENA_pfn_res, directed = FALSE) 
Alv_mac_MF_TB_CN_24hpi_MEGENA_igraph = igraph::graph.data.frame(Alv_mac_normalized_counts_MF_TB_CN_24hpi_MEGENA_pfn_res, directed = FALSE)
Alv_mac_MF_TB_CN_6hpi_MEGENA_igraph = igraph::graph.data.frame(Alv_mac_normalized_counts_MF_TB_CN_6hpi_MEGENA_pfn_res, directed = FALSE) 
Alv_mac_MF_TB_CN_48hpi_MEGENA_igraph = igraph::graph.data.frame(Alv_mac_normalized_counts_MF_TB_CN_48hpi_MEGENA_pfn_res, directed = FALSE)


#MCA clustering is performed to identify multiscale clustering analysis. “MEGENA.output”" is the core output to be used in the down-stream analyses for summarization and plotting of sub-modules#
#I have changed the min mod size to 20, as i discard these anyway but it also interferes with the GO enrichment later. 
Alv_mac_MF_MB_CN_2hpi_MEGENA_output = MEGENA::do.MEGENA(Alv_mac_MF_MB_CN_2hpi_MEGENA_igraph, mod.pval = modulePVal, hub.pval = hubPVal,
                                  remove.unsig = TRUE, min.size = minModSize, max.size = maxModSize,
                                  doPar = parallelize, num.cores = nCores, n.perm = nPerm, save.output = saveOutput)
Alv_mac_MF_MB_CN_24hpi_MEGENA_output = MEGENA::do.MEGENA(Alv_mac_MF_MB_CN_24hpi_MEGENA_igraph, mod.pval = module.pval, hub.pval = hubPVal,
                                  remove.unsig = TRUE, min.size = minModSize, max.size = maxModSize,
                                  doPar = parallelize, num.cores = nCores, n.perm = nPerm, save.output = saveOutput)
Alv_mac_MF_MB_CN_6hpi_MEGENA_output = MEGENA::do.MEGENA(Alv_mac_MF_MB_CN_6hpi_MEGENA_igraph, mod.pval = module.pval, hub.pval = hubPVal,
                                  remove.unsig = TRUE, min.size = minModSize, max.size = maxModSize,
                                  doPar = parallelize, num.cores = nCores, n.perm = nPerm, save.output = saveOutput)
Alv_mac_MF_MB_CN_48hpi_MEGENA_output = MEGENA::do.MEGENA(Alv_mac_MF_MB_CN_48hpi_MEGENA_igraph, mod.pval = 0.1, hub.pval = 0.1,
                                  remove.unsig = TRUE, min.size = minModSize, max.size = maxModSize,
                                  doPar = parallelize, num.cores = nCores, n.perm = nPerm, save.output = saveOutput)
Alv_mac_MF_TB_CN_2hpi_MEGENA_output = MEGENA::do.MEGENA(Alv_mac_MF_TB_CN_2hpi_MEGENA_igraph, mod.pval = module.pval, hub.pval = hubPVal,
                                  remove.unsig = TRUE, min.size = minModSize, max.size = maxModSize,
                                  doPar = parallelize, num.cores = nCores, n.perm = nPerm, save.output = saveOutput)
Alv_mac_MF_TB_CN_24hpi_MEGENA_output = MEGENA::do.MEGENA(Alv_mac_MF_TB_CN_24hpi_MEGENA_igraph, mod.pval = module.pval, hub.pval = hubPVal,
                                  remove.unsig = TRUE, min.size = minModSize, max.size = maxModSize,
                                  doPar = parallelize, num.cores = nCores, n.perm = nPerm, save.output = saveOutput)
Alv_mac_MF_TB_CN_6hpi_MEGENA_output = MEGENA::do.MEGENA(Alv_mac_MF_TB_CN_6hpi_MEGENA_igraph, mod.pval = module.pval, hub.pval = hubPVal,
                                  remove.unsig = TRUE, min.size = minModSize, max.size = maxModSize,
                                  doPar = parallelize, num.cores = nCores, n.perm = nPerm, save.output = saveOutput)
Alv_mac_MF_TB_CN_48hpi_MEGENA_output = MEGENA::do.MEGENA(Alv_mac_MF_TB_CN_48hpi_MEGENA_igraph, mod.pval = module.pval, hub.pval = hubPVal,
                                  remove.unsig = TRUE, min.size = minModSize, max.size = maxModSize,
                                  doPar = parallelize, num.cores = nCores, n.perm = nPerm, save.output = saveOutput)

#I am not sure I actually need this. 			


Alv_mac_MF_MB_CN_2hpi_MEGENA_modules = Alv_mac_MF_MB_CN_2hpi_MEGENA_output$module.output$modules
Alv_mac_MF_MB_CN_2hpi_MEGENA_modules_df = data.frame(Genes = unlist(Alv_mac_MF_MB_CN_2hpi_MEGENA_modules),
                                          Modules = rep(names(Alv_mac_MF_MB_CN_2hpi_MEGENA_modules), sapply(Alv_mac_MF_MB_CN_2hpi_MEGENA_modules, length)))
Alv_mac_MF_MB_CN_2hpi_megena_modules_output = list(modules = Alv_mac_MF_MB_CN_2hpi_MEGENA_modules_df, full = Alv_mac_MF_MB_CN_2hpi_MEGENA_output)

Alv_mac_MF_MB_CN_24hpi_MEGENA_modules = Alv_mac_MF_MB_CN_24hpi_MEGENA_output$module.output$modules
Alv_mac_MF_MB_CN_24hpi_MEGENA_modules_df = data.frame(Genes = unlist(Alv_mac_MF_MB_CN_24hpi_MEGENA_modules),
                                                     Modules = rep(names(Alv_mac_MF_MB_CN_24hpi_MEGENA_modules), sapply(Alv_mac_MF_MB_CN_24hpi_MEGENA_modules, length)))
Alv_mac_MF_MB_CN_24hpi_megena_modules_output = list(modules = Alv_mac_MF_MB_CN_24hpi_MEGENA_modules_df, full = Alv_mac_MF_MB_CN_24hpi_MEGENA_output)

Alv_mac_MF_MB_CN_6hpi_MEGENA_modules = Alv_mac_MF_MB_CN_6hpi_MEGENA_output$module.output$modules
Alv_mac_MF_MB_CN_6hpi_MEGENA_modules_df = data.frame(Genes = unlist(Alv_mac_MF_MB_CN_6hpi_MEGENA_modules),
                                                     Modules = rep(names(Alv_mac_MF_MB_CN_6hpi_MEGENA_modules), sapply(Alv_mac_MF_MB_CN_6hpi_MEGENA_modules, length)))
Alv_mac_MF_MB_CN_6hpi_megena_modules_output = list(modules = Alv_mac_MF_MB_CN_6hpi_MEGENA_modules_df, full = Alv_mac_MF_MB_CN_6hpi_MEGENA_output)

Alv_mac_MF_MB_CN_48hpi_MEGENA_modules = Alv_mac_MF_MB_CN_48hpi_MEGENA_output$module.output$modules
Alv_mac_MF_MB_CN_48hpi_MEGENA_modules_df = data.frame(Genes = unlist(Alv_mac_MF_MB_CN_48hpi_MEGENA_modules),
                                                     Modules = rep(names(Alv_mac_MF_MB_CN_48hpi_MEGENA_modules), sapply(Alv_mac_MF_MB_CN_48hpi_MEGENA_modules, length)))
Alv_mac_MF_MB_CN_48hpi_megena_modules_output = list(modules = Alv_mac_MF_MB_CN_48hpi_MEGENA_modules_df, full = Alv_mac_MF_MB_CN_48hpi_MEGENA_output)

Alv_mac_MF_TB_CN_2hpi_MEGENA_modules = Alv_mac_MF_TB_CN_2hpi_MEGENA_output$module.output$modules
Alv_mac_MF_TB_CN_2hpi_MEGENA_modules_df = data.frame(Genes = unlist(Alv_mac_MF_TB_CN_2hpi_MEGENA_modules),
                                                     Modules = rep(names(Alv_mac_MF_TB_CN_2hpi_MEGENA_modules), sapply(Alv_mac_MF_TB_CN_2hpi_MEGENA_modules, length)))
Alv_mac_MF_TB_CN_2hpi_megena_modules_output = list(modules = Alv_mac_MF_TB_CN_2hpi_MEGENA_modules_df, full = Alv_mac_MF_TB_CN_2hpi_MEGENA_output)

Alv_mac_MF_TB_CN_24hpi_MEGENA_modules = Alv_mac_MF_TB_CN_24hpi_MEGENA_output$module.output$modules
Alv_mac_MF_TB_CN_24hpi_MEGENA_modules_df = data.frame(Genes = unlist(Alv_mac_MF_TB_CN_24hpi_MEGENA_modules),
                                                     Modules = rep(names(Alv_mac_MF_TB_CN_24hpi_MEGENA_modules), sapply(Alv_mac_MF_TB_CN_24hpi_MEGENA_modules, length)))
Alv_mac_MF_TB_CN_24hpi_megena_modules_output = list(modules = Alv_mac_MF_TB_CN_24hpi_MEGENA_modules_df, full = Alv_mac_MF_TB_CN_24hpi_MEGENA_output)

Alv_mac_MF_TB_CN_6hpi_MEGENA_modules = Alv_mac_MF_TB_CN_6hpi_MEGENA_output$module.output$modules
Alv_mac_MF_TB_CN_6hpi_MEGENA_modules_df = data.frame(Genes = unlist(Alv_mac_MF_TB_CN_6hpi_MEGENA_modules),
                                                     Modules = rep(names(Alv_mac_MF_TB_CN_6hpi_MEGENA_modules), sapply(Alv_mac_MF_TB_CN_6hpi_MEGENA_modules, length)))
Alv_mac_MF_TB_CN_6hpi_megena_modules_output = list(modules = Alv_mac_MF_TB_CN_6hpi_MEGENA_modules_df, full = Alv_mac_MF_TB_CN_6hpi_MEGENA_output)

Alv_mac_MF_TB_CN_48hpi_MEGENA_modules = Alv_mac_MF_TB_CN_48hpi_MEGENA_output$module.output$modules
Alv_mac_MF_TB_CN_48hpi_MEGENA_modules_df = data.frame(Genes = unlist(Alv_mac_MF_TB_CN_48hpi_MEGENA_modules),
                                                     Modules = rep(names(Alv_mac_MF_TB_CN_48hpi_MEGENA_modules), sapply(Alv_mac_MF_TB_CN_48hpi_MEGENA_modules, length)))
Alv_mac_MF_TB_CN_48hpi_megena_modules_output = list(modules = Alv_mac_MF_TB_CN_48hpi_MEGENA_modules_df, full = Alv_mac_MF_TB_CN_48hpi_MEGENA_output)



Alv_mac_MF_MB_CN_2hpi_MEGENA_summary.output <- MEGENA.ModuleSummary(Alv_mac_MF_MB_CN_2hpi_MEGENA_output,
                                                                    mod.pvalue = module.pval,hub.pvalue = hub.pval,
                                                                    min.size = 20,max.size = 1000,
                                                                    annot.table = annot.table,id.col = id.col,symbol.col = symbol.col,
                                                                    output.sig = TRUE)

Alv_mac_MF_MB_CN_24hpi_MEGENA_summary.output <- MEGENA.ModuleSummary(Alv_mac_MF_MB_CN_24hpi_MEGENA_output,
                                                                    mod.pvalue = module.pval,hub.pvalue = hub.pval,
                                                                    min.size = 20,max.size = 1000,
                                                                    annot.table = annot.table,id.col = id.col,symbol.col = symbol.col,
                                                                    output.sig = TRUE)

Alv_mac_MF_MB_CN_6hpi_MEGENA_summary.output <- MEGENA.ModuleSummary(Alv_mac_MF_MB_CN_6hpi_MEGENA_output,
                                                                    mod.pvalue = module.pval,hub.pvalue = hub.pval,
                                                                    min.size = 20,max.size = 1000,
                                                                    annot.table = annot.table,id.col = id.col,symbol.col = symbol.col,
                                                                    output.sig = TRUE)

Alv_mac_MF_MB_CN_48hpi_MEGENA_summary.output <- MEGENA.ModuleSummary(Alv_mac_MF_MB_CN_48hpi_MEGENA_output,
                                                                    mod.pvalue = modulePVal,hub.pvalue = 0.1,
                                                                    min.size = 10,max.size = 1000,
                                                                    annot.table = annot.table,id.col = id.col,symbol.col = symbol.col,
                                                                    output.sig = TRUE)

Alv_mac_MF_TB_CN_2hpi_MEGENA_summary.output <- MEGENA.ModuleSummary(Alv_mac_MF_TB_CN_2hpi_MEGENA_output,
                                                                    mod.pvalue = module.pval,hub.pvalue = hub.pval,
                                                                    min.size = 20,max.size = 1000,
                                                                    annot.table = annot.table,id.col = id.col,symbol.col = symbol.col,
                                                                    output.sig = TRUE)

Alv_mac_MF_TB_CN_24hpi_MEGENA_summary.output <- MEGENA.ModuleSummary(Alv_mac_MF_TB_CN_24hpi_MEGENA_output,
                                                                    mod.pvalue = module.pval,hub.pvalue = hub.pval,
                                                                    min.size = 20,max.size = 1000,
                                                                    annot.table = annot.table,id.col = id.col,symbol.col = symbol.col,
                                                                    output.sig = TRUE)

Alv_mac_MF_TB_CN_6hpi_MEGENA_summary.output <- MEGENA.ModuleSummary(Alv_mac_MF_TB_CN_6hpi_MEGENA_output,
                                                                    mod.pvalue = module.pval,hub.pvalue = hub.pval,
                                                                    min.size = 20,max.size = 1000,
                                                                    annot.table = annot.table,id.col = id.col,symbol.col = symbol.col,
                                                                    output.sig = TRUE)

Alv_mac_MF_TB_CN_48hpi_MEGENA_summary.output <- MEGENA.ModuleSummary(Alv_mac_MF_TB_CN_48hpi_MEGENA_output,
                                                                    mod.pvalue = module.pval,hub.pvalue = hub.pval,
                                                                    min.size = 20,max.size = 1000,
                                                                    annot.table = annot.table,id.col = id.col,symbol.col = symbol.col,
                                                                    output.sig = TRUE)

#plotting module hierarchy

#module.table <- summary.output$module.table
#colnames(module.table)[1] <- "id" # first column of module table must be labelled as "id".

#hierarchy.obj <- plot_module_hierarchy(module.table = module.table,label.scaleFactor = 0.15,
#                                       arrow.size = 0.03,node.label.color = "blue")

Alv_mac_MF_MB_CN_2hpi_MEGENA_module_table <- Alv_mac_MF_MB_CN_2hpi_MEGENA_summary.output$module.table
colnames(Alv_mac_MF_MB_CN_2hpi_MEGENA_module_table)[1] <- "id" # first column of module table must be labelled as "id".
Alv_mac_MF_MB_CN_2hpi_MEGENA_hierarchy.obj <- plot_module_hierarchy(module.table = Alv_mac_MF_MB_CN_2hpi_MEGENA_module_table, label.scaleFactor = 0.15,
                                       arrow.size = 0.03,node.label.color = "blue")
#X11();
Alv_mac_MF_MB_CN_2hpi_MEGENA_hierarchy.plot <- print(Alv_mac_MF_MB_CN_2hpi_MEGENA_hierarchy.obj[[1]])

Alv_mac_MF_MB_CN_24hpi_MEGENA_module_table <- Alv_mac_MF_MB_CN_24hpi_MEGENA_summary.output$module.table
colnames(Alv_mac_MF_MB_CN_24hpi_MEGENA_module_table)[1] <- "id" # first column of module table must be labelled as "id".
Alv_mac_MF_MB_CN_24hpi_MEGENA_hierarchy.obj <- plot_module_hierarchy(module.table = Alv_mac_MF_MB_CN_24hpi_MEGENA_module_table, label.scaleFactor = 0.15,
                                                                    arrow.size = 0.03,node.label.color = "#808080")
#X11();
Alv_mac_MF_MB_CN_24hpi_MEGENA_hierarchy.plot <- print(Alv_mac_MF_MB_CN_24hpi_MEGENA_hierarchy.obj[[1]])

Alv_mac_MF_MB_CN_6hpi_MEGENA_module_table <- Alv_mac_MF_MB_CN_6hpi_MEGENA_summary.output$module.table
colnames(Alv_mac_MF_MB_CN_6hpi_MEGENA_module_table)[1] <- "id" # first column of module table must be labelled as "id".
Alv_mac_MF_MB_CN_6hpi_MEGENA_hierarchy.obj <- plot_module_hierarchy(module.table = Alv_mac_MF_MB_CN_6hpi_MEGENA_module_table, label.scaleFactor = 0.15,
                                                                    arrow.size = 0.03,node.label.color = "blue")
#X11();
Alv_mac_MF_MB_CN_6hpi_MEGENA_hierarchy.plot <- print(Alv_mac_MF_MB_CN_6hpi_MEGENA_hierarchy.obj[[1]])

Alv_mac_MF_MB_CN_48hpi_MEGENA_module_table <- Alv_mac_MF_MB_CN_48hpi_MEGENA_summary.output$module.table
colnames(Alv_mac_MF_MB_CN_48hpi_MEGENA_module_table)[1] <- "id" # first column of module table must be labelled as "id".
Alv_mac_MF_MB_CN_48hpi_MEGENA_hierarchy.obj <- plot_module_hierarchy(module.table = Alv_mac_MF_MB_CN_48hpi_MEGENA_module_table, label.scaleFactor = 0.15,
                                                                    arrow.size = 0.03,node.label.color = "blue")
#X11();
Alv_mac_MF_MB_CN_48hpi_MEGENA_hierarchy.plot <- print(Alv_mac_MF_MB_CN_48hpi_MEGENA_hierarchy.obj[[1]])

Alv_mac_MF_TB_CN_2hpi_MEGENA_module_table <- Alv_mac_MF_TB_CN_2hpi_MEGENA_summary.output$module.table
colnames(Alv_mac_MF_TB_CN_2hpi_MEGENA_module_table)[1] <- "id" # first column of module table must be labelled as "id".
Alv_mac_MF_TB_CN_2hpi_MEGENA_hierarchy.obj <- plot_module_hierarchy(module.table = Alv_mac_MF_TB_CN_2hpi_MEGENA_module_table, label.scaleFactor = 0.15,
                                                                    arrow.size = 0.03,node.label.color = "blue")
#X11();
Alv_mac_MF_TB_CN_2hpi_MEGENA_hierarchy.plot <- print(Alv_mac_MF_TB_CN_2hpi_MEGENA_hierarchy.obj[[1]])

Alv_mac_MF_TB_CN_24hpi_MEGENA_module_table <- Alv_mac_MF_TB_CN_24hpi_MEGENA_summary.output$module.table
colnames(Alv_mac_MF_TB_CN_24hpi_MEGENA_module_table)[1] <- "id" # first column of module table must be labelled as "id".
Alv_mac_MF_TB_CN_24hpi_MEGENA_hierarchy.obj <- plot_module_hierarchy(module.table = Alv_mac_MF_TB_CN_24hpi_MEGENA_module_table, label.scaleFactor = 0.15,
                                                                    arrow.size = 0.03,node.label.color = "blue")
#X11();
Alv_mac_MF_TB_CN_24hpi_MEGENA_hierarchy.plot <- print(Alv_mac_MF_TB_CN_24hpi_MEGENA_hierarchy.obj[[1]])

Alv_mac_MF_TB_CN_6hpi_MEGENA_module_table <- Alv_mac_MF_TB_CN_6hpi_MEGENA_summary.output$module.table
colnames(Alv_mac_MF_TB_CN_6hpi_MEGENA_module_table)[1] <- "id" # first column of module table must be labelled as "id".
Alv_mac_MF_TB_CN_6hpi_MEGENA_hierarchy.obj <- plot_module_hierarchy(module.table = Alv_mac_MF_TB_CN_6hpi_MEGENA_module_table, label.scaleFactor = 0.15,
                                                                    arrow.size = 0.03,node.label.color = "blue")
#X11();
Alv_mac_MF_TB_CN_6hpi_MEGENA_hierarchy.plot <- print(Alv_mac_MF_TB_CN_6hpi_MEGENA_hierarchy.obj[[1]])

Alv_mac_MF_TB_CN_48hpi_MEGENA_module_table <- Alv_mac_MF_TB_CN_48hpi_MEGENA_summary.output$module.table
colnames(Alv_mac_MF_TB_CN_48hpi_MEGENA_module_table)[1] <- "id" # first column of module table must be labelled as "id".
Alv_mac_MF_TB_CN_48hpi_MEGENA_hierarchy.obj <- plot_module_hierarchy(module.table = Alv_mac_MF_TB_CN_48hpi_MEGENA_module_table, label.scaleFactor = 0.15,
                                                                    arrow.size = 0.03,node.label.color = "blue")
#X11();
Alv_mac_MF_TB_CN_48hpi_MEGENA_hierarchy.plot <- print(Alv_mac_MF_TB_CN_48hpi_MEGENA_hierarchy.obj[[1]])


#write module summary table to .csv
write.csv(Alv_mac_MF_MB_CN_2hpi_MEGENA_module_table, file="Alv_mac_MF_MB_CN_2hpi_MEGENA_module_table.csv")
write.csv(Alv_mac_MF_MB_CN_24hpi_MEGENA_module_table, file="Alv_mac_MF_MB_CN_24hpi_MEGENA_module_table.csv")
write.csv(Alv_mac_MF_MB_CN_6hpi_MEGENA_module_table, file="Alv_mac_MF_MB_CN_6hpi_MEGENA_module_table.csv")
write.csv(Alv_mac_MF_MB_CN_48hpi_MEGENA_module_table, file="Alv_mac_MF_MB_CN_48hpi_MEGENA_module_table.csv")
write.csv(Alv_mac_MF_TB_CN_2hpi_MEGENA_module_table, file="Alv_mac_MF_TB_CN_2hpi_MEGENA_module_table.csv")
write.csv(Alv_mac_MF_TB_CN_24hpi_MEGENA_module_table, file="Alv_mac_MF_TB_CN_24hpi_MEGENA_module_table.csv")
write.csv(Alv_mac_MF_TB_CN_6hpi_MEGENA_module_table, file="Alv_mac_MF_TB_CN_6hpi_MEGENA_module_table.csv")
write.csv(Alv_mac_MF_TB_CN_48hpi_MEGENA_module_table, file="Alv_mac_MF_TB_CN_48hpi_MEGENA_module_table.csv")

#3 column empval adjusted for networks

write.csv(Alv_mac_normalized_counts_MF_MB_CN_2hpi_MEGENA_pfn_res, file="Alv_mac_network_table_MB2hpi.csv")
write.csv(Alv_mac_normalized_counts_MF_MB_CN_24hpi_MEGENA_pfn_res, file="Alv_mac_network_table_MB24hpi.csv")
write.csv(Alv_mac_normalized_counts_MF_MB_CN_6hpi_MEGENA_pfn_res, file="Alv_mac_network_table_MB6hpi.csv")
write.csv(Alv_mac_normalized_counts_MF_MB_CN_48hpi_MEGENA_pfn_res, file="Alv_mac_network_table_MB48hpi.csv")
write.csv(Alv_mac_normalized_counts_MF_TB_CN_2hpi_MEGENA_pfn_res, file="Alv_mac_network_table_TB2hpi.csv")
write.csv(Alv_mac_normalized_counts_MF_TB_CN_24hpi_MEGENA_pfn_res, file="Alv_mac_network_table_TB24hpi.csv")
write.csv(Alv_mac_normalized_counts_MF_TB_CN_6hpi_MEGENA_pfn_res, file="Alv_mac_network_table_TB6hpi.csv")
write.csv(Alv_mac_normalized_counts_MF_TB_CN_48hpi_MEGENA_pfn_res, file="Alv_mac_network_table_TB48hpi.csv")

#write specific module to .csv.
head(Alv_mac_MF_MB_CN_2hpi_MEGENA_summary.output$modules)
Alv_mac_MF_MB_CN_2hpi_c1_3 <- Alv_mac_MF_MB_CN_2hpi_MEGENA_summary.output$modules$c1_3 #or whatever module you want from Alv_mac_MF_TB_CN_48hpi_MEGENA_summary.output$modules$module.table 
write.csv(Alv_mac_MF_MB_CN_2hpi_c1_3, file = "Alv_mac_MF_MB_CN_2hpi_c1_3.csv")

#What I did then was cat all the module files using *C1* as the identifier. Removed "GENE_I.D:" and "BGD:". Used biomart to annotate and remove LOC's. 

######################## plotting a specific sub-module. #############################

#output.summary 		output from summary function, "MEGENA.ModuleSummary".
#PFN 					igraph object retaining PFN topology.
#subset.module 			A character vector for list of module names to plot. Default = NULL plots all modules in output.summary
#col.names 				A character vector for list of colors to be used for coloring children modules.
#gene.set 				A list object containing signatures for customized coloring of nodes in resulting network plot.
#color.code 			A character vector with matched length to "gene.set", to specify colors for each signature.
#label.hubs.only		TRUE/FALSE to show labels for significant hub genes only, or all genes. Default is TRUE.
#hubLabel.col 			Label color for hubs. Default is "red"
#show.legend 			TRUE/FALSE for showing node legend on the bottom of the figure.
#hubLabel.sizeProp		A multiplicative factor to adjust hub label sizes with respect to node size values. Default is 0.5
#show.topn.hubs 		Maximal number of hubs to label on module subnetwork. Default is 10.
#node.sizeProp 			A multiplicative factor to adjust node sizes with respect to 90th percentile degree node size. Default is 13
#label.size				Prop A multiplicative factor to adjust node label sizes with respect to 90th percentile degree node size. Default is 13
#label.scaleFactor		Overall scale factor to control the final size of node labels appearing in figure. Default is 10. 
#label.alpha 			Transparency value ranging from 0 (transparent) to 1 (solid). Default is 0.5.
#layout 				Network layout algorithm to apply. Options are: "kamada.kawai", "fruchterman.reingold".
#output.plot 			logical value. output.plot = TRUE generates figure files under folder, "modulePlot".
#out.dir 				if output.plot = TRUE, then out.dir is created and resulting figures are exported to .png files to the folder


MB_CN_24_sub_module_c1_6 <- plot_module(output.summary = Alv_mac_MF_MB_CN_24hpi_MEGENA_summary.output, PFN = Alv_mac_MF_MB_CN_24hpi_MEGENA_igraph, subset.module = "c1_6",
            layout = "kamada.kawai", label.hubs.only = FALSE,
            gene.set = NULL, color.code =  c("#ff3333", "grey"),
            output.plot = FALSE, out.dir = "modulePlot", col.names = c("#ff3333"),
            hubLabel.col = "black",  show.topn.hubs = 10, show.legend = FALSE, ,hubLabel.sizeProp = 0.5, node.sizeProp = 13)
			
MB_CN_24_sub_module_c1_6


MB_CN_48_sub_modules <- plot_module(output.summary = Alv_mac_MF_MB_CN_48hpi_MEGENA_summary.output, PFN = Alv_mac_MF_MB_CN_48hpi_MEGENA_igraph,
            layout = "fruchterman.reingold", label.hubs.only =  FALSE,
            gene.set = NULL, color.code =  c("#4d94ff", "grey"),
            output.plot = FALSE, out.dir = "modulePlot", label.scaleFactor = 3,
            hubLabel.col = "#4d94ff", hubLabel.sizeProp = 6, show.topn.hubs = 10, show.legend = FALSE)
			
MB_CN_48_sub_modules
