########################################################################################
##################################      DGCA     #######################################
########################################################################################
setwd("/home/workspace/thall/networks/correlation_networks/Alv_mac/network_output")

#Load in the normalized counts, first samples are control, second set is the infected
#In the case of Kirstens data, the only thing that has to be changed here is the Wx number. 
Alv_mac_normalized_counts_MF<- read.csv("Alv_mac_normalized_counts_ARD.csv", header = TRUE, row.names = 1)

#The design matrix is created based on data from a separate file that maps sample names to cell types.
design_mat <- read_csv("Design_DGCA_Alv_mac.csv")

#We can map the ID's here if not done so before
rownames <- as.data.frame(mapIds(org.Bt.eg.db, row.names(Alv_mac_normalized_counts_MF), 'SYMBOL', 'ENTREZID'))
row.names(Alv_mac_normalized_counts_MF) <- make.names(rownames$`mapIds(org.Bt.eg.db, row.names(Alv_mac_normalized_counts_MF), "SYMBOL", "ENTREZID")`, unique = TRUE)

#Here I am setting the control and infected to an integer of 10, and assing them to each of the samples in normalized dataset.
#This replaces the need to make a desgin matrix in excell as above. 
control_samples = 10; infected_samples = 10 
cell_type = c(rep("control", control_samples), rep("infected", infected_samples))
colnames(design_mat) = c("control", "infected")
str(design_mat)
design_mat = makeDesign(cell_type)

#48 hours is missing a sample, so the matrix must reflect this
control_samples = 9; infected_samples = 9 
cell_type = c(rep("control", control_samples), rep("infected", infected_samples))
colnames(design_mat) = c("control", "infected")
str(design_mat)
design_mat_48hpi = makeDesign(cell_type)


#before filtration, sometimes  the data needs a bit of sorting. With the Alv_mac data, I needed to sort TB and MB into groups. I did this by exporting columns from the master counts file. 
Alv_mac_normalized_counts_MF_MB_CN <- Alv_mac_normalized_counts_MF[,-grep("TB",colnames(Alv_mac_normalized_counts_MF))] #This deletes TB samples 
Alv_mac_normalized_counts_MF_TB_CN <- Alv_mac_normalized_counts_MF[,-grep("MB",colnames(Alv_mac_normalized_counts_MF))] #This deletes MB samples 
#and then seperated into timepoints 
Alv_mac_normalized_counts_MF_MB_CN_2hpi <- Alv_mac_normalized_counts_MF_MB_CN[,grep("2H",colnames(Alv_mac_normalized_counts_MF_MB_CN))] 
Alv_mac_normalized_counts_MF_MB_CN_6hpi <- Alv_mac_normalized_counts_MF_MB_CN[,grep("6H",colnames(Alv_mac_normalized_counts_MF_MB_CN))] 
Alv_mac_normalized_counts_MF_MB_CN_24hpi <- Alv_mac_normalized_counts_MF_MB_CN[,grep("24H",colnames(Alv_mac_normalized_counts_MF_MB_CN))] 
Alv_mac_normalized_counts_MF_MB_CN_48hpi <- Alv_mac_normalized_counts_MF_MB_CN[,grep("48H",colnames(Alv_mac_normalized_counts_MF_MB_CN))] 
Alv_mac_normalized_counts_MF_TB_CN_2hpi <- Alv_mac_normalized_counts_MF_TB_CN[,grep("2H",colnames(Alv_mac_normalized_counts_MF_TB_CN))] 
Alv_mac_normalized_counts_MF_TB_CN_6hpi <- Alv_mac_normalized_counts_MF_TB_CN[,grep("6H",colnames(Alv_mac_normalized_counts_MF_TB_CN))] 
Alv_mac_normalized_counts_MF_TB_CN_24hpi <- Alv_mac_normalized_counts_MF_TB_CN[,grep("24H",colnames(Alv_mac_normalized_counts_MF_TB_CN))] 
Alv_mac_normalized_counts_MF_TB_CN_48hpi <- Alv_mac_normalized_counts_MF_TB_CN[,grep("48H",colnames(Alv_mac_normalized_counts_MF_TB_CN))] 



#It is often the case that the genes with the lowest average expression levels, or the lowest variance in expression levels,
#are less likely to have interesting and/or biologically relevant differences in correlation between conditions.
#The first way that the input data can be filtered is by removing genes with a low average -- or more precisely, a low measure of central tendency. 
#The second way that the input data can be filtered is by removing genes that have a low value in a dispersion measure.
#Here I do both. if sequential = FALSE, the same genes removed by the central tendency filtering could potentially be removed by the dispersion filtering.
Alv_mac_normalized_counts_MF_MB_CN_2hpi_filtered = filterGenes(Alv_mac_normalized_counts_MF_MB_CN_2hpi, sequential = TRUE, 
                                                               filterTypes = c("central","dispersion"), filterCentralType = "median",  filterCentralPercentile = 0.3, filterDispersionType = "cv", 
                                                               filterDispersionPercentile = 0.3)
Alv_mac_normalized_counts_MF_MB_CN_24hpi_filtered = filterGenes(Alv_mac_normalized_counts_MF_MB_CN_24hpi, sequential = TRUE, 
                                                               filterTypes = c("central","dispersion"), filterCentralType = "median",  filterCentralPercentile = 0.3, filterDispersionType = "cv", 
                                                               filterDispersionPercentile = 0.3)
Alv_mac_normalized_counts_MF_MB_CN_6hpi_filtered = filterGenes(Alv_mac_normalized_counts_MF_MB_CN_6hpi, sequential = TRUE, 
                                                               filterTypes = c("central","dispersion"), filterCentralType = "median",  filterCentralPercentile = 0.3, filterDispersionType = "cv", 
                                                               filterDispersionPercentile = 0.3)
Alv_mac_normalized_counts_MF_MB_CN_48hpi_filtered = filterGenes(Alv_mac_normalized_counts_MF_MB_CN_48hpi, sequential = TRUE, 
                                                               filterTypes = c("central","dispersion"), filterCentralType = "median",  filterCentralPercentile = 0.3, filterDispersionType = "cv", 
                                                               filterDispersionPercentile = 0.3)

Alv_mac_normalized_counts_MF_TB_CN_2hpi_filtered = filterGenes(Alv_mac_normalized_counts_MF_TB_CN_2hpi, sequential = TRUE, 
                                                               filterTypes = c("central","dispersion"), filterCentralType = "median",  filterCentralPercentile = 0.3, filterDispersionType = "cv", 
                                                               filterDispersionPercentile = 0.3)
Alv_mac_normalized_counts_MF_TB_CN_24hpi_filtered = filterGenes(Alv_mac_normalized_counts_MF_TB_CN_24hpi, sequential = TRUE, 
                                                                filterTypes = c("central","dispersion"), filterCentralType = "median",  filterCentralPercentile = 0.3, filterDispersionType = "cv", 
                                                                filterDispersionPercentile = 0.3)
Alv_mac_normalized_counts_MF_TB_CN_6hpi_filtered = filterGenes(Alv_mac_normalized_counts_MF_TB_CN_6hpi, sequential = TRUE, 
                                                               filterTypes = c("central","dispersion"), filterCentralType = "median",  filterCentralPercentile = 0.3, filterDispersionType = "cv", 
                                                               filterDispersionPercentile = 0.3)
Alv_mac_normalized_counts_MF_TB_CN_48hpi_filtered = filterGenes(Alv_mac_normalized_counts_MF_TB_CN_48hpi, sequential = TRUE, 
                                                                filterTypes = c("central","dispersion"), filterCentralType = "median",  filterCentralPercentile = 0.3, filterDispersionType = "cv", 
                                                                filterDispersionPercentile = 0.3)

#you can then log transform the data, which is often a good idea for RNA-seq data
Alv_mac_normalized_counts_MF_MB_CN_2hpi_log = log(Alv_mac_normalized_counts_MF_MB_CN_2hpi_filtered + 1)
Alv_mac_normalized_counts_MF_MB_CN_24hpi_log = log(Alv_mac_normalized_counts_MF_MB_CN_24hpi_filtered + 1)
Alv_mac_normalized_counts_MF_MB_CN_6hpi_log = log(Alv_mac_normalized_counts_MF_MB_CN_6hpi_filtered + 1)
Alv_mac_normalized_counts_MF_MB_CN_48hpi_log = log(Alv_mac_normalized_counts_MF_MB_CN_48hpi_filtered + 1)
Alv_mac_normalized_counts_MF_TB_CN_2hpi_log = log(Alv_mac_normalized_counts_MF_TB_CN_2hpi_filtered + 1)
Alv_mac_normalized_counts_MF_TB_CN_24hpi_log = log(Alv_mac_normalized_counts_MF_TB_CN_24hpi_filtered + 1)
Alv_mac_normalized_counts_MF_TB_CN_6hpi_log = log(Alv_mac_normalized_counts_MF_TB_CN_6hpi_filtered + 1)
Alv_mac_normalized_counts_MF_TB_CN_48hpi_log = log(Alv_mac_normalized_counts_MF_TB_CN_48hpi_filtered + 1)

#To run the full differential correlation analysis and extract all of the top differentially correlated pairs, run this:
#NOTE: for the purpose of saving time, I have not run permutations, which if adjust = "perm", will run the 
#correlation analysis once, shuffle the samples and run it again as many times as nPerm = whatever. This takes alot of time. 
Alv_mac_normalized_counts_MF_MB_CN_2hpi_ddcor_res_perm_5 = ddcorAll(inputMat = Alv_mac_normalized_counts_MF_MB_CN_2hpi_log, design = design_mat,
                                                                    compare = c("control", "infected"),
                                                                    adjust = "perm", heatmapPlot = FALSE, nPerm = nPerm, nPairs = 3500, corrType = "pearson")

Alv_mac_normalized_counts_MF_MB_CN_24hpi_ddcor_res_perm_5 = ddcorAll(inputMat = Alv_mac_normalized_counts_MF_MB_CN_24hpi_log, design = design_mat,
                                                                    compare = c("control", "infected"),
                                                                    adjust = "perm", heatmapPlot = FALSE, nPerm = nPerm, nPairs = 3500, corrType = "pearson")

Alv_mac_normalized_counts_MF_MB_CN_6hpi_ddcor_res_perm_5 = ddcorAll(inputMat = Alv_mac_normalized_counts_MF_MB_CN_6hpi_log, design = design_mat,
                                                                    compare = c("control", "infected"),
                                                                    adjust = "perm", heatmapPlot = FALSE, nPerm = nPerm, nPairs = 3500, corrType = "pearson")

Alv_mac_normalized_counts_MF_MB_CN_48hpi_ddcor_res_perm_5 = ddcorAll(inputMat = Alv_mac_normalized_counts_MF_MB_CN_48hpi_log, design = design_mat_48hpi,
                                                                    compare = c("control", "infected"),
                                                                    adjust = "perm", heatmapPlot = FALSE, nPerm = nPerm, nPairs = 3500, corrType = "pearson")
#Spearman test																	
Alv_mac_normalized_counts_MF_MB_CN_48hpi_ddcor_res_perm_15_spearman = ddcorAll(inputMat = Alv_mac_normalized_counts_MF_MB_CN_48hpi_log, design = design_mat_48hpi,
                                                                   compare = c("control", "infected"),
                                                                   adjust = "perm", heatmapPlot = FALSE, nPerm = 15, nPairs = 3500, corrType = "spearman")

Alv_mac_normalized_counts_MF_TB_CN_2hpi_ddcor_res_perm_5 = ddcorAll(inputMat = Alv_mac_normalized_counts_MF_TB_CN_2hpi_log, design = design_mat,
                                                                    compare = c("control", "infected"),
                                                                    adjust = "perm", heatmapPlot = FALSE, nPerm = nPerm, nPairs = 3500, corrType = "pearson")

Alv_mac_normalized_counts_MF_TB_CN_24hpi_ddcor_res_perm_5 = ddcorAll(inputMat = Alv_mac_normalized_counts_MF_TB_CN_24hpi_log, design = design_mat,
                                                                     compare = c("control", "infected"),
                                                                     adjust = "perm", heatmapPlot = FALSE, nPerm = nPerm, nPairs = 3500, corrType = "pearson")

Alv_mac_normalized_counts_MF_TB_CN_6hpi_ddcor_res_perm_5 = ddcorAll(inputMat = Alv_mac_normalized_counts_MF_TB_CN_6hpi_log, design = design_mat,
                                                                    compare = c("control", "infected"),
                                                                    adjust = "perm", heatmapPlot = FALSE, nPerm = nPerm, nPairs = 3500, corrType = "pearson")

Alv_mac_normalized_counts_MF_TB_CN_48hpi_ddcor_res_perm_5 = ddcorAll(inputMat = Alv_mac_normalized_counts_MF_TB_CN_48hpi_log, design = design_mat_48hpi,
                                                                     compare = c("control", "infected"),
                                                                     adjust = "perm", heatmapPlot = FALSE, nPerm = nPerm, nPairs = 3500, corrType = "pearson")
																	 
																	 
#larger network test 
Alv_mac_normalized_counts_MF_MB_CN_24hpi_ddcor_res_perm_5_10k_pairs = ddcorAll(inputMat = Alv_mac_normalized_counts_MF_MB_CN_24hpi_log, design = design_mat,
                                                                    compare = c("control", "infected"),
                                                                    adjust = "perm", heatmapPlot = FALSE, nPerm = nPerm, nPairs = 10000, corrType = "pearson")
																	
																	
Alv_mac_normalized_counts_MF_TB_CN_24hpi_ddcor_res_perm_5_10k_pairs = ddcorAll(inputMat = Alv_mac_normalized_counts_MF_TB_CN_24hpi_log, design = design_mat,
                                                                     compare = c("control", "infected"),
                                                                     adjust = "perm", heatmapPlot = FALSE, nPerm = nPerm, nPairs = 10000, corrType = "pearson")


Alv_mac_normalized_counts_MF_MB_CN_24hpi_ddcor_res_perm_5_10k_pairs = ddcorAll(inputMat = Alv_mac_normalized_counts_MF_MB_CN_24hpi_log, design = design_mat,
                                                                    compare = c("control", "infected"),
                                                                    adjust = "perm", heatmapPlot = FALSE, nPerm = nPerm,  corrType = "pearson")
																	
																	
Alv_mac_normalized_counts_MF_TB_CN_24hpi_ddcor_res_perm_5_10k_pairs = ddcorAll(inputMat = Alv_mac_normalized_counts_MF_TB_CN_24hpi_log, design = design_mat,
                                                                     compare = c("control", "infected"),
                                                                     adjust = "perm", heatmapPlot = FALSE, nPerm = nPerm,  corrType = "pearson")
																	 
																	 
																 
Alv_mac_normalized_counts_MF_MB_CN_24hpi_ddcor_res_perm_5_10k_pairs = Alv_mac_normalized_counts_MF_MB_CN_24hpi_ddcor_res_perm_5_10k_pairs[Alv_mac_normalized_counts_MF_MB_CN_24hpi_ddcor_res_perm_5_10k_pairs$empPVals < pval_gene_thresh, ]
																	 
Alv_mac_normalized_counts_MF_TB_CN_24hpi_ddcor_res_perm_5_10k_pairs = Alv_mac_normalized_counts_MF_TB_CN_24hpi_ddcor_res_perm_5_10k_pairs[Alv_mac_normalized_counts_MF_TB_CN_24hpi_ddcor_res_perm_5_10k_pairs$empPVals < pval_gene_thresh, ]																	 
																	 
#then write the networks for cytoscape 
write.csv(Alv_mac_normalized_counts_MF_MB_CN_2hpi_ddcor_res_perm_5, file="Alv_mac_MF_MB_CN_2hpi_ddcor_res_perm_5.csv")
write.csv(Alv_mac_normalized_counts_MF_MB_CN_24hpi_ddcor_res_perm_5, file="Alv_mac_MF_MB_CN_24hpi_ddcor_res_perm_5.csv")
write.csv(Alv_mac_normalized_counts_MF_MB_CN_6hpi_ddcor_res_perm_5, file="Alv_mac_MF_MB_CN_6hpi_ddcor_res_perm_5.csv")
write.csv(Alv_mac_normalized_counts_MF_MB_CN_48hpi_ddcor_res_perm_5, file="Alv_mac_MF_MB_CN_48hpi_ddcor_res_perm_5.csv")
write.csv(Alv_mac_normalized_counts_MF_TB_CN_2hpi_ddcor_res_perm_5, file="Alv_mac_MF_TB_CN_2hpi_ddcor_res_perm_5.csv")
write.csv(Alv_mac_normalized_counts_MF_TB_CN_24hpi_ddcor_res_perm_5, file="Alv_mac_MF_TB_CN_24hpi_ddcor_res_perm_5.csv")
write.csv(Alv_mac_normalized_counts_MF_TB_CN_6hpi_ddcor_res_perm_5, file="Alv_mac_MF_TB_CN_6hpi_ddcor_res_perm_5.csv")
write.csv(Alv_mac_normalized_counts_MF_TB_CN_48hpi_ddcor_res_perm_5, file="Alv_mac_MF_TB_CN_48hpi_ddcor_res_perm_5.csv")
