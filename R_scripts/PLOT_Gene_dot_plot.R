################Dot plots for Genes (TNXB)#######################
setwd("/home/workspace/thall/Analysis/Bovine_network_paper/GWAS_integration/R_analysis/New_GWAS/TNXB_mapping")
TNXB_CN_qvals <- read.delim(file = "MB24_mod20_classes.csv__SNP_set_1e+05_bases.txt n=469114_new_qval.txt")
TNXB_SNPS_only <- read.delim(file = "TNXB.csv__SNP_set_1e+05_bases.txt")
TNXB_Sig <- read.csv(file = "TNXB_SIG.csv")
colnames(TNXB_SNPS_only)[1] <- "SNP"
TNXB_SNPS_WITH_QVALS <- merge(TNXB_SNPS_only, TNXB_CN_qvals, by = "SNP", all.x = TRUE)
colnames(TNXB_SNPS_WITH_QVALS)[1] <- "ID"
TNXB_SNPS_WITH_SIGQVALS <- merge(TNXB_SNPS_WITH_QVALS, TNXB_Sig, by = "ID", all.x = TRUE)
GWAS_file = read.csv("/home/workspace/thall/Analysis/Human_Bovine_comparison_paper/GWAS_integration/GWAS_remapping_data/HOFR_ARS.csv")
colnames(TNXB_SNPS_WITH_SIGQVALS)[1] <- "ID"
TNXB_FINAL = merge(TNXB_SNPS_WITH_SIGQVALS, GWAS_file, by = "ID", all.x = TRUE)

sp2<-ggplot(TNXB_Plot_data, aes(x=Position, y=qval)) + 
geom_point() +
scale_y_continuous(trans = "reverse")
sp2
