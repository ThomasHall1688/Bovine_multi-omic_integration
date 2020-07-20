######################################################################################
################################ GWAS remapping ######################################
######################################################################################
#The "All_sires" GWAS SNPs are aligned to UMD 3.1.1 and contain no ID's, only pvalues, chromsome locations and substitution effect.
#In order to use this with the Alv_mac data we need to update these SNPs to the new genome. 
#This involves cross referencing the GWAS SNP locations to the old genome using the complete 1000bulls genome run 6, which contains 42920228 snps, which gives us the SNP IDs
#Then, use those SNP IDs with ARS1.2PlusY_BQSR.vcf, which contains 110270188 SNPs aligned to the ARS_UCD genome. This will give us the GWAS SNPS with updated locations. 
#Read in data. All three datasets differ in format and arrangement so we must clean 
#load data into Rstudio

setwd("/home/workspace/thall/Analysis/Networks/GWAS_integration/R_analysis/New_GWAS/GWAS_remapping_data")

#Here we can specify what GWAS we want to remap.  v
GWAS_SNPs <- read.delim("All_sires.txt", header = F, sep = ",", fill = TRUE)
GWAS_SNPs <- read.delim("HOFR_sires.txt", header = F, sep = ",", fill = TRUE)
GWAS_SNPs <- read.delim("CH_sires.txt", header = F, sep = ",", fill = TRUE)
GWAS_SNPs <- read.delim("LM_sires.txt", header = F, sep = ",", fill = TRUE)



GWAS_SNPs <- read.delim("HOFR_sires.txt", header = F, sep = ",", fill = TRUE)
ARS_SNPs <- read.delim("ARS1.2PlusY_BQSR.vcf", header = F, sep = "\t", fill = TRUE)
UMD_SNPs <- read.delim("1000bulls_v6_annotated_snps.tab", header = F, sep = "\t", fill = TRUE)
UMD_INDELS <- read.delim("1000bulls_v6_annotated_indels.tab", header = F, sep = "\t", fill = TRUE)

#First ARS_SNPs. It contains meta data in the first rows, so we remove and set the headers. We change #CHROM to CHROM for merging later. 
ARS_SNPs<- tail(ARS_SNPs, -4)
colnames(ARS_SNPs) <- as.character(unlist(ARS_SNPs[1,]))
ARS_SNPs = ARS_SNPs[-1, ]
colnames(ARS_SNPs)[1] <- "CHROM"
ARS_SNPs <- ARS_SNPs[,-c(6,7,8)]
colnames(ARS_SNPs) <- c("CHROM", "POS", "ID", "REF", "ALT")

#Now we clean GWAS_SNPs. Column 2,3 can be removed, as 2 and 3 are repeats of 1, which is in the format we need to convert ARS and UMD to. Also set headers for clarity
GWAS_SNPs <- GWAS_SNPs[,-c(2,3)]
colnames(GWAS_SNPs) <- c("CHROMPOS", "P_value", "Substitution_effect")

#Now the UMD indels
UMD_INDELS  <- UMD_INDELS[, c(1,2,4,5,29)]
colnames(UMD_INDELS ) <-  c("CHROM", "POS", "REF", "ALT", "ID")
UMD_INDELS = UMD_INDELS[-1, ]
UMD_INDELS <- UMD_INDELS[ , c("CHROM", "POS", "ID", "REF", "ALT")]

#Finally UMD_SNPs. It would seem there are redundant columns for the purposes of the integration. I beleive this data is missing for IP reasons. 
#We reomve these and rename  X.CHROM to CHROM for merging later. 
UMD_SNPs <- UMD_SNPs[, c(1,2,4,5,39)]
colnames(UMD_SNPs) <-  c("CHROM", "POS", "REF", "ALT", "ID")
UMD_SNPs = UMD_SNPs[-1, ]
#Rearrange ID so its in the same place as ARS
UMD_SNPs <- UMD_SNPs[ , c("CHROM", "POS", "ID", "REF", "ALT")]


#For merging, we want to create a common column. For this I have picked a combination of chromsome and position, which should give me unique values for each row. 
#First, we need to add Chr to the CHROM column in UMD and ARS 
UMD_INDELS$CHROM <- sub("^", "Chr", UMD_INDELS$CHROM)
UMD_SNPs$CHROM <- sub("^", "Chr", UMD_SNPs$CHROM)

#Now merge the columns CHROM and POS to create CHROMPOS, for merging the tables, and remove whitspace
UMD_INDELS$CHROMPOS <- paste(UMD_INDELS$CHROM, ":", UMD_INDELS$POS)
UMD_INDELS$CHROMPOS <- gsub("\\s+", "", UMD_INDELS$CHROMPOS)
UMD_INDELS <- UMD_INDELS[ , c("CHROMPOS", "CHROM", "POS", "ID", "REF", "ALT")]

UMD_SNPs$CHROMPOS <- paste(UMD_SNPs$CHROM, ":", UMD_SNPs$POS)
UMD_SNPs$CHROMPOS <- gsub("\\s+", "", UMD_SNPs$CHROMPOS)
UMD_SNPs <- UMD_SNPs[ , c("CHROMPOS", "CHROM", "POS", "ID", "REF", "ALT")]

#Now append the indel and snp list 
UMD_SNPs <- rbind(UMD_SNPs, UMD_INDELS)

#ARS has (i think) previous coordinates as RS ID's. I can potentially use this. However its in a wrong format, with a fulll stop instead of 

ARS_SNPs$ID <- chartr(".", ":", ARS_SNPs$ID)

#Now we merge based on the commmon values in CHROMPOS between GWAS_SNPs and UMD_SNPs

GWAS_SNPs_WITH_UMD_RSID <- merge(UMD_SNPs, GWAS_SNPs, by='CHROMPOS')

#in column ID, we have blanks where no RSID exists. However, the CHROMPOS column is the SNP ID in ARS_SNPs. So merging ID and CHROMPOS gives us the full list to merge with ARS_SNPs

GWAS_SNPs_WITH_UMD_RSID <- GWAS_SNPs_WITH_UMD_RSID[ , c("CHROMPOS", "ID", "CHROM", "POS", "REF", "ALT", "P_value", "Substitution_effect")]
#For some totally unknown reason, when you get rid of the . in the blank id spaces, it also gets rid of the r's. 
GWAS_SNPs_WITH_UMD_RSID$ID<- sub(".", "", GWAS_SNPs_WITH_UMD_RSID$ID)
GWAS_SNPs_WITH_UMD_RSID$ID <- ifelse(GWAS_SNPs_WITH_UMD_RSID$ID == "", GWAS_SNPs_WITH_UMD_RSID$CHROMPOS, GWAS_SNPs_WITH_UMD_RSID$ID)
#Now we put back the r's. 
GWAS_SNPs_WITH_UMD_RSID$ID <- sub("s", "rs", GWAS_SNPs_WITH_UMD_RSID$ID)
GWAS_SNPs_WITH_ARS_RSID <- merge(GWAS_SNPs_WITH_UMD_RSID, ARS_SNPs, by='ID')

#Take valid columns
ARS_RSID  <- GWAS_SNPs_WITH_ARS_RSID[, c(1,5,6,7,9,10)]
colnames(ARS_RSID) <- c("ID", "Ref", "Alt", "P_value", "Chromosome", "Position")
ARS_RSID <- ARS_RSID[ , c("ID", "Chromosome", "Position", "Ref", "Alt", "P_value")]

#When the new file is generated, we apply FDR, as we only have unadjusted p-values
library(qvalue)
Q_value <- qvalue(HOFR_ARS$P_value, fdr.level = NULL, pi0 = 1)
HOFR_ARS$q_value  <- Q_value$qvalues

#Here we can now disdinguish what GWAS we want to write 

write.csv(ARS_RSID, file = "All_sires_ARS.csv", row.names = FALSE)
write.csv(ARS_RSID, file = "HOFR_ARS.csv", row.names = FALSE)
write.csv(ARS_RSID, file = "CH_ARS.csv", row.names = FALSE)
write.csv(ARS_RSID, file = "LM_ARS.csv", row.names = FALSE)

#for getting sub effect
LM_sub <- read_csv("LM_sub.csv")
HOFR_sub <- read_csv("HOFR_sub.csv")
CH_sub <- read_csv("CH_sub.csv")

GWAS_SNPs <- read.delim("HOFR_sires.txt", header = F, sep = ",", fill = TRUE)
GWAS_SNPs <- GWAS_SNPs[,-c(2,3)]
colnames(GWAS_SNPs) <- c("CHROMPOS", "P_value", "Substitution_effect")
GWAS_SNPs_WITH_UMD_RSID <- merge(UMD_SNPs, GWAS_SNPs, by='CHROMPOS')
GWAS_SNPs_WITH_UMD_RSID <- GWAS_SNPs_WITH_UMD_RSID[ , c("CHROMPOS", "ID", "CHROM", "POS", "REF", "ALT", "P_value", "Substitution_effect")]
GWAS_SNPs_WITH_UMD_RSID$ID<- sub(".", "", GWAS_SNPs_WITH_UMD_RSID$ID)
GWAS_SNPs_WITH_UMD_RSID$ID <- ifelse(GWAS_SNPs_WITH_UMD_RSID$ID == "", GWAS_SNPs_WITH_UMD_RSID$CHROMPOS, GWAS_SNPs_WITH_UMD_RSID$ID)
GWAS_SNPs_WITH_UMD_RSID$ID <- sub("s", "rs", GWAS_SNPs_WITH_UMD_RSID$ID)
GWAS_SNPs_WITH_ARS_RSID <- merge(GWAS_SNPs_WITH_UMD_RSID, ARS_SNPs, by='ID')
HOFR_sub <- merge(HOFR_sub, GWAS_SNPs_WITH_ARS_RSID, by = "ID", all.x = TRUE)
write.csv(HOFR_sub, file = "HOFR_sub_done.csv", row.names = FALSE)


GWAS_SNPs <- read.delim("CH_sires.txt", header = F, sep = ",", fill = TRUE)
GWAS_SNPs <- GWAS_SNPs[,-c(2,3)]
colnames(GWAS_SNPs) <- c("CHROMPOS", "P_value", "Substitution_effect")
GWAS_SNPs_WITH_UMD_RSID <- merge(UMD_SNPs, GWAS_SNPs, by='CHROMPOS')
GWAS_SNPs_WITH_UMD_RSID <- GWAS_SNPs_WITH_UMD_RSID[ , c("CHROMPOS", "ID", "CHROM", "POS", "REF", "ALT", "P_value", "Substitution_effect")]
GWAS_SNPs_WITH_UMD_RSID$ID<- sub(".", "", GWAS_SNPs_WITH_UMD_RSID$ID)
GWAS_SNPs_WITH_UMD_RSID$ID <- ifelse(GWAS_SNPs_WITH_UMD_RSID$ID == "", GWAS_SNPs_WITH_UMD_RSID$CHROMPOS, GWAS_SNPs_WITH_UMD_RSID$ID)
GWAS_SNPs_WITH_UMD_RSID$ID <- sub("s", "rs", GWAS_SNPs_WITH_UMD_RSID$ID)
GWAS_SNPs_WITH_ARS_RSID <- merge(GWAS_SNPs_WITH_UMD_RSID, ARS_SNPs, by='ID')
CH_sub <- merge(CH_sub, GWAS_SNPs_WITH_ARS_RSID, by = "ID", all.x = TRUE)
write.csv(CH_sub, file = "CH_sub_done.csv", row.names = FALSE)


GWAS_SNPs <- read.delim("LM_sires.txt", header = F, sep = ",", fill = TRUE)
GWAS_SNPs <- GWAS_SNPs[,-c(2,3)]
colnames(GWAS_SNPs) <- c("CHROMPOS", "P_value", "Substitution_effect")
GWAS_SNPs_WITH_UMD_RSID <- merge(UMD_SNPs, GWAS_SNPs, by='CHROMPOS')
GWAS_SNPs_WITH_UMD_RSID <- GWAS_SNPs_WITH_UMD_RSID[ , c("CHROMPOS", "ID", "CHROM", "POS", "REF", "ALT", "P_value", "Substitution_effect")]
GWAS_SNPs_WITH_UMD_RSID$ID<- sub(".", "", GWAS_SNPs_WITH_UMD_RSID$ID)
GWAS_SNPs_WITH_UMD_RSID$ID <- ifelse(GWAS_SNPs_WITH_UMD_RSID$ID == "", GWAS_SNPs_WITH_UMD_RSID$CHROMPOS, GWAS_SNPs_WITH_UMD_RSID$ID)
GWAS_SNPs_WITH_UMD_RSID$ID <- sub("s", "rs", GWAS_SNPs_WITH_UMD_RSID$ID)
GWAS_SNPs_WITH_ARS_RSID <- merge(GWAS_SNPs_WITH_UMD_RSID, ARS_SNPs, by='ID')
LM_sub <- merge(LM_sub, GWAS_SNPs_WITH_ARS_RSID, by = "ID", all.x = TRUE)
write.csv(LM_sub, file = "LM_sub_done.csv", row.names = FALSE)


#This is now our new GWAS dataset, going from 10506081 SNPS to 8808822 (83.84%). Next we obtain SNPs near our genes of interest (All sires)
#This is now our new GWAS dataset, going from 15017692 SNPS to 12740315 (84.83%). Next we obtain SNPs near our genes of interest (HOFR)
#This is now our new GWAS dataset, going from 17250600 SNPS to 14583567 (84.53%). Next we obtain SNPs near our genes of interest (CH)
#This is now our new GWAS dataset, going from 17267260 SNPS to 14586972 (84.47%). Next we obtain SNPs near our genes of interest (LM)
