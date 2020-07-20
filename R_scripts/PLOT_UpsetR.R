######################################UpsetR plots#########################################

#UpsetR takes a list of genes, with each column named after the DE list, with 1's (present gene) and 0's. I made a combined list of 
setwd("C:/Users/Thomas Hall/Dropbox/New Alv_mac analysis/Analysis/DE analysis/Data/upsetr")
library(readr)
MB_24hpi_DE_genes <- read_csv("MB_24hpi_DE_genes.csv")
MB_2hpi_DE_genes <- read_csv("MB_2hpi_DE_genes.csv")
MB_48hpi_DE_genes <- read_csv("MB_48hpi_DE_genes.csv")
MB_6hpi_DE_genes <- read_csv("MB_6hpi_DE_genes.csv")
Bovine_gene_list <- read_csv("Bovine_gene_list.csv")

MB_2hpi_DE_genes_FDR_0.05 <- as.data.frame(subset(MB_2hpi_DE_genes, padj < 0.05))
MB_6hpi_DE_genes_FDR_0.05 <- as.data.frame(subset(MB_6hpi_DE_genes, padj < 0.05))
MB_24hpi_DE_genes_FDR_0.05 <- as.data.frame(subset(MB_24hpi_DE_genes, padj < 0.05))
MB_48hpi_DE_genes_FDR_0.05 <- as.data.frame(subset(MB_48hpi_DE_genes, padj < 0.05))

colnames(MB_2hpi_DE_genes_FDR_0.05)[1] <- "ENTREZID"
colnames(MB_6hpi_DE_genes_FDR_0.05)[1] <- "ENTREZID"
colnames(MB_24hpi_DE_genes_FDR_0.05)[1] <- "ENTREZID"
colnames(MB_48hpi_DE_genes_FDR_0.05)[1] <- "ENTREZID"

write.csv(MB_2hpi_DE_genes_FDR_0.05, file = "MB_2hpi_DE_genes_FDR_0.05.csv", row.names = FALSE)
write.csv(MB_6hpi_DE_genes_FDR_0.05, file = "MB_6hpi_DE_genes_FDR_0.05.csv", row.names = FALSE)
write.csv(MB_24hpi_DE_genes_FDR_0.05, file = "MB_24hpi_DE_genes_FDR_0.05.csv", row.names = FALSE)
write.csv(MB_48hpi_DE_genes_FDR_0.05, file = "MB_48hpi_DE_genes_FDR_0.05.csv", row.names = FALSE)

MB_2hpi_DE_genes_FDR_0.05$MB_2hpi <- rep(1,nrow(MB_2hpi_DE_genes_FDR_0.05))
MB_24hpi_DE_genes_FDR_0.05$MB_24hpi <- rep(1,nrow(MB_24hpi_DE_genes_FDR_0.05))
MB_6hpi_DE_genes_FDR_0.05$MB_6hpi <- rep(1,nrow(MB_6hpi_DE_genes_FDR_0.05))
MB_48hpi_DE_genes_FDR_0.05$MB_48hpi <- rep(1,nrow(MB_48hpi_DE_genes_FDR_0.05))

MB_6hpi <- MB_6hpi_DE_genes_FDR_0.05[, c(1, 8)]
MB_2hpi <- MB_2hpi_DE_genes_FDR_0.05[, c(1, 8)]
MB_24hpi <- MB_24hpi_DE_genes_FDR_0.05[, c(1, 10)]
MB_48hpi <- MB_48hpi_DE_genes_FDR_0.05[, c(1, 12)]

UpsetR_figure <- merge(Bovine_gene_list, MB_2hpi, by = "ENTREZID", all.x = TRUE)
UpsetR_figure <- merge(UpsetR_figure, MB_6hpi, by = "ENTREZID", all.x = TRUE)
UpsetR_figure <- merge(UpsetR_figure, MB_24hpi, by = "ENTREZID", all.x = TRUE)
UpsetR_figure <- merge(UpsetR_figure, MB_48hpi, by = "ENTREZID", all.x = TRUE)
UpsetR_figure[is.na(UpsetR_figure)] <- 0

#simple
library(UpSetR)

upset(UpsetR_figure, sets = c("MB_2hpi", "MB_6hpi", "MB_24hpi", "MB_48hpi"), 
      sets.bar.color = c("#ff6666", "#809fff", "#838786", "#838786" ), order.by = "freq", empty.intersections = "on",
      matrix.color = "#3d3d5c", main.bar.color = "#52527a" )
							
