########## Heatmaps ##########

#Use the vector object from the venn diagrams and the RLD transformation 
MB2_pheatmap <- pheatmap(assay(rld)[ DE_MB_2hr.vector,],cluster_rows =F)
MB6_pheatmap <- pheatmap(assay(rld)[ DE_MB_6hr.vector,],cluster_rows =F)
MB24_pheatmap <- pheatmap(assay(rld)[ DE_MB_24hr.vector,],cluster_rows =F)
MB48_pheatmap <- pheatmap(assay(rld)[ DE_MB_48hr.vector,],cluster_rows =F)
TB2_pheatmap <- pheatmap(assay(rld)[ DE_TB_2hr.vector,],cluster_rows =F)
TB6_pheatmap <- pheatmap(assay(rld)[ DE_TB_6hr.vector,],cluster_rows =F)
TB24_pheatmap <- pheatmap(assay(rld)[ DE_TB_24hr.vector,],cluster_rows =F)
TB48_pheatmap <- pheatmap(assay(rld)[ DE_TB_48hr.vector,],cluster_rows =F)
