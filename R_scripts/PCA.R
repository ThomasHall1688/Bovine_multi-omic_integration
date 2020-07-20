#for PCA of all groups, for condition and TT
All_PCA <- plotPCA(rld, intgroup=c("condition"))
All_PCA_con <- plotPCA(rld, intgroup=c("condition", "TT"))

#for the individual timepoints for each condition
MB2_PCA <- plotPCA(MB2_rld, "TT")
MB6_PCA <- plotPCA(MB6_rld, "TT")
MB24_PCA <- plotPCA(MB24_rld, "TT")
MB48_PCA <- plotPCA(MB48_rld, "TT")
TB2_PCA <- plotPCA(TB2_rld, "TT")
TB6_PCA <- plotPCA(TB6_rld, "TT")
TB24_PCA <- plotPCA(TB24_rld, "TT")
TB48_PCA <- plotPCA(TB48_rld, "TT")
