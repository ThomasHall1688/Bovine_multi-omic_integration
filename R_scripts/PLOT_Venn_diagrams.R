########## Venn diagrams ##########

library("VennDiagram")     
library("ggplot2")         
library("gridExtra")
library("grid")
library("lattice")

#make a contrast timepoint with hr of interest 
MB24 <- results(dds, contrast=c("TT","MB24","CN24"))
#turn into a dataframe
MB_24hr <- as.data.frame(MB24)
#remove genes that fail cooks distance threshold (NAs)
MB_24hr <- na.omit(MB_24hr)
#Convert the gene ids in col 1 into rownames 
rownames(MB_24hr)<-MB_24hr[,1]
#I think this just creates a vector with all the row names? Makes sense as this is a gene comparison.
DE_MB_24hr.vector <- c(rownames(MB_24hr))

#for LFC > abs(1)
DE_MB_24hr.vector <- c(rownames(MB24_padj_LFC))



MB_2hr <- as.data.frame(MB2)
MB_2hr <- na.omit(MB_2hr)
rownames(MB_2hr)<-MB_2hr[,1]
DE_MB_2hr.vector <- c(rownames(MB_2hr))

#for LFC > abs(1)
DE_MB_2hr.vector <- c(rownames(MB2_padj_LFC))

MB_6hr <- as.data.frame(MB6)
MB_6hr <- na.omit(MB_6hr)
rownames(MB_6hr)<-MB_6hr[,1]
DE_MB_6hr.vector <- c(rownames(MB_6hr))

#for LFC > abs(1)
DE_MB_6hr.vector <- c(rownames(MB6_padj_LFC))


MB_48hr <- as.data.frame(MB48)
MB_48hr <- na.omit(MB_48hr)
rownames(MB_48hr)<-MB_48hr[,1]
DE_MB_48hr.vector <- c(rownames(MB_48hr))

#for LFC > abs(1)
DE_MB_48hr.vector <- c(rownames(MB48_padj_LFC))


TB_2hr <- as.data.frame(TB2)
TB_2hr <- na.omit(TB_2hr)
rownames(TB_2hr)<-TB_2hr[,1]
DE_TB_2hr.vector <- c(rownames(TB_2hr))

#for LFC > abs(1)
DE_TB_2hr.vector <- c(rownames(TB2_padj_LFC))


TB_6hr <- as.data.frame(TB6)
TB_6hr <- na.omit(TB_6hr)
rownames(TB_6hr)<-TB_6hr[,1]
DE_TB_6hr.vector <- c(rownames(TB_6hr))

#for LFC > abs(1)
DE_TB_6hr.vector <- c(rownames(TB6_padj_LFC))


TB_24hr <- as.data.frame(TB24)
TB_24hr <- na.omit(TB_24hr)
rownames(TB_24hr)<-TB_24hr[,1]
DE_TB_24hr.vector <- c(rownames(TB_24hr))

#for LFC > abs(1)
DE_TB_24hr.vector <- c(rownames(TB24_padj_LFC))


TB_48hr <- as.data.frame(TB48)
TB_48hr <- na.omit(TB_48hr)
rownames(TB_48hr)<-TB_48hr[,1]
DE_TB_48hr.vector <- c(rownames(TB_48hr))

#for LFC > abs(1)
DE_TB_48hr.vector <- c(rownames(TB48_padj_LFC))


venn.plot.bovis <- venn.diagram(list(A = DE_MB_2hr.vector,
                                    B = DE_MB_6hr.vector,
                                    C = DE_MB_24hr.vector,
                                    D = DE_MB_48hr.vector),
                               filename        = NULL,
                               fill            = c("#c6e1ef",
                                                   "#a1d7f4",
                                                   "#2aa7ea",
                                                   "#0c15cc"),
                               euler.d=FALSE,
                               scaled=FALSE,
                               sub=substitute( paste(bolditalic('M. bovis'))),
                               sub.fontfamily = "Arial",
                               sub.cex=1.2,
                               sub.pos=c(0.5,1),
                               lwd=1,
                               alpha           = rep(0.50,4),
                               label.col       = "#003333",
                               cex             = 0.85,
                               fontfamily      = "Arial",
                               category.names  = c("2hr",
                                                   "6hr",
                                                   "24hr",
                                                   "48hr"),
                               cat.pos=c(-5,10,10,0),
                               cat.col         = "black",
                               cat.cex         = 0.95,
                               cat.fontfamily  = "Arial",
                               cat.fontface    = 2,
                               rotation.degree = 360,
                               margin          = 0,
                               height          = 10,
                               width           = 4,
                               units           = 'cm',
                               compression     = 'lzw',
                               resolution      = 1200)

venn.plot.TB <- venn.diagram(list(A = DE_TB_2hr.vector,
                                     B = DE_TB_6hr.vector,
                                     C = DE_TB_24hr.vector,
                                     D = DE_TB_48hr.vector),
                                filename        = NULL,
                                fill            = c("#eabbc4",
                                                    "#ed8e9f",
                                                    "#ef3456",
                                                    "#96031d"),
                                euler.d=FALSE,
                                scaled=FALSE,
                                lwd=1,
                                alpha           = rep(0.50,4),
                                sub=substitute( paste(bolditalic('M. tuberculosis'))),
                                sub.fontfamily = "Arial",
                                sub.cex=1.2,
                                sub.pos=c(0.5,1),
                                label.col       = "#003333",
                                cex             = 0.85,
                                fontfamily      = "Arial",
                                category.names  = c("2hr",
                                                    "6hr",
                                                    "24hr",
                                                    "48hr"),
                                cat.pos=c(-5,10,10,0),
                                cat.col         = "black",
                                cat.cex         = 0.9,
                                cat.fontfamily  = "Arial",
                                cat.fontface    = 2,
                                rotation.degree = 360,
                                margin          = 0,
                                height          = 10,
                                width           = 4,
                                units           = 'cm',
                                compression     = 'lzw',
                                resolution      = 1200)
								
								
venn.plot<-grid.arrange(gTree(children=venn.plot.bovis),
             gTree(children=venn.plot.TB),
             ncol = 2,
             widths = c(1,1),
             heights = c(1,1))
			 
			 
#Venn diagrams @ FDR 0.05		


rownames(MB2_FDR_0.05_genes)<-MB2_FDR_0.05_genes[,1]
DE_MB_2hr.vector <- c(rownames(MB2_FDR_0.05))

rownames(MB6_FDR_0.05_genes)<-MB6_FDR_0.05_genes[,1]
DE_MB_6hr.vector <- c(rownames(MB6_FDR_0.05_genes))

rownames(MB24_FDR_0.05_genes)<-MB24_FDR_0.05_genes[,1]
DE_MB_24hr.vector <- c(rownames(MB24_FDR_0.05_genes))

rownames(MB48_FDR_0.05_genes)<-MB48_FDR_0.05_genes[,1]
DE_MB_48hr.vector <- c(rownames(MB48_FDR_0.05_genes))

rownames(TB2_FDR_0.05_genes)<-TB2_FDR_0.05_genes[,1]
DE_TB_2hr.vector <- c(rownames(TB2_FDR_0.05_genes))

rownames(TB6_FDR_0.05_genes)<-TB6_FDR_0.05_genes[,1]
DE_TB_6hr.vector <- c(rownames(TB6_FDR_0.05_genes))

rownames(TB24_FDR_0.05_genes)<-TB24_FDR_0.05_genes[,1]
DE_TB_24hr.vector <- c(rownames(TB24_FDR_0.05_genes))

rownames(TB48_FDR_0.05_genes)<-TB48_FDR_0.05_genes[,1]
DE_TB_48hr.vector <- c(rownames(TB48_FDR_0.05_genes))		

#both venn
#TT font 
windowsFonts(Times=windowsFont("TT Times New Roman"))

venn.plot.bovis_FDR_0.05 <- venn.diagram(list(A = DE_MB_2hr.vector,
                                    B = DE_MB_6hr.vector,
                                    C = DE_MB_24hr.vector,
                                    D = DE_MB_48hr.vector),
                               filename        = NULL,
                               fill            = c("#c6e1ef",
                                                   "#a1d7f4",
                                                   "#2aa7ea",
                                                   "#0c15cc"),
                               euler.d=FALSE,
                               scaled=FALSE,
                               sub=substitute( paste(bolditalic('M. bovis'))),
                               sub.cex=1.2,
                               sub.pos=c(0.5,1),
                               lwd=1,
                               alpha           = rep(0.50,4),
                               label.col       = "#003333",
                               cex             = 0.85,
                               category.names  = c("2hr",
                                                   "6hr",
                                                   "24hr",
                                                   "48hr"),
                               cat.pos=c(-5,10,10,0),
                               cat.col         = "black",
                               cat.cex         = 0.95,
                               cat.fontface    = 2,
                               rotation.degree = 360,
                               margin          = 0,
                               height          = 10,
                               width           = 4,
                               units           = 'cm',
                               compression     = 'lzw',
                               resolution      = 1200)

venn.plot.TB_FDR_0.05 <- venn.diagram(list(A = DE_TB_2hr.vector,
                                     B = DE_TB_6hr.vector,
                                     C = DE_TB_24hr.vector,
                                     D = DE_TB_48hr.vector),
                                filename        = NULL,
                                fill            = c("#eabbc4",
                                                    "#ed8e9f",
                                                    "#ef3456",
                                                    "#96031d"),
                                euler.d=FALSE,
                                scaled=FALSE,
                                lwd=1,
                                alpha           = rep(0.50,4),
                                sub=substitute( paste(bolditalic('M. tuberculosis'))),
                                sub.cex=1.2,
                                sub.pos=c(0.5,1),
                                label.col       = "#003333",
                                cex             = 0.85,
                                category.names  = c("2hr",
                                                    "6hr",
                                                    "24hr",
                                                    "48hr"),
                                cat.pos=c(-5,10,10,0),
                                cat.col         = "black",
                                cat.cex         = 0.9,
                                cat.fontface    = 2,
                                rotation.degree = 360,
                                margin          = 0,
                                height          = 10,
                                width           = 4,
                                units           = 'cm',
                                compression     = 'lzw',
                                resolution      = 1200)

						
venn.plot_FDR_0.05<-grid.arrange(gTree(children=venn.plot.bovis_FDR_0.05),
             gTree(children=venn.plot.TB_FDR_0.05),
             ncol = 2,
             widths = c(1,1),
             heights = c(1,1))		
			 
#M. bovis only 
#A - D is where I can substitute DE genes with common genes from networks. 
venn.plot.bovis_FDR_0.05 <- venn.diagram(list(A = DE_MB_2hr.vector,
                                    B = DE_MB_6hr.vector,
                                    C = DE_MB_24hr.vector,
                                    D = DE_MB_48hr.vector),
                               filename        = NULL,
fill            = c("#d9ffb3",
                    "#afd1a2",
                    "#95b488",
                    "#239023"),
                               euler.d=FALSE,
                               scaled=FALSE,
                               sub=substitute( paste(bolditalic('M. bovis'))),
#                              sub.fontfamily = "Arial",
                               sub.cex=1.2,
                               sub.pos=c(0.5,1),
                               lwd=1,
                               alpha           = rep(0.50,4),
                               label.col       = "#003333",
                               cex             = 0.85,
#                              fontfamily      = "Arial",
                               category.names  = c("2hr",
                                                   "6hr",
                                                   "24hr",
                                                   "48hr"),
                               cat.pos=c(-5,10,10,0),
                               cat.col         = "black",
                               cat.cex         = 0.95,
#                              cat.fontfamily  = "Arial",
                               cat.fontface    = 2,
                               rotation.degree = 360,
                               margin          = 0,
                               height          = 10,
                               width           = 4,
                               units           = 'cm',
                               compression     = 'lzw',
                               resolution      = 1900)
							   
venn.plot_FDR_0.05<-grid.arrange(gTree(children=venn.plot.bovis_FDR_0.05),
             ncol = 2,
             widths = c(1,1),
             heights = c(1,1))
			 
			 
#other possible colour combinations. Substitute the fill option. 
#greens
fill            = c("#d9ffb3",
                    "#afd1a2",
                    "#95b488",
                    "#239023"),
#purples
fill            = c("#f0c3ec",
                    "#cb9ccc",
                    "#ac8aaa",
                    "#602060"),
#browns
fill            = c("#ffffcc",
                    "#d3c0af",
                    "#cca36e",
                    "#e69900"),
#reds
fill            = c("#ffb3b3",
                    "#ff99bb",
                    "#ff3333",
                    "#660000")
					
