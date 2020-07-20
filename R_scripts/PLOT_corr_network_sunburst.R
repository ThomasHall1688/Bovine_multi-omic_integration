############### Sunburst plot #################
#PLots the module hierarchy as a sunburst figure. The Log transform colouring also indicated the module adjusted p.value 

#No colour 

sbobj_MB2 = draw_sunburst_wt_fill(module.df = Alv_mac_MF_MB_CN_2hpi_MEGENA_summary.output$module.table, feat.col = NULL,id.col = "module.id",parent.col = "module.parent")
sbobj_MB2

sbobj_MB24 = draw_sunburst_wt_fill(module.df = Alv_mac_MF_MB_CN_24hpi_MEGENA_summary.output$module.table, feat.col = NULL,id.col = "module.id",parent.col = "module.parent")
sbobj_MB24

sbobj_MB6 = draw_sunburst_wt_fill(module.df = Alv_mac_MF_MB_CN_6hpi_MEGENA_summary.output$module.table, feat.col = NULL,id.col = "module.id",parent.col = "module.parent")
sbobj_MB6

sbobj_MB48 = draw_sunburst_wt_fill(module.df = Alv_mac_MF_MB_CN_48hpi_MEGENA_summary.output$module.table, feat.col = NULL,id.col = "module.id",parent.col = "module.parent")
sbobj_MB48

sbobj_TB2 = draw_sunburst_wt_fill(module.df = Alv_mac_MF_TB_CN_2hpi_MEGENA_summary.output$module.table, feat.col = NULL,id.col = "module.id",parent.col = "module.parent")
sbobj_TB2

sbobj_TB24 = draw_sunburst_wt_fill(module.df = Alv_mac_MF_TB_CN_24hpi_MEGENA_summary.output$module.table, feat.col = NULL,id.col = "module.id",parent.col = "module.parent")
sbobj_TB24

sbobj_TB6 = draw_sunburst_wt_fill(module.df = Alv_mac_MF_TB_CN_6hpi_MEGENA_summary.output$module.table, feat.col = NULL,id.col = "module.id",parent.col = "module.parent")
sbobj_TB6

sbobj_TB48 = draw_sunburst_wt_fill(module.df = Alv_mac_MF_TB_CN_48hpi_MEGENA_summary.output$module.table, feat.col = NULL,id.col = "module.id",parent.col = "module.parent")
sbobj_TB48

#get some coloring (with log transform option)

mdf_MB24 <- Alv_mac_MF_MB_CN_24hpi_MEGENA_summary.output$module.table
mdf_MB24$heat.pvalue = runif(nrow(mdf),0,0.1)
sbobj_MB24_colour = draw_sunburst_wt_fill(module.df = mdf_MB24,feat.col = "heat.pvalue",log.transform = TRUE,
                               fill.type = "continuous",
                               fill.scale = scale_fill_gradient2(low = "white",mid = "white",high = "red",
                                                                 midpoint = -log10(0.05),na.value = "white"), 
                               id.col = "module.id",parent.col = "module.parent")
sbobj_MB24_colour							   


mdf_MB48 <- Alv_mac_MF_MB_CN_48hpi_MEGENA_summary.output$module.table
mdf_MB48$heat.pvalue = runif(nrow(mdf),0,0.1)
sbobj_MB48_colour = draw_sunburst_wt_fill(module.df = mdf_MB48,feat.col = "heat.pvalue",log.transform = TRUE,
                               fill.type = "continuous",
                               fill.scale = scale_fill_gradient2(low = "white",mid = "white",high = "red",
                                                                 midpoint = -log10(0.05),na.value = "white"), 
                               id.col = "module.id",parent.col = "module.parent")
sbobj_MB48_colour



# get discrete coloring done

mdf_MB24$category = factor(sample(x = c("A","B"),size = nrow(mdf),replace = TRUE))
sbobj_MB24_discrete = draw_sunburst_wt_fill(module.df = mdf,feat.col = "category",
								fill.type = "discrete",
								fill.scale = scale_fill_manual(values = c("A" = "red","B" = "blue")), 
								id.col = "module.id",parent.col = "module.parent")
sbobj_MB24_discrete
