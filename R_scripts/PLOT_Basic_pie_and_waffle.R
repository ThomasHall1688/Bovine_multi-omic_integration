############################Pie charts####################
#This is for the correlation classes for the CN network

MB24_classes <- data.frame(Alv_mac_MF_MB_CN_24hpi_ddcor_res_perm_5$Classes)
MB48_classes <- data.frame(Alv_mac_MF_MB_CN_24hpi_ddcor_res_perm_5$Classes)

#Now remove correlations that had no change ie 0/0 in MB48
MB24_classes[- grep("0/0", MB48_classes$Correlation_classes.MB48),]

colnames(MB24_classes) <-  c("Classes")
colnames(MB48_classes) <-  c("Classes")

#Make a minimal theme (optional)

blank_theme <- theme_minimal()+
  theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  plot.title=element_text(size=14, face="bold")
  )
  
#make the pie charts and conclude they are a bad way to represent data. 

bp_24 <- ggplot(MB24_classes, aes(x="", fill=Classes))+
    geom_bar(width = 1)+
    coord_polar("y")+
    scale_fill_brewer(palette="Blues")+
    blank_theme +
    theme(axis.text.x=element_blank()) 
bp_24

bp_48 <- ggplot(MB48_classes, aes(x=factor(1), fill=Classes))+
    geom_bar(width = 1)+
    coord_polar("y")+
    scale_fill_brewer(palette="Reds")+
    blank_theme +
    theme(axis.text.x=element_blank()) 
bp_48


#################Waffle charts################

devtools::install_github("hrbrmstr/waffle")
library(waffle)

mydata <- c(``=20, `B`=32, `0`=32, `AB`=16)
waffle(MB24_classes, title = "Yummy waffle pie!")
