#################Line plot by group###############
#A more stable, streamlined and updated version is in my Human_bovine_comparison repository
setwd("C:/Users/Thomas Hall/Dropbox/New Alv_mac analysis/Utilities/For paper/Summaries")

library(ggplot2)
library(ggthemes)

#HOFR
rm(Summary)
Summary <- read.csv("HOFR_summary.csv")
Summary <- Summary[order(Summary$SNPS_total),] 
Summary$Distance <- as.character(Summary$Distance)
Summary$probability_0.1 <- -log10(Summary$probability_0.1)
#Summary$probability_0.1 <- abs(Summary$probability_0.1)
#Summary$probability_0.05 <- Summary$probability_0.05 + 0.001 
#Summary$probability_0.1 <- Summary$probability_0.1 * 0.001
#Summary$probability_0.05 <- log10(Summary$probability_0.05)
#Summary$probability_0.05  <- abs(Summary$probability_0.05)

Line <- ggplot(data = Summary,  aes(x= reorder(Distance, +SNPS_total), y=probability_0.1, group = Group, color = Group)) + 
  geom_line(size = 1.5, alpha = 0.9, show.legend = TRUE) +
  geom_point(size = 3.5, show.legend = TRUE) +
  ylim(NA, 4) +
  #	theme_clean() + scale_color_colorblind() +
  theme(axis.text=element_text(size=15, colour = "#595959"), axis.title=element_text(size=16, face = "bold", , colour = "#595959" )) +
  labs(y = "Probability(-log10)", x = "Distance from gene", title = "Holstein Friesian GWAS") 

plot(Line) + theme_hc()+ scale_colour_hc()

#save as PDF 9 x 15 inches


#CH
rm(Summary)
Summary <- read.csv("CH_summary.csv")
Summary <- Summary[order(Summary$SNPS_total),] 
Summary$Distance <- as.character(Summary$Distance)
Summary$probability_0.1 <- -log10(Summary$probability_0.1)
#Summary$probability_0.1 <- abs(Summary$probability_0.1)
#Summary$probability_0.05 <- Summary$probability_0.05 + 0.001 
#Summary$probability_0.1 <- Summary$probability_0.1 * 0.001
#Summary$probability_0.05 <- log10(Summary$probability_0.05)
#Summary$probability_0.05  <- abs(Summary$probability_0.05)

Line <- ggplot(data = Summary,  aes(x= reorder(Distance, +SNPS_total), y=probability_0.1, group = Group, color = Group)) + 
  geom_line(size = 1.5, alpha = 0.9, show.legend = TRUE) +
  geom_point(size = 3.5, show.legend = TRUE) +
  ylim(NA, 4) +
  #	theme_clean() + scale_color_colorblind() +
  theme(axis.text=element_text(size=15, colour = "#595959"), axis.title=element_text(size=16, face = "bold", , colour = "#595959" )) +
  labs(y = "Probability(-log10)", x = "Distance from gene", title = "Charolais GWAS") 
plot(Line) + theme_hc()+ scale_colour_hc()


#LM
rm(Summary)
Summary <- read.csv("LM_summary.csv")
Summary <- Summary[order(Summary$SNPS_total),] 
Summary$Distance <- as.character(Summary$Distance)
Summary$probability_0.1 <- -log10(Summary$probability_0.1)
#Summary$probability_0.1 <- abs(Summary$probability_0.1)
#Summary$probability_0.05 <- Summary$probability_0.05 + 0.001 
#Summary$probability_0.1 <- Summary$probability_0.1 * 0.001
#Summary$probability_0.05 <- log10(Summary$probability_0.05)
#Summary$probability_0.05  <- abs(Summary$probability_0.05)

Line <- ggplot(data = Summary,  aes(x= reorder(Distance, +SNPS_total), y=probability_0.1, group = Group, color = Group)) + 
  geom_line(size = 1.5, alpha = 0.9, show.legend = TRUE) +
  geom_point(size = 3.5, show.legend = TRUE) +
  ylim(NA, 4) +
  #	theme_clean() + scale_color_colorblind() +
  theme(axis.text=element_text(size=15, colour = "#595959"), axis.title=element_text(size=16, face = "bold", , colour = "#595959" )) +
  labs(y = "Probability(-log10)", x = "Distance from gene", title = "Limousin GWAS") 
plot(Line) + theme_hc()+ scale_colour_hc()

#save as PDF 9 x 15 inches


#The averages for bar chart

Average_bar <- ggplot() +
			geom_bar(data = Summary, aes(x= reorder(Summary$Distance, +Summary$SNPS_total), y=Summary$SNPS_total), stat="summary", fun.y = "mean", fill = "#B0CF92")+
			scale_y_continuous(position = "right") + 
			labs(y = "", x = "", title = "") +
			theme(text=element_text(size=16, family="Calibri"), axis.text.x = element_text(hjust=1)) + 
			theme_bw()
plot(Average_bar)

#themes
plot(Line) + theme_calc()+ scale_colour_calc()
plot(Line) + theme_hc()+ scale_colour_hc()
plot(Line) + theme_stata() + scale_color_stata() 
plot(Line) + theme_clean() +scale_color_colorblind()
plot(Line) + theme_fivethirtyeight() +scale_color_colorblind()
plot(Average_bar) + theme_calc()+ scale_colour_calc()
plot(Average_bar) + theme_hc()+ scale_colour_hc()
plot(Average_bar) + theme_stata() + scale_color_stata() 
plot(Average_bar) + theme_clean() +scale_color_colorblind()
plot(Average_bar) + theme_fivethirtyeight() +scale_color_colorblind()
