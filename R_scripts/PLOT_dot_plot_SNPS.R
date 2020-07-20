#################Dotplots for SNPS###################
#takes 2 columns: gene and position. 

Pathway_histo_pi3k_Akt <- read.csv(file = "Pathway_histo_pi3k-Akt.csv")

e <- ggplot(Pathway_histo_pi3k_Akt, aes(x = Gene, y = Distance, fill = Gene))



#genes
e + geom_dotplot(binaxis = "y", stackdir = "center", aes(color = Gene), method = "histodot", dotsize = 0.5, stackratio = 0.9, show.legend = FALSE) + 
#	theme_minimal() +
	theme_clean() +	
	theme(axis.text.y = element_text(face = "italic"), axis.text=element_text(size=15, colour = "#595959"), axis.title=element_text(size=16, face = "bold", , colour = "#595959" )) +
	labs(y = "Distance from gene", x = "Gene", title = "PI3K-AKT") +
	coord_flip() +
	ylim(-100000, 100000)

#GWAS
e + geom_dotplot(binaxis = "y", stackdir = "centerwhole", aes(color = GWAS), fill = "white", method = "histodot", dotsize = 0.15, stackratio = 1, show.legend = FALSE) + 
#	theme_minimal() +
	theme_clean() +	
	theme(axis.text.y = element_text(face = "italic")) +
	labs(y = "Distance from gene", x = "GWAS", title = "SNP Distribution") +
	coord_flip() +
	ylim(-100000, 100000)

#violin genes
e + geom_violin(trim = FALSE) + 
	geom_dotplot(binaxis = "y", stackdir ="center", dotsize = 0.5, show.legend = FALSE) + 
	stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1),  geom = "pointrange", color = "red", shape = 11) +
	labs(y = "Distance from gene", x = "Gene", title = "PI3K-AKT") +
	theme(axis.text.x = element_text(face = "italic")) +
	ylim(-100000, 100000)


#HOFR
e <- ggplot(HOFR_dot, aes(x = GWAS, y = HOFR_dot$`Distance from Gene`))
e + geom_dotplot(binaxis = "y", stackdir = "centerwhole", color = "#B93442", fill = "#B93442", method = "histodot", dotsize = 0.3, stackratio = 1, show.legend = FALSE) + 
#	theme_minimal() +
	theme_clean() +	
	theme(axis.text.y = element_text(face = "italic")) +
	labs(y = "Distance from gene", x = "GWAS", title = "SNP Distribution") +
	coord_flip() +
	ylim(-100000, 100000)
	
#CH
e <- ggplot(CH_dot, aes(x = GWAS, y = CH_dot$`Distance from Gene`))
e + geom_dotplot(binaxis = "y", stackdir = "centerwhole", color = "#3453B9", fill = "#3453B9", method = "histodot", dotsize = 0.3, stackratio = 1, show.legend = FALSE) + 
#	theme_minimal() +
	theme_clean() +	
	theme(axis.text.y = element_text(face = "italic")) +
	labs(y = "Distance from gene", x = "GWAS", title = "SNP Distribution") +
	coord_flip() +
	ylim(-100000, 100000)
	
#LM
e <- ggplot(LM_dot, aes(x = GWAS, y = LM_dot$`Distance from Gene`))
e + geom_dotplot(binaxis = "y", stackdir = "centerwhole", color = "#3BAF43", fill = "#3BAF43", method = "histodot", dotsize = 0.3, stackratio = 1, show.legend = FALSE) + 
#	theme_minimal() +
	theme_clean() +	
	theme(axis.text.y = element_text(face = "italic")) +
	labs(y = "Distance from gene", x = "GWAS", title = "SNP Distribution") +
	coord_flip() +
	ylim(-100000, 100000)
