########## Plot counts of a gene for each time point ##########
#will plot the gene of interest (2nd line, gene = "",) as a bar chart with error bars
#Fix the positions of the the control and infected (use only when needed)
positions <- c("MB", "CN")

#create matrix of counts and samples based on TT grouping
p <- plotCounts(dds, gene=which.min(MB2_res$padj), intgroup="TT", returnData = TRUE)
#select time point of interest
p <- p %>% filter(TT == "MB2" | TT == "CN2")
#calculate stats for error bars.
p.summ <- p %>% group_by(TT) %>% summarize(Mean = mean(count), Min = min(count), Max = max(count), sd = sd(count), n = n(), se = sd(count)/sqrt(n()))
#create bar chart 
Gene_plot = ggplot(p.summ, aes(TT, Mean)) + 
     geom_bar(stat = "identity", width = 0.5, fill = c("#4d94ff", "#ff3333")) + 
     geom_errorbar(aes(ymin = Mean - sd, ymax = Mean + sd), width=0.2) 
Gene_plot + labs(y="Normalised counts", x = "Condition", title = expression(italic("Gene_name")))

p <- plotCounts(dds, gene=which.min(MB6_res$padj), intgroup="TT", returnData = TRUE)
p <- p %>% filter(TT == "MB6" | TT == "CN6")
p.summ <- p %>% group_by(TT) %>% summarize(Mean = mean(count), Min = min(count), Max = max(count), sd = sd(count), n = n(), se = sd(count)/sqrt(n()))
Gene_plot = ggplot(p.summ, aes(TT, Mean)) + 
     geom_bar(stat = "identity", width = 0.5, fill = c("#4d94ff", "#ff3333")) + 
     geom_errorbar(aes(ymin = Mean - sd, ymax = Mean + sd), width=0.2) 
Gene_plot + labs(y="Normalised counts", x = "Condition", title = expression(italic("Gene_name")))

p <- plotCounts(dds, gene=which.min(MB24_res$padj), intgroup="TT", returnData = TRUE)
p <- p %>% filter(TT == "MB24" | TT == "CN24")
p.summ <- p %>% group_by(TT) %>% summarize(Mean = mean(count), Min = min(count), Max = max(count), sd = sd(count), n = n(), se = sd(count)/sqrt(n()))
Gene_plot = ggplot(p.summ, aes(TT, Mean)) + 
     geom_bar(stat = "identity", width = 0.5, fill = c("#4d94ff", "#ff3333")) + 
     geom_errorbar(aes(ymin = Mean - sd, ymax = Mean + sd), width=0.2) 
Gene_plot + labs(y="Normalised counts", x = "Condition", title = expression(italic("Gene_name")))

p <- plotCounts(dds, gene=which.min(MB48_res$padj), intgroup="TT", returnData = TRUE)
p <- p %>% filter(TT == "MB48" | TT == "CN48")
p.summ <- p %>% group_by(TT) %>% summarize(Mean = mean(count), Min = min(count), Max = max(count), sd = sd(count), n = n(), se = sd(count)/sqrt(n()))
Gene_plot = ggplot(p.summ, aes(TT, Mean)) + 
     geom_bar(stat = "identity", width = 0.5, fill = c("#4d94ff", "#ff3333")) + 
     geom_errorbar(aes(ymin = Mean - sd, ymax = Mean + sd), width=0.2) 
Gene_plot + labs(y="Normalised counts", x = "Condition", title = expression(italic("Gene_name")))

p <- plotCounts(dds, gene=which.min(TB2_res$padj), intgroup="TT", returnData = TRUE)
p <- p %>% filter(TT == "TB2" | TT == "CN2")
p.summ <- p %>% group_by(TT) %>% summarize(Mean = mean(count), Min = min(count), Max = max(count), sd = sd(count), n = n(), se = sd(count)/sqrt(n()))
Gene_plot = ggplot(p.summ, aes(TT, Mean)) + 
     geom_bar(stat = "identity", width = 0.5, fill = c("#4d94ff", "#ff3333")) + 
     geom_errorbar(aes(ymin = Mean - sd, ymax = Mean + sd), width=0.2) 
Gene_plot + labs(y="Normalised counts", x = "Condition", title = expression(italic("Gene_name")))

p <- plotCounts(dds, gene=which.min(TB6_res$padj), intgroup="TT", returnData = TRUE)
p <- p %>% filter(TT == "TB6" | TT == "CN6")
p.summ <- p %>% group_by(TT) %>% summarize(Mean = mean(count), Min = min(count), Max = max(count), sd = sd(count), n = n(), se = sd(count)/sqrt(n()))
Gene_plot = ggplot(p.summ, aes(TT, Mean)) + 
     geom_bar(stat = "identity", width = 0.5, fill = c("#4d94ff", "#ff3333")) + 
     geom_errorbar(aes(ymin = Mean - sd, ymax = Mean + sd), width=0.2) 
Gene_plot + labs(y="Normalised counts", x = "Condition", title = expression(italic("Gene_name")))

p <- plotCounts(dds, gene=which.min(TB24_res$padj), intgroup="TT", returnData = TRUE)
p <- p %>% filter(TT == "TB24" | TT == "CN24")
p.summ <- p %>% group_by(TT) %>% summarize(Mean = mean(count), Min = min(count), Max = max(count), sd = sd(count), n = n(), se = sd(count)/sqrt(n()))
Gene_plot = ggplot(p.summ, aes(TT, Mean)) + 
     geom_bar(stat = "identity", width = 0.5, fill = c("#4d94ff", "#ff3333")) + 
     geom_errorbar(aes(ymin = Mean - sd, ymax = Mean + sd), width=0.2) 
Gene_plot + labs(y="Normalised counts", x = "Condition", title = expression(italic("Gene_name")))

p <- plotCounts(dds, gene=which.min(TB48_res$padj), intgroup="TT", returnData = TRUE)
p <- p %>% filter(TT == "TB48" | TT == "CN48")
p.summ <- p %>% group_by(TT) %>% summarize(Mean = mean(count), Min = min(count), Max = max(count), sd = sd(count), n = n(), se = sd(count)/sqrt(n()))
Gene_plot = ggplot(p.summ, aes(TT, Mean)) + 
     geom_bar(stat = "identity", width = 0.5, fill = c("#4d94ff", "#ff3333")) + 
     geom_errorbar(aes(ymin = Mean - sd, ymax = Mean + sd), width=0.2) 
Gene_plot + labs(y="Normalised counts", x = "Condition", title = expression(italic("Gene_name")))
