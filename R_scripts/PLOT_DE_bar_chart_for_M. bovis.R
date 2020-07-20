########## DE bar chart for M. bovis ##########

#Create vectors with desirable variables for graphing
Time.vec<-c(rep("02hr",2),rep("06hr",2),
            rep("24hr",2),rep("48hr",2))
Variable.vec<-c("Up","Down","Up","Down","Up","Down","Up","Down")
Variable.condition.vec<-c("delta_up","delta_down","delta_up","delta_down","delta_up",
                          "delta_down","delta_up","delta_down")

						 
#Going to count how many genes are up and down so,
#set up empty vector with blank entries to store the gene counts
values.vec<-rep("x",8)

#here, we need the LF2C, so we use the LFC datasets
DE_MB_2hr_logFC <- as.vector(MB2_FDR_0.05["log2FoldChange"])
DE_MB_2hr_logFC<-subset(DE_MB_2hr_logFC, abs(log2FoldChange) > 1)
View(DE_MB_2hr_logFC)
values.vec[1] <-sum(DE_MB_2hr_logFC > 1)
values.vec[2] <-sum(DE_MB_2hr_logFC < 1)

DE_MB_6hr_logFC <- as.vector(MB6_FDR_0.05["log2FoldChange"])
DE_MB_6hr_logFC<-subset(DE_MB_6hr_logFC, abs(log2FoldChange) > 1)
View(DE_MB_6hr_logFC)
values.vec[3] <-sum(DE_MB_6hr_logFC > 1)
values.vec[4] <-sum(DE_MB_6hr_logFC < 1)

DE_MB_24hr_logFC <- as.vector(MB24_FDR_0.05["log2FoldChange"])
DE_MB_24hr_logFC<-subset(DE_MB_24hr_logFC, abs(log2FoldChange) > 1)
View(DE_MB_24hr_logFC)
values.vec[5] <-sum(DE_MB_24hr_logFC > 1)
values.vec[6] <-sum(DE_MB_24hr_logFC < 1)

DE_MB_48hr_logFC <- as.vector(MB48_FDR_0.05["log2FoldChange"])
DE_MB_48hr_logFC<-subset(DE_MB_48hr_logFC, abs(log2FoldChange) > 1)
View(DE_MB_48hr_logFC)
values.vec[7] <-sum(DE_MB_48hr_logFC > 1)
values.vec[8] <-sum(DE_MB_48hr_logFC < 1)

bar_data.raw<-data.frame(a=character(),b=character(),c=numeric(), d=character())
bar_data<-rbind(bar_data.raw, data.frame(a=Time.vec, b=Variable.vec, d=as.numeric(values.vec), e=Variable.condition.vec))
colnames(bar_data)<-c("Time","Variable","Value","Variable.condition")

#Make custom labels for legend of graph to include both italicised and plain text
label_5<-expression(paste(italic(""),"Increased expression"))
label_6<-expression(paste(italic(""),"Decreased expression"))

q<-ggplot(bar_data, aes(Time), ylim(-2500:2500)) + 
  geom_bar(data = subset(bar_data, Variable == "Up"), 
           aes(y = Value, fill = Variable.condition), stat = "identity", position = "dodge",colour="black",size=0.4,width=0.7) +
  scale_fill_manual(values=c("#4d94ff", "#ff3333"), 
                    name=" ",
                    breaks=c("delta_up", "delta_down"), #define the
                    #breaks so that you can relabel
                    labels=c(label_5, label_6)) +
  geom_bar(data = subset(bar_data, Variable == "Down"), #colours are bovis up, tb up, bovis down, tb down
           aes(y = -Value, fill = Variable.condition), stat = "identity", position = "dodge",colour="black",size=0.4,width=0.7) + 
  geom_hline(yintercept = 0,colour = "black") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.text.align = 0)  #aligning the legend labels to legend boxes

#Defining custom subsets for adding stars to the graph
x_up<-(data = subset(bar_data, Variable == "Up"))
x_down<-(data = subset(bar_data, Variable == "Down"))
x_delta_2_up<-subset(x_up[1,])
x_delta_2_down<-(x_down[1,])
x_delta_6_up<-(x_up[2,])
x_delta_6_down<-(x_down[2,])
x_delta_24_up<-(x_up[3,])
x_delta_24_down<-(x_down[3,])
x_delta_48_up<-(x_up[4,])
x_delta_48_down<-(x_down[4,])

q +
  geom_text(data = subset(bar_data, Variable == "Up"), 
            aes(Time, Value, group=Variable.condition, label=Value),
            position = position_dodge(width=0.9), vjust = -0.25, size=3.5) +
  geom_text(data = subset(bar_data, Variable == "Down"), 
            aes(Time, -Value, group=Variable.condition, label=Value),
            position = position_dodge(width=0.9), vjust = 1.25, size=3.5) +
  geom_text(data = x_delta_48_up, 
            aes(Time, Value, group=Variable.condition, label="*"),
            position = position_dodge(width=0.9), vjust = -0.9, hjust= 0.4, size=4) +
  geom_text(data = x_delta_48_down, 
            aes(Time, -Value, group=Variable.condition, label="*"),
            position = position_dodge(width=0.9), vjust = 2.5, hjust= 0.4, size=4) +
  coord_cartesian(ylim = c(-1500, 1500)) +
  scale_x_discrete(name="Time post-infection", breaks=c("02hr","06hr","24hr","48hr"),
                   labels=c("2hr","6hr","24hr","48hr")) + #getting rid of the 0 in 02hr and 06hr
  scale_y_continuous("Number of differentially expressed genes") +
  theme(legend.text=element_text(size=9),legend.key.size=unit(0.4,"cm")) + #changing size of legend
  theme(axis.title.x=element_text(size=11)) +
  theme(axis.title.y=element_text(size=11)) +
  theme(aspect.ratio = 1.3) +
  theme(legend.position="bottom", legend.box = "horizontal") + #horizontal legend at bottom of graph
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size=0.2)) +
  annotate("text", x = 1, y = -1450, size=2.5, label = "* FDR < 0.05")
#change background and axis colour of graph
