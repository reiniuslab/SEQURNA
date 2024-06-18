library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pracma)
library(ggrepel)
library(viridis)
library(data.table)
library(ggplot2)
library(cowplot)
library(janitor)
library(ggpubr)

qualimap_RNAseq <- read.table("merge_qualimap_RNAseq_trunc_final.csv", sep=",", header = F)
qualimap_RNAseq_t <- t(qualimap_RNAseq)
qualimap_RNAseq_t <- janitor::row_to_names(qualimap_RNAseq_t, row_number =1)
rownames(qualimap_RNAseq_t) <- NULL
qualimap_RNAseq_t <- as.data.table(qualimap_RNAseq_t)
qualimap_RNAseq_t[,2:8] <- lapply(qualimap_RNAseq_t[,2:8], as.numeric)

qualimap_RNAseq_t$exonic <- gsub("\\s*\\([^\\)]+\\)","",as.character(qualimap_RNAseq_t$exonic))
qualimap_RNAseq_t$intronic <- gsub("\\s*\\([^\\)]+\\)","",as.character(qualimap_RNAseq_t$intronic))
qualimap_RNAseq_t$intergenic <- gsub("\\s*\\([^\\)]+\\)","",as.character(qualimap_RNAseq_t$intergenic))

qualimap_RNAseq_t[,11:17] <- lapply(qualimap_RNAseq_t[,11:17], as.numeric)

metadata=read.table("data/SS2_RNA_bulk_BarcodeAnnotations.csv",sep=",",header = T)

pq_RNA_data <- qualimap_RNAseq_t %>% left_join(select(metadata, sample_full, sample_id, inhibitor_type), by = c("parameters" = "sample_name"))
pq_RNA_data <- pq_RNA_data[complete.cases(pq_RNA_data[ , 31]), ]

pq_RNA_data <- as.data.table(pq_RNA_data)
pq_RNA_data[, inhibitor_type := factor(inhibitor_type, levels = unique(sort(as.numeric(inhibitor_type))))]
pq_RNA_data[is.na(inhibitor_type), inhibitor_type := "RRI"]



exonic <- ggplot(pq_RNA_data,aes(x=(inhibitor_type), y=exonic/1000))+geom_boxplot()+ scale_fill_grey() +
  theme_bw()+theme(axis.text = element_text(size=12),axis.title = element_text(size=12)) +
  theme(legend.position = "none") + theme_cowplot() +
  geom_hline(lty=2, yintercept = pq_RNA_data[inhibitor_type == "RRI", median(`exonic`/1000)]) +
  ylab("% reads exonic") + xlab("inhibitor_type Type") + theme(axis.text.x = element_text(angle = 90))+
  theme(axis.title.x=element_blank())

intronic <- ggplot(pq_RNA_data,aes(x=(inhibitor_type), y=intronic/1000))+geom_boxplot()+ scale_fill_grey() +
  theme_bw()+theme(axis.text = element_text(size=12),axis.title = element_text(size=12)) +
  theme(legend.position = "none") + theme_cowplot() +
  geom_hline(lty=2, yintercept = pq_RNA_data[inhibitor_type == "RRI", median(`intronic`/1000)]) +
  ylab("% intronic") + xlab("inhibitor_type Type") + theme(axis.text.x = element_text(angle = 90))+
  theme(axis.title.x=element_blank())

intergenic <- ggplot(pq_RNA_data,aes(x=(inhibitor_type), y=intergenic/1000))+geom_boxplot()+ scale_fill_grey() +
  theme_bw()+theme(axis.text = element_text(size=12),axis.title = element_text(size=12)) +
  theme(legend.position = "none") + theme_cowplot() +
  geom_hline(lty=2, yintercept = pq_RNA_data[inhibitor_type == "RRI", median(`intergenic`/1000)]) +
  ylab("% intergenic") + xlab("inhibitor_type Type") + theme(axis.text.x = element_text(angle = 90))+
  theme(axis.title.x=element_blank())

readsaligned <- ggplot(pq_RNA_data,aes(x=(inhibitor_type), y=readsaligned))+geom_boxplot()+ scale_fill_grey() +
  theme_bw()+theme(axis.text = element_text(size=12),axis.title = element_text(size=12)) +
  theme(legend.position = "none") + theme_cowplot() +
  geom_hline(lty=2, yintercept = pq_RNA_data[inhibitor_type == "RRI", median(`readsaligned`)]) +
  ylab("reads aligned") + xlab("inhibitor_type Type") + theme(axis.text.x = element_text(angle = 90))+
  theme(axis.title.x=element_blank())



figure <- plot_grid(align = "v", exonic,intronic,intergenic,readsaligned)

title_fig <- ggdraw() + draw_label("Qualimap stats bulk MTTF RRI", fontface='bold')

plot_grid(title_fig, figure, ncol=1, rel_heights=c(0.1, 1))

#8x6



pq_RNA_data_stack <- subset(pq_RNA_data, select=c("parameters","exonic","intronic","intergenic",
                                                  "inhibitor_type"))

pq_RNA_data_stack.melt <- melt(pq_RNA_data_stack, id=c("parameters","inhibitor_type"))

#pq_RNA_data_stack.melt_1 <- pq_RNA_data_stack.melt %>% 
#  group_by(inhibitor_type, variable) %>% 
#  summarise(sd = sd(value),
#            n = n(),
#            se = sd / sqrt(n))


#pq_RNA_data_stack.melt_2 <- pq_RNA_data_stack.melt %>% 
#  group_by(inhibitor_type, variable) %>% 
#  summarise(mean = mean(value))


#pq_RNA_data_stack.melt_3 <- left_join(pq_RNA_data_stack.melt_1, pq_RNA_data_stack.melt_2, 
#                                      by = c("inhibitor_type", "variable") )


pq_RNA_data_stack.melt_fraction <- transform(pq_RNA_data_stack.melt, new_value = value / 100000)



stacked_color <- c( "#942c80", "#e85362" , "#fecd90" )

p <- ggbarplot(pq_RNA_data_stack.melt_fraction, 
               x = "inhibitor_type", y = "new_value", add = "mean_se",
               error.plot = "errorbar", fill = "variable", legend.title= "Genomic Alignment" ,
               add.params = list(group = "variable") , ylim=c(0,1)) + 
  theme(axis.text.x = element_text(angle = 90))  +
  labs(x="Sample", y="Reads")

p


