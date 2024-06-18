library(data.table)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pracma)
library(ggrepel)
library(viridis)
library(data.table)
library(cowplot)
library(janitor)
library(ggpubr)
library(rstatix)

##

dat.star <- read.table("multiqc_star.txt",header = T)
dat.star$Sample <- gsub("_trunc", "", dat.star$Sample)
dat.star <- setDT(dat.star)

dat.star[, inhibitor_type := lapply(strsplit(Sample, "_"), "[[",2)]
dat.star[, inhibitor_type := factor(inhibitor_type, levels = unique(sort(as.numeric(inhibitor_type))))]
dat.star[is.na(inhibitor_type), inhibitor_type := "RRI"]


p1 <- ggplot(dat.star, aes(x=inhibitor_type, y = uniquely_mapped_percent)) + geom_boxplot() + geom_hline(lty=2, yintercept = dat.star[inhibitor_type == "RRI", median(uniquely_mapped_percent)]) + expand_limits(y=c(0,1)) + theme_cowplot() + labs(y = "%Uniquely Mapped")+ theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))+
  theme(axis.title.x=element_blank())
p2 <- ggplot(dat.star, aes(x=inhibitor_type, y = unmapped_tooshort_percent)) + geom_boxplot() + geom_hline(lty=2, yintercept = dat.star[inhibitor_type == "RRI", median(unmapped_tooshort_percent)]) + expand_limits(y=c(0,1)) + theme_cowplot() + labs(y = "%Too Short")+ theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))+
  theme(axis.title.x=element_blank())
p3 <- ggplot(dat.star, aes(x=inhibitor_type, y = unmapped_other_percent)) + geom_boxplot() + geom_hline(lty=2, yintercept = dat.star[inhibitor_type == "RRI", median(unmapped_other_percent)]) + expand_limits(y=c(0,1)) + theme_cowplot() + labs(y = "%Unmapped (Other)")+ theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))+
  theme(axis.title.x=element_blank())

p4 <- ggplot(dat.star, aes(x=inhibitor_type, y = mismatch_rate)) + geom_boxplot() + geom_hline(lty=2, yintercept = dat.star[inhibitor_type == "RRI", median(mismatch_rate)]) + expand_limits(y=c(0,1)) + theme_cowplot() + labs(y = "Rate of Mismatch")+ theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))+
  theme(axis.title.x=element_blank())
p5 <- ggplot(dat.star, aes(x=inhibitor_type, y = insertion_length)) + geom_boxplot() + geom_hline(lty=2, yintercept = dat.star[inhibitor_type == "RRI", median(insertion_length)]) + expand_limits(y=c(0,1)) + labs(y = "Insertion Length")+ theme_cowplot() + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))+
  theme(axis.title.x=element_blank())
p6 <- ggplot(dat.star, aes(x=inhibitor_type, y = deletion_length)) + geom_boxplot() + geom_hline(lty=2, yintercept = dat.star[inhibitor_type == "RRI", median(deletion_length)]) + expand_limits(y=c(0,1)) + theme_cowplot() + labs(y = "Deletion Length")+ theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))+
  theme(axis.title.x=element_blank())

title_fig <- ggdraw() + draw_label("Mapping stats bulk", fontface='bold')         


figure <- plot_grid(align = "v", p1,p2,p3,p4,p5,p6)

plot_grid(title_fig, figure, ncol=1, rel_heights=c(0.1, 1))

#8x5


pwc <- dat.star.keep %>%
  pairwise_t_test(uniquely_mapped_percent ~ inhibitor_type, p.adjust.method = "bonferroni")

write.csv(pwc, "pwc.csv")

          