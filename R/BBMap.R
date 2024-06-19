library(data.table)
library(dplyr)
library(cowplot)
library(ggplot2)
library(stringr)
library(readr)
library(tidyverse)
library(tidyr)
library(pracma)
library(ggrepel)
library(viridis)
library(cowplot)
library(janitor)
library(ggpubr)


##### BBMap

match <- read_tsv("trunc_match_merge_final")
rownames(match) <- NULL
match_t <- t(match)
match_t <- as.data.frame(match_t)
names(match_t)<-sapply(str_remove_all(colnames(match_t),"V"),"[")

match_t <- setDT(match_t, keep.rownames = TRUE)[]
match_t.melt <- melt(match_t, id.vars = "rn")

match_t.melt[, inhibitor_type := lapply(strsplit(rn, "_"), "[[",2)]
match_t.melt[, inhibitor_type := factor(inhibitor_type, levels = unique(sort(as.numeric(inhibitor_type))))]
match_t.melt[is.na(inhibitor_type), inhibitor_type := "RRI"]



matchp <- ggplot(match_t.melt, aes(x=variable, y=(value))) +
  geom_line(data=match_t.melt[inhibitor_type == "RRI", mean(value), by = "variable"], col="black", 
            inherit.aes = F, aes(x=variable, y=V1, group=1)) +
  stat_summary(fun.y="mean", geom="line", show.legend = F, aes(group=inhibitor_type, col=inhibitor_type)) +
  stat_summary(fun.data="mean_cl_boot", geom="ribbon", alpha=0.25, show.legend = F, aes(group=inhibitor_type, fill=inhibitor_type)) +
  facet_wrap(~inhibitor_type, ncol = 6) +
  labs(x="Base position of read", y="Fraction matching genome") +
  coord_cartesian(ylim=c(0,1)) +
  theme_cowplot() + geom_vline(xintercept = 23, linetype="dashed", 
                                 color = "grey30", size=0.4) + 
  theme(strip.background = element_blank(), aspect.ratio = 1, 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.border = element_rect(colour = "black") )



ins <- read_tsv("trunc_ins_merge_final")
rownames(ins) <- NULL
ins_t <- t(ins)
ins_t <- as.data.frame(ins_t)
names(ins_t)<-sapply(str_remove_all(colnames(ins_t),"V"),"[")

ins_t <- setDT(ins_t, keep.rownames = TRUE)[]
ins_t.melt <- melt(ins_t, id.vars = "rn")

ins_t.melt[, inhibitor_type := lapply(strsplit(rn, "_"), "[[",2)]
ins_t.melt[, inhibitor_type := factor(inhibitor_type, levels = unique(sort(as.numeric(inhibitor_type))))]
ins_t.melt[is.na(inhibitor_type), inhibitor_type := "RRI"]


insp <- ggplot(ins_t.melt, aes(x=variable, y=(value))) +
  geom_line(data=ins_t.melt[inhibitor_type == "RRI", mean(value), by = "variable"], col="black", 
            inherit.aes = FALSE, aes(x=variable, y=V1, group=1)) +
  stat_summary(fun.y="mean", geom="line", show.legend = F, aes(group=inhibitor_type, col=inhibitor_type)) +
  stat_summary(fun.data="mean_cl_boot", geom="ribbon", alpha=0.25, show.legend = F, aes(group=inhibitor_type, fill=inhibitor_type)) +
  facet_wrap(~inhibitor_type, ncol = 6) +
  labs(x="Base position of read", y="Fraction with insertion") +
  coord_cartesian(ylim=c(0,0.5)) +
  theme_cowplot() + geom_vline(xintercept = 23, linetype="dashed", 
                              color = "grey30", size=0.4) + 
  theme(strip.background = element_blank(), aspect.ratio = 1, 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.border = element_rect(colour = "black") )


del <- read_tsv("trunc_del_merge_final")
rownames(del) <- NULL
del_t <- t(del)
del_t <- as.data.frame(del_t)
names(del_t)<-sapply(str_remove_all(colnames(del_t),"V"),"[")

del_t <- setDT(del_t, keep.rownames = TRUE)[]
del_t.melt <- melt(del_t, id.vars = "rn")

del_t.melt[, inhibitor_type := lapply(strsplit(rn, "_"), "[[",2)]
del_t.melt[, inhibitor_type := factor(inhibitor_type, levels = unique(sort(as.numeric(inhibitor_type))))]
del_t.melt[is.na(inhibitor_type), inhibitor_type := "RRI"]



delp <- ggplot(del_t.melt, aes(x=variable, y=(del_t.melt$value))) +
  geom_line(data=del_t.melt[inhibitor_type == "RRI", mean(value), by = "variable"], col="black", 
            inherit.aes = FALSE, aes(x=variable, y=V1, group=1)) +
  stat_summary(fun.y="mean", geom="line", show.legend = F, aes(group=inhibitor_type, col=inhibitor_type)) +
  stat_summary(fun.data="mean_cl_boot", geom="ribbon", alpha=0.25, show.legend = F, aes(group=inhibitor_type, fill=inhibitor_type)) +
  facet_wrap(~inhibitor_type, ncol = 6) +
  labs(x="Base position of read", y="Fraction with deletion") +
  coord_cartesian(ylim=c(0,0.1)) +
  theme_cowplot() + geom_vline(xintercept = 23, linetype="dashed", 
                               color = "grey30", size=0.4) + 
  theme(strip.background = element_blank(), aspect.ratio = 1, 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.border = element_rect(colour = "black") )


sub <- read_tsv("./trunc_sub_merge_final")
rownames(sub) <- NULL
sub_t <- t(sub)
sub_t <- as.data.frame(sub_t)
names(sub_t)<-sapply(str_remove_all(colnames(sub_t),"V"),"[")

sub_t <- setDT(sub_t, keep.rownames = TRUE)[]
sub_t.melt <- melt(sub_t, id.vars = "rn")

sub_t.melt[, inhibitor_type := lapply(strsplit(rn, "_"), "[[",2)]
sub_t.melt[, inhibitor_type := factor(inhibitor_type, levels = unique(sort(as.numeric(inhibitor_type))))]
sub_t.melt[is.na(inhibitor_type), inhibitor_type := "RRI"]


subp <- ggplot(sub_t.melt, aes(x=variable, y=(sub_t.melt$value))) +
  geom_line(data=sub_t.melt[inhibitor_type == "RRI", mean(value), by = "variable"], col="black", 
            inherit.aes = FALSE, aes(x=variable, y=V1, group=1)) +
  stat_summary(fun.y="mean", geom="line", show.legend = F, aes(group=inhibitor_type, col=inhibitor_type)) +
  stat_summary(fun.data="mean_cl_boot", geom="ribbon", alpha=0.25, show.legend = F, aes(group=inhibitor_type, fill=inhibitor_type)) +
  facet_wrap(~inhibitor_type, ncol = 6) +
  labs(x="Base position of read", y="Fraction with substitution") +
  coord_cartesian(ylim=c(0,1)) +
  theme_cowplot() + geom_vline(xintercept = 23, linetype="dashed", 
                               color = "grey30", size=0.4) + 
  theme(strip.background = element_blank(), aspect.ratio = 1, 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.border = element_rect(colour = "black") )


BBMap_p <- plot_grid(align = "v", matchp,insp,delp,subp)

title_fig <- ggdraw() + draw_label("BBMap stats SS2 bulk", fontface='bold')

plot_grid(title_fig, BBMap_p, ncol=1, rel_heights=c(0.1, 1))

#10inch x 6inch



