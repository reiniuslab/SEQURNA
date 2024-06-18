library(data.table)
library(dplyr)
library(cowplot)
library(ggplot2)
library(plyr)

##don't forget to fix the plotprofile file in excel before (delete things and add Percentile as colname)
dat <- fread("plotprofile_SS2.tab", head = T)

dat.melt <- melt(dat, id.vars = "bins")
dat.melt <- setDT(dat.melt)
dat.melt[, norm.coverage := (value-min(value))/(max(value)-min(value)), by="bins"]


dat.melt[, inhibitor_type := lapply(strsplit(bins, "_"), "[[",2)]
dat.melt[, inhibitor_type := factor(inhibitor_type, levels = unique(sort(as.numeric(inhibitor_type))))]
dat.melt[is.na(inhibitor_type), inhibitor_type := "RRI"]

count(dat.melt, "inhibitor_type")


ggplot(dat.melt, aes(x=variable, y=norm.coverage)) +
  geom_line(data=dat.melt[inhibitor_type == "RRI", mean(norm.coverage), 
        by = "variable"], col="black", inherit.aes = F, 
        aes(x=variable, y=V1, group=1)) +
  stat_summary(fun.y="mean", geom="line", show.legend = F, 
        aes(group=inhibitor_type, col=inhibitor_type)) +
  stat_summary(fun.data="mean_cl_boot", geom="ribbon", 
        alpha=0.25, show.legend = F, aes(group=inhibitor_type, fill=inhibitor_type)) +
  facet_wrap(~inhibitor_type, ncol = 6) +
  labs(x="Gene body percentile (5'−>3')", y="Coverage") +
  ## coord_cartesian(ylim=c(0,1)) +
  theme_cowplot() +
  theme(strip.background = element_blank(), aspect.ratio = 1, 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        panel.border = element_rect(colour = "black") ) + ggtitle("geneBodyCoverage SS2 bulk")
