library(tidyverse)
library(dplyr)
library(readr)
library(data.table)
library(plyr)
library(rlist)
library(readxl)
library(RColorBrewer)
library(viridis)
library(plotly)
library(ggplot2)
library(cowplot)
library(purrr)
library(gridExtra)
library(pracma)
library(scales)
library(ggpubr)
library(rstatix)
library(janitor)
library(reshape)

#Exported bioanalyzer "result" files in .csv format
##delete rows and columns using excel to get Time vs sample

#if you haven't merged the files yet
files <- list.files(pattern = 'Bioanalyzer*.*csv', full.names = T) %>%
map(read_csv, skip = 1)
files_bind <- list.cbind(files)
write.table(files_bind, file='SS2_bulk_cDNA.tsv', quote=FALSE, sep='\t', col.names = NA)


##data analysis

SS2_bB <-read.table(file = 'SS2_bulk_cDNA_combined.csv', sep = ';', header = TRUE)
SS2_bB <- type_convert(SS2_bB,trim_ws=TRUE,col_types = cols(Pesos=col_integer(),Alturas=col_double()),locale = locale(decimal_mark = ","))


SS2_B_meta <- read.table('data/SS2_RNA_bulk_BarcodeAnnotations.csv')

SS2_bB_melt <- melt(SS2_bB, id.vars = "Time")

SS2_bB_melt_join <- full_join(SS2_bB_melt, SS2_B_meta, by = c("variable" = "sample_name"))
colnames(SS2_bB_melt_join) [4] <- "Sample"

SS2_bB_melt_join <- as.data.table(SS2_bB_melt_join)
SS2_bB_melt_join <-  filter(SS2_bB_melt_join, Sample != "Ladder")
SS2_bB_melt_join[, inhibitor_type := factor(inhibitor_type, levels = unique(sort(as.numeric(inhibitor_type))))]
SS2_bB_melt_join[is.na(inhibitor_type), inhibitor_type := "RRI"]

###Area under curve for primer dimer/degraded RNA (<75)

SS2_subsetformean = subset(SS2_bB_melt_join, Time < 75)
SS2_subsetformean <- subset(SS2_subsetformean, Time > 45)
SS2_subsetformean$value <- pmax(SS2_subsetformean$value,0)

SS2bulk_df_sum <- SS2_subsetformean[, value2 := sum(value), by = c("Sample", "Time")]
SS2bulk_df_mean <- SS2bulk_df_sum[, mean(value2), by = c("Sample", "Time")]

SS2bulk_df_mean=SS2bulk_df_mean[!is.na(SS2bulk_df_mean$V1),]
unique_groups = as.character(unique(SS2bulk_df_mean$Sample))

##all samples
all_samples = as.character(unique(SS2_subsetformean$variable))
AUC_result_all = c()

for(i in 1:length(all_samples))
{
  df_subset_all = SS2bulk_df_sum[SS2bulk_df_sum$variable %in% all_samples[i],]
  AUC_all = trapz(as.numeric(df_subset_all$Time),as.numeric(df_subset_all$value))
  AUC_result_all[i] = AUC_all
  names(AUC_result_all)[i] = all_samples[i]
}

Area_SS2bulk_all <- as.data.table(AUC_result_all)
samples_all <- as.data.table(all_samples)
Area_SS2bulk_all <- cbind(samples_all,Area_SS2bulk_all)
Area_SS2bulk_all <- left_join(Area_SS2bulk_all, SS2bulk_df_sum, by= c("all_samples" = "variable"))

Area_SS2bulk_plot_all <- Area_SS2bulk_all[!duplicated(Area_SS2bulk_all$all_samples), ]

Area_SS2bulk_plot_all <- as.data.table(Area_SS2bulk_plot_all)
Area_SS2bulk_plot_all <- as.data.frame(Area_SS2bulk_plot_all)
Area_SS2bulk_plot_all[,2:7] <- Area_SS2bulk_plot_all[,2:7] %>% mutate_if(is.character,as.numeric)

write.table(Area_SS2bulk_plot_all[,c(1, 2)], file='SS2_bulk_primerdimer.tsv', quote=FALSE, sep='\t', row.names=F)


Area_SS2bulk_plot_all_2 <- Area_SS2bulk_plot_all %>%
  group_by(Sample) %>%
  summarise(sd = sd(AUC_result_all),
            n = n(),
            se = sd / sqrt(n))


Area_SS2bulk_plot_all_3 <- Area_SS2bulk_plot_all %>%
  group_by(Sample) %>%
  summarise(mean = mean(AUC_result_all))

Area_SS2bulk_plot_all_final <- left_join(Area_SS2bulk_plot_all_2, Area_SS2bulk_plot_all_3,
                                         by= "Sample")


Area_SS2bulk_plot_all_final_true <- Area_SS2bulk_plot_all_final_true %>% mutate_if(is.character,as.numeric)


### area under curve for PEAK (intact cDNA) - use for graph
SS2_subsetformean_peak = subset(SS2_bB_melt_join, Time < 110)
SS2_subsetformean_peak <- subset(SS2_subsetformean_peak, Time > 75)
SS2_subsetformean_peak$value <- pmax(SS2_subsetformean_peak$value,0)


SS2bulk_df_sum_peak <- SS2_subsetformean_peak[, value2 := sum(value), by = c("Sample", "Time")]
SS2bulk_df_mean_peak <- SS2bulk_df_sum_peak[, mean(value2), by = c("Sample", "Time")]

SS2bulk_df_mean_peak=SS2bulk_df_mean_peak[!is.na(SS2bulk_df_mean_peak$V1),]
unique_groups = as.character(unique(SS2bulk_df_mean_peak$Sample))


all_samples_peak = as.character(unique(SS2_subsetformean_peak$variable))
AUC_result_all_peak = c()

for(i in 1:length(all_samples_peak))
{
  df_subset_all_peak = SS2bulk_df_sum_peak[SS2bulk_df_sum_peak$variable %in% all_samples_peak[i],]
  AUC_all_peak = trapz(df_subset_all_peak$Time,df_subset_all_peak$value)
  AUC_result_all_peak[i] = AUC_all_peak
  names(AUC_result_all_peak)[i] = all_samples_peak[i]
}

Area_SS2bulk_all_peak <- as.data.table(AUC_result_all_peak)
samples_all_peak <- as.data.table(all_samples_peak)
Area_SS2bulk_all_peak <- cbind(samples_all_peak,Area_SS2bulk_all_peak)
Area_SS2bulk_all_peak <- left_join(Area_SS2bulk_all_peak, SS2bulk_df_sum_peak, by= c("all_samples_peak" = "variable"))

Area_SS2bulk_plot_all_peak <- Area_SS2bulk_all_peak[!duplicated(Area_SS2bulk_all_peak$all_samples_peak), ]

Area_SS2bulk_plot_all_peak <- as.data.table(Area_SS2bulk_plot_all_peak)
Area_SS2bulk_plot_all_peak <- as.data.frame(Area_SS2bulk_plot_all_peak)
Area_SS2bulk_plot_all_peak[,2:7] <- Area_SS2bulk_plot_all_peak[,2:7] %>% mutate_if(is.character,as.numeric)

write.table(Area_SS2bulk_plot_all_peak[,c(1, 2)], file='SS2_bulk_rawpeak.tsv', quote=FALSE, sep='\t', row.names=F)

Area_SS2bulk_plot_all_peak_2 <- Area_SS2bulk_plot_all_peak %>%
  group_by(Sample) %>%
  summarise(sd = sd(AUC_result_all_peak),
            n = n(),
            se = sd / sqrt(n))


Area_SS2bulk_plot_all_peak_3 <- Area_SS2bulk_plot_all_peak %>%
  group_by(Sample) %>%
  summarise(mean = mean(AUC_result_all_peak))

Area_SS2bulk_plot_all_final_peak <- left_join(Area_SS2bulk_plot_all_peak_2, Area_SS2bulk_plot_all_peak_3,
                                         by= "Sample")

Area_SS2bulk_plot_all_final_true_peak <- Area_SS2bulk_plot_all_final_true_peak %>% mutate_if(is.character,as.numeric)


ggplot(Area_SS2bulk_plot_all_final_true_peak, aes(x=Sample, color=Sample)) +
  geom_point(aes(y=mean, color=Sample)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(.9)) +
  labs(x="Inhibitor concentration", y="Fluorescence Units (Area)") +
  theme_cowplot() + labs(color = "Inhibitor") + geom_vline(aes(xintercept = as.factor(100)),
                                         linetype = 4, colour = "black")

#6x4.5


###


### area under curve for total cDNA (everything) - use to calculate %primer dimer
SS2_subsetformean_totlib = subset(SS2_bB_melt_join, Time < 110)
SS2_subsetformean_totlib <- subset(SS2_subsetformean_totlib, Time > 45)
SS2_subsetformean_totlib$value <- pmax(SS2_subsetformean_totlib$value,0)

SS2bulk_df_sum_totlib <- SS2_subsetformean_totlib[, value2 := sum(value), by = c("Sample", "Time")]
SS2bulk_df_mean_totlib <- SS2bulk_df_sum_totlib[, mean(value2), by = c("Sample", "Time")]

SS2bulk_df_mean_totlib=SS2bulk_df_mean_totlib[!is.na(SS2bulk_df_mean_totlib$V1),]
unique_groups = as.character(unique(SS2bulk_df_mean_totlib$Sample))


all_samples_totlib = as.character(unique(SS2_subsetformean_totlib$variable))
AUC_result_all_totlib = c()

for(i in 1:length(all_samples_totlib))
{
  df_subset_all_totlib = SS2bulk_df_sum_totlib[SS2bulk_df_sum_totlib$variable %in% all_samples_totlib[i],]
  AUC_all_totlib = trapz(df_subset_all_totlib$Time,df_subset_all_totlib$value)
  AUC_result_all_totlib[i] = AUC_all_totlib
  names(AUC_result_all_totlib)[i] = all_samples_totlib[i]
}

Area_SS2bulk_all_totlib <- as.data.table(AUC_result_all_totlib)
samples_all_totlib <- as.data.table(all_samples_totlib)
Area_SS2bulk_all_totlib <- cbind(samples_all_totlib,Area_SS2bulk_all_totlib)
Area_SS2bulk_all_totlib <- left_join(Area_SS2bulk_all_totlib, SS2bulk_df_sum_totlib, by= c("all_samples_totlib" = "variable"))

Area_SS2bulk_plot_all_totlib <- Area_SS2bulk_all_totlib[!duplicated(Area_SS2bulk_all_totlib$all_samples_totlib), ]

Area_SS2bulk_plot_all_totlib <- as.data.table(Area_SS2bulk_plot_all_totlib)
Area_SS2bulk_plot_all_totlib <- as.data.frame(Area_SS2bulk_plot_all_totlib)
Area_SS2bulk_plot_all_totlib[,2:7] <- Area_SS2bulk_plot_all_totlib[,2:7] %>% mutate_if(is.character,as.numeric)

write.table(Area_SS2bulk_plot_all_totlib[,c(1, 2)], file='SS2_bulk_totalarea.tsv', quote=FALSE, sep='\t', row.names=F)

Area_SS2bulk_plot_all_totlib_2 <- Area_SS2bulk_plot_all_totlib %>%
  group_by(Sample) %>%
  summarise(sd = sd(AUC_result_all_totlib),
            n = n(),
            se = sd / sqrt(n))


Area_SS2bulk_plot_all_totlib_3 <- Area_SS2bulk_plot_all_totlib %>%
  group_by(Sample) %>%
  summarise(mean = mean(AUC_result_all_totlib))

Area_SS2bulk_plot_all_final_totlib <- left_join(Area_SS2bulk_plot_all_totlib_2, Area_SS2bulk_plot_all_totlib_3,
                                              by= "Sample")


ggplot(Area_SS2bulk_plot_all_final_true_totlib, aes(x=Sample, color=Sample)) +
  geom_point(aes(y=mean, color=Sample)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(.9)) +
  labs(x="Inhibitor concentration", y="Fluorescence Units (Area)") +
  theme_cowplot() + labs(color = "Inhibitor") + geom_vline(aes(xintercept = as.factor(100)),
                                         linetype = 4, colour = "black")



###put it all together
## primer-dimer
Area_SS2bulk_plot_all_final_true
Area_SS2bulk_plot_all_final_true_4merge <- Area_SS2bulk_plot_all_final_true %>%
  dplyr::rename(
    PD_se = se,
    PD_n = n,
    PD_sd = sd,
    PD_mean = mean
  )

## intact cDNA
Area_SS2bulk_plot_all_final_true_peak
Area_SS2bulk_plot_all_final_true_peak_4merge <- Area_SS2bulk_plot_all_final_true_peak %>%
  dplyr::rename(
    peak_se = se,
    peak_n = n,
    peak_sd = sd,
    peak_mean = mean
  )

##total library
Area_SS2bulk_plot_all_final_true_totlib
Area_SS2bulk_plot_all_final_true_totlib_4merge <- Area_SS2bulk_plot_all_final_true_totlib %>%
  dplyr::rename(
    tot_se = se,
    tot_n = n,
    tot_sd = sd,
    tot_mean = mean
  )


####put it all together

Area_SS2bulk_plot_all #primerdimer
Area_SS2bulk_plot_all_peak #intact cDNA
Area_SS2bulk_plot_all_totlib #all


SS2bulk_summarize_4plot <- cbind(Area_SS2bulk_plot_all,
                                  Area_SS2bulk_plot_all_peak[,2, drop=FALSE],
                                  Area_SS2bulk_plot_all_totlib[,2, drop=FALSE])

SS2bulk_summarize_4plot <- full_join(SS2bulk_summarize_4plot,
                                      Area_SS2bulk_plot_all_final_true_totlib_4merge, by = "Sample")


SS2bulk_summarize_4plot <-  transform(SS2bulk_summarize_4plot, percentpeak = AUC_result_all / AUC_result_all_totlib)


SS2bulk_summarize_4plot_2 <- SS2bulk_summarize_4plot %>%
  group_by(Sample) %>%
  summarise(mean = mean(percentpeak))

SS2bulk_summarize_4plot_3 <- SS2bulk_summarize_4plot %>%
  group_by(Sample) %>%
  summarise(sd = sd(percentpeak),
            n = n(),
            se = sd / sqrt(n))

SS2bulk_summarize_4plot_FINAL <- full_join(SS2bulk_summarize_4plot_2,
                                           SS2bulk_summarize_4plot_3, by = "Sample")


##Actual graph for %primer dimer
ggplot(SS2bulk_summarize_4plot_FINAL, aes(x=Sample, color=Sample)) +
  geom_point(aes(y=mean, color=Sample)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(.9)) +
  labs(x="Inhibitor concentration", y="% Non-cDNA product") +
  theme_cowplot() + scale_color_manual(values=scale_color_area_true) +
  labs(color = "Inhibitor") + geom_vline(aes(xintercept = as.factor(100)),
                                         linetype = 4, colour = "black")

