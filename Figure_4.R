# script to compute delta MiP POC content in each export column


# libraries 
library(ggplot2)
library(readr)
library(tidyverse)
library(dplyr)
library(castr)
library(gridExtra)
library(patchwork)
library("ggpubr")

#import POC data (MiP and MaP)

ecopart <- read_tsv("/Users/alexandreaccardo/Library/CloudStorage/OneDrive-Personnel/Stage M2/GIT/benguela/Ecopart_diagnostics_data_655.tsv")
ecopart <- ecopart %>%
  filter(grepl("a", Profile))


profiles <- unique(ecopart$Profile)
CYCLE_NUMBER <- c(1:184)

# bins (optionnal)
ecopart <- ecopart %>% mutate(binned_depth = round(`Pressure [dbar]`/5)*5 + 2.5)
names(ecopart)[names(ecopart) == "binned_depth"] <- "Depth"
ecopart$Date_Time <- as.Date(ecopart$Date_Time, format = '%Y-%m-%d')

# first export column (MaP)


# during
start_date <- '2021-10-01'
end_date <- '2021-10-17'

column_1_during <- subset(ecopart, Date_Time >= start_date & Date_Time <= end_date) #%>%
  select(Profile, `Pressure [dbar]`, Mip_POC_cont_mgC_m3) %>%
  group_by(`Pressure [dbar]`) %>%
  summarise(mean_mip_POC = mean(Mip_POC_cont_mgC_m3, na.rm = TRUE),
            sd_mip_POC = sd(Mip_POC_cont_mgC_m3, na.rm = TRUE)) %>%
  mutate(mean_mip_POC_smooth = slide(mean_mip_POC, k = 7, fun = mean, n=5, na.rm = TRUE),
         sd_mip_POC_smooth = slide(mean_mip_POC, k=7, fun = sd, n=1, na.rm = TRUE))

# one month before 
start_date <- '2021-09-01'
end_date <- '2021-10-01'

column_1_before <- subset(ecopart, Date_Time >= start_date & Date_Time <= end_date) %>%
  select(Profile, `Pressure [dbar]`, Mip_POC_cont_mgC_m3) %>%
  group_by(`Pressure [dbar]`) %>%
  summarise(mean_mip_POC = mean(Mip_POC_cont_mgC_m3, na.rm = TRUE),
            sd_mip_POC = sd(Mip_POC_cont_mgC_m3, na.rm = TRUE)) %>%
  mutate(mean_mip_POC_smooth = slide(mean_mip_POC, k = 7, fun = mean, n=5, na.rm = TRUE),
         sd_mip_POC_smooth = slide(mean_mip_POC, k = 7, fun = sd, n=1, na.rm = TRUE))


# one month after  
start_date <- '2021-10-17'
end_date <- '2021-11-17'

column_1_after <- subset(ecopart, Date_Time >= start_date & Date_Time <= end_date) %>%
  select(Profile, `Pressure [dbar]`, Mip_POC_cont_mgC_m3) %>%
  group_by(`Pressure [dbar]`) %>%
  summarise(mean_mip_POC = mean(Mip_POC_cont_mgC_m3, na.rm = TRUE),
            sd_mip_POC = sd(Mip_POC_cont_mgC_m3, na.rm = TRUE)) %>%
  mutate(mean_mip_POC_smooth = slide(mean_mip_POC, k = 7, fun = mean, n=5, na.rm = TRUE),
         sd_mip_POC_smooth = slide(mean_mip_POC, k = 7, fun = sd, n=1, na.rm = TRUE))


outside <- merge(column_1_after, column_1_before, by = "Pressure [dbar]") 

outside <- outside %>%
  mutate(mean_mip_POC = rowMeans(outside[, c("mean_mip_POC_smooth.x", "mean_mip_POC_smooth.y")]))


#delta POC 

before_after <- outside %>%
  select(`Pressure [dbar]`, mean_mip_POC)

names(before_after)[names(before_after) == "mean_mip_POC"] <- "mean_mip_POC_before_after"

delta1 <- column_1_during %>% 
  select(`Pressure [dbar]`, mean_mip_POC_smooth) %>%
  merge(before_after, by = "Pressure [dbar]") %>%
  mutate(delta_POC = mean_mip_POC_smooth - mean_mip_POC_before_after,
         ratio_POC = mean_mip_POC_smooth / mean_mip_POC_before_after)

# first map column 
p1 <- ggplot() +
  # before and after mean 
  geom_point(data = column_1_before, aes(x = mean_mip_POC, y = `Pressure [dbar]`), color = "#ef8a62", alpha = 0.5) +
  geom_point(data = column_1_after, aes(x = mean_mip_POC, y = `Pressure [dbar]`), color = "#ef8a62", alpha = 0.5) +
  geom_path(data = outside, mapping = aes(y = `Pressure [dbar]`, x = mean_mip_POC), color = "#ef8a62", size = 1)+
  # during
  geom_point(data = column_1_during, aes(x = mean_mip_POC, y = `Pressure [dbar]`), color = "#1f78b4", alpha = 0.5) +
  geom_path(data = column_1_during, mapping = aes(y = `Pressure [dbar]`, x = mean_mip_POC_smooth), color = "#1f78b4", size = 1)+
  #geom_path(data = column_1_during, mapping = aes(y = `Pressure [dbar]`, x = mean_map_POC_smooth-sd_map_POC_smooth), color = 'yellow')+
  #geom_path(data = column_1_during, mapping = aes(y = `Pressure [dbar]`, x = mean_map_POC_smooth+sd_map_POC_smooth), color = 'yellow')+
  theme_bw()+
  scale_y_reverse(limits = c(700,0))+
  scale_x_continuous(position = "top", limits = c(0,5))+
  labs(y = 'Depth (m)', x = '')+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=20,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"))
p1
# second export column (MaP)


# during
start_date <- '2021-12-01'
end_date <- '2021-12-19'

column_1_during <- subset(ecopart, Date_Time >= start_date & Date_Time <= end_date) %>%
  select(Profile, `Pressure [dbar]`, Mip_POC_cont_mgC_m3) %>%
  group_by(`Pressure [dbar]`) %>%
  summarise(mean_mip_POC = mean(Mip_POC_cont_mgC_m3, na.rm = TRUE),
            sd_mip_POC = sd(Mip_POC_cont_mgC_m3, na.rm = TRUE)) %>%
  mutate(mean_mip_POC_smooth = slide(mean_mip_POC, k = 7, fun = mean, n=5, na.rm = TRUE),
         sd_mip_POC_smooth = slide(mean_mip_POC, k = 7, fun = sd, n=1, na.rm = TRUE))


# one month before 
start_date <- '2021-11-01'
end_date <- '2021-12-01'

column_1_before <- subset(ecopart, Date_Time >= start_date & Date_Time <= end_date) %>%
  select(Profile, `Pressure [dbar]`, Mip_POC_cont_mgC_m3) %>%
  group_by(`Pressure [dbar]`) %>%
  summarise(mean_mip_POC = mean(Mip_POC_cont_mgC_m3, na.rm = TRUE),
            sd_mip_POC = sd(Mip_POC_cont_mgC_m3, na.rm = TRUE)) %>%
  mutate(mean_mip_POC_smooth = slide(mean_mip_POC, k = 7, fun = mean, n=5, na.rm = TRUE),
         sd_mip_POC_smooth = slide(mean_mip_POC, k = 7, fun = sd, n=1, na.rm = TRUE))


# one month after  
start_date <- '2021-12-19'
end_date <- '2022-01-19'

column_1_after <- subset(ecopart, Date_Time >= start_date & Date_Time <= end_date) %>%
  select(Profile, `Pressure [dbar]`, Mip_POC_cont_mgC_m3) %>%
  group_by(`Pressure [dbar]`) %>%
  summarise(mean_mip_POC = mean(Mip_POC_cont_mgC_m3, na.rm = TRUE),
            sd_mip_POC = sd(Mip_POC_cont_mgC_m3, na.rm = TRUE)) %>%
  mutate(mean_mip_POC_smooth = slide(mean_mip_POC, k = 7, fun = mean, n=5, na.rm = TRUE),
         sd_mip_POC_smooth = slide(mean_mip_POC, k = 7, fun = sd, n=1, na.rm = TRUE))

outside <- merge(column_1_after, column_1_before, by = "Pressure [dbar]")

outside <- outside %>%
  mutate(mean_mip_POC = rowMeans(outside[, c("mean_mip_POC_smooth.x", "mean_mip_POC_smooth.y")]))

#delta POC 

before_after <- outside %>%
  select(`Pressure [dbar]`, mean_mip_POC)

names(before_after)[names(before_after) == "mean_mip_POC"] <- "mean_mip_POC_before_after"

delta2 <- column_1_during %>% 
  select(`Pressure [dbar]`, mean_mip_POC_smooth) %>%
  merge(before_after, by = "Pressure [dbar]") %>%
  mutate(delta_POC = mean_mip_POC_smooth - mean_mip_POC_before_after,
         ratio_POC = mean_mip_POC_smooth / mean_mip_POC_before_after)




p2 <- ggplot() +
  # before and after mean 
  geom_point(data = column_1_before, aes(x = mean_mip_POC, y = `Pressure [dbar]`), color = "#ef8a62", alpha = 0.5) +
  geom_point(data = column_1_after, aes(x = mean_mip_POC, y = `Pressure [dbar]`), color = "#ef8a62", alpha = 0.5) +
  geom_path(data = outside, mapping = aes(y = `Pressure [dbar]`, x = mean_mip_POC, color = "Outside"), size = 1)+
  # during
  geom_point(data = column_1_during, aes(x = mean_mip_POC, y = `Pressure [dbar]`), color = "#1f78b4", alpha = 0.5) +
  geom_path(data = column_1_during, mapping = aes(y = `Pressure [dbar]`, x = mean_mip_POC_smooth, color = "Inside"), size = 1)+
  #geom_path(data = column_1_during, mapping = aes(y = `Pressure [dbar]`, x = mean_map_POC_smooth-sd_map_POC_smooth), color = 'yellow')+
  #geom_path(data = column_1_during, mapping = aes(y = `Pressure [dbar]`, x = mean_map_POC_smooth+sd_map_POC_smooth), color = 'yellow')+
  scale_color_manual(values = c(Outside = "#ef8a62", Inside = "#1f78b4"), name = '')+
  theme_bw()+
  scale_y_reverse(limits = c(700,0))+
  scale_x_continuous(position = "top", limits = c(0,5))+
  labs(y = '', x = 'MiP POC content (mgC/m3)')+
  theme(legend.position = "bottom",
        legend.text = element_text(colour = "black", size = 20),
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=20,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black")) 
p2
# third export column (MaP)


# during
start_date <- '2022-03-01'
end_date <- '2022-03-28'

column_1_during <- subset(ecopart, Date_Time >= start_date & Date_Time <= end_date) %>%
  select(Profile, `Pressure [dbar]`, Mip_POC_cont_mgC_m3) %>%
  group_by(`Pressure [dbar]`) %>%
  summarise(mean_mip_POC = mean(Mip_POC_cont_mgC_m3, na.rm = TRUE),
            sd_mip_POC = sd(Mip_POC_cont_mgC_m3, na.rm = TRUE)) %>%
  mutate(mean_mip_POC_smooth = slide(mean_mip_POC, k = 7, fun = mean, n=5, na.rm = TRUE),
         sd_mip_POC_smooth = slide(mean_mip_POC, k = 7, fun = sd, n=1, na.rm = TRUE))

# one month before 
start_date <- '2022-02-01'
end_date <- '2022-03-01'

column_1_before <- subset(ecopart, Date_Time >= start_date & Date_Time <= end_date) %>%
  select(Profile, `Pressure [dbar]`, Mip_POC_cont_mgC_m3) %>%
  group_by(`Pressure [dbar]`) %>%
  summarise(mean_mip_POC = mean(Mip_POC_cont_mgC_m3, na.rm = TRUE),
            sd_mip_POC = sd(Mip_POC_cont_mgC_m3, na.rm = TRUE)) %>%
  mutate(mean_mip_POC_smooth = slide(mean_mip_POC, k = 7, fun = mean, n=5, na.rm = TRUE),
         sd_mip_POC_smooth = slide(mean_mip_POC, k = 7, fun = sd, n=1, na.rm = TRUE))


# one month after  
start_date <- '2022-03-28'
end_date <- '2022-04-28'

column_1_after <- subset(ecopart, Date_Time >= start_date & Date_Time <= end_date) %>%
  select(Profile, `Pressure [dbar]`, Mip_POC_cont_mgC_m3) %>%
  group_by(`Pressure [dbar]`) %>%
  summarise(mean_mip_POC = mean(Mip_POC_cont_mgC_m3, na.rm = TRUE),
            sd_mip_POC = sd(Mip_POC_cont_mgC_m3, na.rm = TRUE)) %>%
  mutate(mean_mip_POC_smooth = slide(mean_mip_POC, k = 7, fun = mean, n=5, na.rm = TRUE),
         sd_mip_POC_smooth = slide(mean_mip_POC, k = 7, fun = sd, n=1, na.rm = TRUE))
outside <- merge(column_1_after, column_1_before, by = "Pressure [dbar]")

outside <- outside %>%
  mutate(mean_mip_POC = rowMeans(outside[, c("mean_mip_POC_smooth.x", "mean_mip_POC_smooth.y")]))

#delta POC 

before_after <- outside %>%
  select(`Pressure [dbar]`, mean_mip_POC)

names(before_after)[names(before_after) == "mean_mip_POC"] <- "mean_mip_POC_before_after"

delta3 <- column_1_during %>% 
  select(`Pressure [dbar]`, mean_mip_POC_smooth) %>%
  merge(before_after, by = "Pressure [dbar]") %>%
  mutate(delta_POC = mean_mip_POC_smooth - mean_mip_POC_before_after, 
         ratio_POC = mean_mip_POC_smooth / mean_mip_POC_before_after)

p3 <- ggplot() +
  # before and after mean 
  geom_point(data = column_1_before, aes(x = mean_mip_POC, y = `Pressure [dbar]`), color = "#ef8a62", alpha = 0.5) +
  geom_point(data = column_1_after, aes(x = mean_mip_POC, y = `Pressure [dbar]`), color = "#ef8a62", alpha = 0.5) +
  geom_path(data = outside, mapping = aes(y = `Pressure [dbar]`, x = mean_mip_POC), color = "#ef8a62", size = 1)+
  # during
  geom_point(data = column_1_during, aes(x = mean_mip_POC, y = `Pressure [dbar]`), color = "#1f78b4", alpha = 0.5) +
  geom_path(data = column_1_during, mapping = aes(y = `Pressure [dbar]`, x = mean_mip_POC_smooth), color = "#1f78b4", size = 1)+
  #geom_path(data = column_1_during, mapping = aes(y = `Pressure [dbar]`, x = mean_map_POC_smooth-sd_map_POC_smooth), color = 'yellow')+
  #geom_path(data = column_1_during, mapping = aes(y = `Pressure [dbar]`, x = mean_map_POC_smooth+sd_map_POC_smooth), color = 'yellow')+
  theme_bw()+
  scale_y_reverse(limits = c(700,0))+
  scale_x_continuous(position = "top", limits = c(0,5))+
  labs(y = '', x = '')+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=20,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"))

p3
List <- list(p1,p2,p3)
#Plot
Plot <- wrap_plots(List,ncol = 3,nrow = 1)
Plot




# delta POC 

ggplot() +
  geom_path(data = delta1, mapping = aes(y = `Pressure [dbar]`, x = delta_POC, color = "MiP export 1"), size = 1)+
  geom_path(data = delta2, mapping = aes(y = `Pressure [dbar]`, x = delta_POC, color = "MiP export 2"), size = 1)+
  geom_path(data = delta3, mapping = aes(y = `Pressure [dbar]`, x = delta_POC, color = "MiP export 3"), size = 1)+
  geom_vline(xintercept = 0, 
             linetype = "dashed", color = "black", size = 1.5)+
  scale_color_manual(values = c('MiP export 1' = "#33a02c", 'MiP export 2' = "#fb9a99", 'MiP export 3' = "#1f78b4"), name = '')+
  theme_bw()+
  scale_y_reverse(limits = c(700,0))+
  scale_x_continuous(position = "top", limits = c(-0.35,1.5))+
  labs(y = 'Depth (m)', x = 'Delta POC (mgC/m3)')+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=20,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"),
        legend.position = "bottom",
        legend.text = element_text(colour = "black", size = 20))

# ratio POC 


ratio_MiP <- ggplot() +
  geom_path(data = delta1, mapping = aes(y = `Pressure [dbar]`, x = ratio_POC, color = 'Export 1'), size = 1)+
  geom_path(data = delta2, mapping = aes(y = `Pressure [dbar]`, x = ratio_POC, color = 'Export 2'), size = 1)+
  geom_path(data = delta3, mapping = aes(y = `Pressure [dbar]`, x = ratio_POC, color = 'Export 3'), size = 1)+
  geom_vline(xintercept = 1, 
             linetype = "dashed", color = "black", size = 1.5)+
  geom_vline(xintercept = mean(delta1$ratio_POC), 
             linetype = "dashed", color = '#33a02c', size = 1)+
  geom_vline(xintercept = mean(delta2$ratio_POC), 
             linetype = "dashed", color = '#fb9a99', size = 1)+
  geom_vline(xintercept = mean(delta3$ratio_POC), 
             linetype = "dashed", color = '#1f78b4', size = 1)+
  scale_color_manual(values = c('Export 1' = "#33a02c", 'Export 2' = "#fb9a99", 'Export 3' = "#1f78b4"), name = '')+
  theme_bw()+
  scale_y_reverse(limits = c(700,0))+
  scale_x_continuous(position = "top") +#, limits = c(-0.35,0.35))+
  labs(y = '', x = 'Ratio (Oustide / Inside)')+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=20,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"),
        legend.position = "bottom",
        legend.text = element_text(colour = "black", size = 20))
ratio_MiP

wilcox.test(delta1$mean_map_POC_smooth, delta2$mean_map_POC_before_after, alternative = 'greater', paired = FALSE)