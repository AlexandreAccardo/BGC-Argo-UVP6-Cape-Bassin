# Compute proportion of morphotypes depending on depth # 

# library 
library(readr)
library(tidyverse)
library(ggplot2)


# import morphotypes concentrations 
# first try with all individuals and after with only selected ones

morpho_conc <- read_csv('C:/Users/lemoi/Desktop/GIT/k-means/clusters_concentrations_all.csv') %>% 
  select(profile, depth, Cluster, n, conc) %>%
  filter(depth < 1000)
morpho_conc$Cluster <-  ifelse(morpho_conc$Cluster == "cluster 1", "Small",
                             ifelse(morpho_conc$Cluster == "cluster 2", "Bright",
                                    ifelse(morpho_conc$Cluster == "cluster 3", "Elongated",
                                           ifelse(morpho_conc$Cluster == "cluster 4", "Aggregates",NA))))


# compute layer depth 
morpho_conc$layer <- 
  ifelse(morpho_conc$depth < 100, "0-100m",
  ifelse(morpho_conc$depth >= 100 & morpho_conc$depth < 300, "100-300m",
  ifelse(morpho_conc$depth >= 300 & morpho_conc$depth < 600, "300-600m",
  ifelse(morpho_conc$depth >= 600 & morpho_conc$depth < 1000, "600-1000m", NA))))

# compute the mean concentration in each layer depth 

morpho_conc_mean <- morpho_conc %>%
  group_by(profile, Cluster, layer) %>%
  summarise(n = sum(n),
            conc = mean(conc, na.rm = TRUE))
# import date, longitude, latitude for each profile 
profile <-  read_csv("C:/Users/lemoi/Desktop/GIT/benguela/Data_created/list_of_profiles.csv")

morpho_conc_mean <- merge(morpho_conc_mean, profile, by = 'profile')

# convert into date format 
morpho_conc_mean$date <- as.Date(as.character(morpho_conc_mean$date), "%Y%m%d")


# select only profiles that are in export columns and compute the mean with all profile that are in the column and represent it 


# import results from STARS method

events_date <- read_csv("C:/Users/lemoi/Desktop/GIT/STARS_method/events_date.csv")

# first export column
start_date <- '2021-10-05'
end_date <- '2021-10-20'

column_1 <- subset(morpho_conc_mean, date >= start_date & date <= end_date)


column_1 <- column_1 %>%
  group_by(layer, Cluster) %>%
  summarise(n = sum(n),
            conc = mean(conc, na.rm = TRUE))%>%
  mutate(total_count = sum(n),
         conc_total = sum(conc),
         proportion = conc*100/conc_total)

# one month before first export column
start_date <- '2021-10-17'
end_date <- '2021-11-01'

column_1_before <- subset(morpho_conc_mean, date >= start_date & date <= end_date)


column_1_before <- column_1_before %>%
  group_by(layer, Cluster) %>%
  summarise(n = sum(n),
            conc = mean(conc, na.rm = TRUE))%>%
  mutate(total_count = sum(n),
         conc_total = sum(conc),
         proportion = conc*100/conc_total)

# contingence table for Khi2 test

table_1 <- column_1 %>%
  select(layer, Cluster, n) %>%
  pivot_wider(names_from = layer, values_from = n) %>%
  select(-Cluster) # cluster are order in the alphabetic order (Aggregates, Bright, Elongated, Small)


chisq.test(table_1, correct = FALSE)
chisq.test(table_1,correct=F)$expected # Cochran criterion verified

# There is a difference between layers but now I need to know where 
# khi2 test 2-2 with adjusted p-values 
layer1vs2 <- chisq.test(table_1[, 1:2], correct = FALSE)
p1 <- layer1vs2$p.value

layer2vs3 <- chisq.test(table_1[, 2:3], correct = FALSE)
p2 <- layer2vs3$p.value

layer3vs4 <- chisq.test(table_1[, 3:4], correct = FALSE)
p3 <- layer3vs4$p.value

layer1vs4 <- chisq.test(table_1[, c(1,4)], correct = FALSE)
p4 <- layer1vs4$p.value

layer2vs4 <- chisq.test(table_1[, c(2,4)], correct = FALSE)
p5 <- layer2vs4$p.value

layer1vs3 <- chisq.test(table_1[, c(1,3)], correct = FALSE)
p6 <- layer1vs3$p.value


p_adjusted <- p.adjust(c(p1, p2, p3),method="holm")
p_adjusted


# create the pie chart using ggplot
ggplot(data = column_1, aes(x = "", y = proportion, fill = Cluster)) +
  geom_bar(width = 1, stat = "identity", color= 'white') +
  coord_polar(theta = "y") +
  scale_fill_manual('Morphotypes',values=c("#CD8BFF", "#00BFC4", "#7AAF00", "#F8766D"))+
  facet_wrap(~layer, ncol = 2) +
  labs(x = NULL, y = NULL, fill = "Cluster") +
  #ggtitle("Proportion of moprhotypes (all individuals) depending on depth in column 1") +
  geom_text(aes(label = paste0(round(proportion), "%")), 
            position = position_stack(vjust = 0.5), color = 'black')+
  theme_void()+
  theme(plot.title.position = "plot",
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"))

plot1 <- ggplot(column_1, aes(fill=Cluster, y=as.factor(layer), x=proportion)) + 
  scale_fill_manual('',values=c("#CD8BFF", "#00BFC4", "#7AAF00", "#F8766D"))+
  geom_bar(position="stack", stat="identity")+
  geom_text(aes(label = paste0(round(column_1$proportion), "%")), 
            position = position_stack(vjust = 0.5), color = 'white', fontface = 'bold')+
  scale_y_discrete(limits=rev)+
  theme_classic()+
  labs(y = 'Layer depth (m)', x = '')+
  theme(legend.position = "none",
    panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=20,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"))

plot1

# second export column
start_date <- '2021-11-25'
end_date <- '2021-12-22'

column_2 <- subset(morpho_conc_mean, date >= start_date & date <= end_date)
column_2 <- column_2 %>%
  group_by(layer, Cluster) %>%
  summarise(n = sum(n),
            conc = mean(conc, na.rm = TRUE))%>%
  mutate(total_count = sum(n),
         conc_total = sum(conc),
         proportion = conc*100/conc_total)

# contingence tabel for Khi2 test

table_2 <- column_2 %>%
  select(layer, Cluster, n) %>%
  pivot_wider(names_from = layer, values_from = n) %>%
  select(-Cluster) # cluster are order in the alphabetic order (Aggregates, Bright, Elongated, Small)


chisq.test(table_2, correct = FALSE)
chisq.test(table_2,correct=F)$expected # Cochran criterion verified

# There is a difference between layers but now I need to know where 
# khi2 test 2-2 with adjusted p-values 
layer1vs2 <- chisq.test(table_2[, 1:2], correct = FALSE)
p1 <- layer1vs2$p.value

layer2vs3 <- chisq.test(table_2[, 2:3], correct = FALSE)
p2 <- layer2vs3$p.value

layer3vs4 <- chisq.test(table_2[, 3:4], correct = FALSE)
p3 <- layer3vs4$p.value

layer1vs4 <- chisq.test(table_2[, c(1,4)], correct = FALSE)
p4 <- layer1vs4$p.value

layer2vs4 <- chisq.test(table_2[, c(2,4)], correct = FALSE)
p5 <- layer2vs4$p.value

layer1vs3 <- chisq.test(table_2[, c(1,3)], correct = FALSE)
p6 <- layer1vs3$p.value


p_adjusted <- p.adjust(c(p1, p2, p3),method="holm")
p_adjusted

# create the pie chart using ggplot
ggplot(data = column_2, aes(x = "", y = proportion, fill = Cluster)) +
  geom_bar(width = 1, stat = "identity", color= 'white') +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Pastel1") +
  facet_wrap(~layer, ncol = 2) +
  labs(x = NULL, y = NULL, fill = "Cluster") +
  #ggtitle("Proportion of moprhotypes (all individuals) depending on depth in column 2") +
  geom_text(aes(label = paste0(round(proportion), "%")), 
            position = position_stack(vjust = 0.5), color = 'black')+
  theme_void()+
  theme(plot.title.position = "plot",
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"))

plot2 <- ggplot(column_2, aes(fill=Cluster, y=as.factor(layer), x=proportion)) + 
  scale_fill_manual('',values=c("#CD8BFF", "#00BFC4", "#7AAF00", "#F8766D"))+
  geom_bar(position="stack", stat="identity")+
  geom_text(aes(label = paste0(round(column_2$proportion), "%")), 
            position = position_stack(vjust = 0.5), color = 'white', fontface = 'bold')+
  scale_y_discrete(limits=rev)+
  theme_classic()+
  labs(y = '', x = 'Proportion (%)')+
  theme(legend.position = "top",
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=20,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"),
        axis.text.y = element_blank(), 
        legend.text = element_text(colour = "black", size = 20))+
  guides(color = guide_legend(nrow = 2))

plot2
# third export column
start_date <- '2022-03-10'
end_date <- '2022-04-15'

column_3 <- subset(morpho_conc_mean, date >= start_date & date <= end_date)
column_3 <- column_3 %>%
  group_by(layer, Cluster) %>%
  summarise(n = sum(n),
            conc = mean(conc, na.rm = TRUE))%>%
  mutate(total_count = sum(n),
         conc_total = sum(conc),
         proportion = conc*100/conc_total)

# contingence tabel for Khi2 test

table_3 <- column_3 %>%
  select(layer, Cluster, n) %>%
  pivot_wider(names_from = layer, values_from = n) %>%
  select(-Cluster) # cluster are order in the alphabetic order (Aggregates, Bright, Elongated, Small)


chisq.test(table_3, correct = FALSE)
chisq.test(table_3,correct=F)$expected # Cochran criterion verified

# There is a difference between layers but now I need to know where 
# khi2 test 2-2 with adjusted p-values 
layer1vs2 <- chisq.test(table_3[, 1:2], correct = FALSE)
p1 <- layer1vs2$p.value

layer2vs3 <- chisq.test(table_3[, 2:3], correct = FALSE)
p2 <- layer2vs3$p.value

layer3vs4 <- chisq.test(table_3[, 3:4], correct = FALSE)
p3 <- layer3vs4$p.value

layer1vs4 <- chisq.test(table_3[, c(1,4)], correct = FALSE)
p4 <- layer1vs4$p.value

layer2vs4 <- chisq.test(table_3[, c(2,4)], correct = FALSE)
p5 <- layer2vs4$p.value

layer1vs3 <- chisq.test(table_3[, c(1,3)], correct = FALSE)
p6 <- layer1vs3$p.value


p_adjusted <- p.adjust(c(p1, p2, p3),method="holm")
p_adjusted

# create the pie chart using ggplot
ggplot(data = column_3, aes(x = "", y = proportion, fill = Cluster)) +
  geom_bar(width = 1, stat = "identity", color= 'white') +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Pastel1") +
  facet_wrap(~layer, ncol = 2) +
  labs(x = NULL, y = NULL, fill = "Cluster") +
  #ggtitle("Proportion of moprhotypes (all individuals) depending on depth in column 3") +
  geom_text(aes(label = paste0(round(proportion), "%")), 
            position = position_stack(vjust = 0.5), color = 'black')+
  theme_void()+
  theme(plot.title.position = "plot",
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"))

plot3 <- ggplot(column_3, aes(fill=Cluster, y=as.factor(layer), x=proportion)) + 
  scale_fill_manual('',values=c("#CD8BFF", "#00BFC4", "#7AAF00", "#F8766D"))+
  geom_bar(position="stack", stat="identity")+
  geom_text(aes(label = paste0(round(column_3$proportion), "%")), 
            position = position_stack(vjust = 0.5), color = 'white', fontface = 'bold')+
  scale_y_discrete(limits=rev)+
  theme_classic()+
  labs(y = '', x = '')+
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=20,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"), 
        axis.text.y = element_blank(), 
        legend.text = element_text(colour = "black", size = 20))
plot3

grid.arrange(plot1, plot2, plot3, nrow = 1, ncol = 3)
ggarrange(plot1, plot2, plot3,
          ncol = 3, nrow = 1, 
          common.legend = TRUE, legend = "top", 
          widths = c(1, 1, 1))
library(gridExtra)
library(patchwork)
library("ggpubr")

List <- list(plot1,plot2,plot3)
#Plot
Plot <- wrap_plots(List,ncol = 3,nrow = 1)
Plot
############################################################################
 
# now try with only selected individuals

morpho_conc <- read_csv('C:/Users/lemoi/Desktop/GIT/k-means/clusters_concentrations_q25.csv') %>% 
  select(profile, depth, Cluster, n, conc) %>%
  filter(depth < 1000)

morpho_conc$Cluster <-  ifelse(morpho_conc$Cluster == "cluster 1", "Bright",
                               ifelse(morpho_conc$Cluster == "cluster 2", "Aggregates",
                                      ifelse(morpho_conc$Cluster == "cluster 3", "Small",
                                             ifelse(morpho_conc$Cluster == "cluster 4", "Elongated",NA))))

# compute layer depth 
morpho_conc$layer <- 
  ifelse(morpho_conc$depth < 100, "0-100m",
         ifelse(morpho_conc$depth >= 100 & morpho_conc$depth < 300, "100-300m",
                ifelse(morpho_conc$depth >= 300 & morpho_conc$depth < 600, "300-600m",
                       ifelse(morpho_conc$depth >= 600 & morpho_conc$depth < 1000, "600-1000m", NA))))

# compute the mean concentration in each layer depth 

morpho_conc_mean <- morpho_conc %>%
  group_by(profile, Cluster, layer) %>%
  summarise(n = sum(n),
            conc = mean(conc, na.rm = TRUE))
# import date, longitude, latitude for each profile 
profile <-  read_csv("C:/Users/lemoi/Desktop/GIT/benguela/Data_created/list_of_profiles.csv")

morpho_conc_mean <- merge(morpho_conc_mean, profile, by = 'profile')

# convert into date format 
morpho_conc_mean$date <- as.Date(as.character(morpho_conc_mean$date), "%Y%m%d")

# select only profiles that are in export columns and compute the mean with all profile that are in the column and represent it 


# import results from STARS method

events_date <- read_csv("C:/Users/lemoi/Desktop/GIT/STARS_method/events_date.csv")

# first export column
start_date <- '2021-10-08'
end_date <- '2021-10-23'

column_1 <- subset(morpho_conc_mean, date >= start_date & date <= end_date)
column_1 <- column_1 %>%
  group_by(layer, Cluster) %>%
  summarise(n = sum(n),
            conc = mean(conc, na.rm = TRUE))%>%
  mutate(total_count = sum(n),
         conc_total = sum(conc),
         proportion = conc*100/conc_total)

# plot proportion according to depth 

plot1 <- ggplot(column_1, aes(fill=Cluster, y=as.factor(layer), x=proportion)) + 
  scale_fill_manual('',values=c("#CD8BFF", "#00BFC4", "#7AAF00", "#F8766D"))+
  geom_bar(position="stack", stat="identity")+
  geom_text(aes(label = paste0(round(column_1$proportion), "%")), 
            position = position_stack(vjust = 0.5), color = 'white', fontface = 'bold')+
  scale_y_discrete(limits=rev)+
  theme_classic()+
  labs(y = 'Layer depth (m)', x = '')+
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=20,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"))

plot1

# second export column
start_date <- '2021-11-25'
end_date <- '2021-12-22'

column_2 <- subset(morpho_conc_mean, date >= start_date & date <= end_date)
column_2 <- column_2 %>%
  group_by(layer, Cluster) %>%
  summarise(n = sum(n),
            conc = mean(conc, na.rm = TRUE))%>%
  mutate(total_count = sum(n),
         conc_total = sum(conc),
         proportion = conc*100/conc_total)

# plot proportion according to depth
plot2 <- ggplot(column_2, aes(fill=Cluster, y=as.factor(layer), x=proportion)) + 
  scale_fill_manual('',values=c("#CD8BFF", "#00BFC4", "#7AAF00", "#F8766D"))+
  geom_bar(position="stack", stat="identity")+
  geom_text(aes(label = paste0(round(column_2$proportion), "%")), 
            position = position_stack(vjust = 0.5), color = 'white', fontface = 'bold')+
  scale_y_discrete(limits=rev)+
  theme_classic()+
  labs(y = '', x = 'Proportion (%)')+
  theme(legend.position = "top",
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=20,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"),
        axis.text.y = element_blank(), 
        legend.text = element_text(colour = "black", size = 20))+
  guides(color = guide_legend(nrow = 2))

plot2

# third export column
start_date <- '2022-03-10'
end_date <- '2022-04-15'

column_3 <- subset(morpho_conc_mean, date >= start_date & date <= end_date)
column_3 <- column_3 %>%
  group_by(layer, Cluster) %>%
  summarise(n = sum(n),
            conc = mean(conc, na.rm = TRUE))%>%
  mutate(total_count = sum(n),
         conc_total = sum(conc),
         proportion = conc*100/conc_total)

# plot proportion according to depth
plot3 <- ggplot(column_3, aes(fill=Cluster, y=as.factor(layer), x=proportion)) + 
  scale_fill_manual('',values=c("#CD8BFF", "#00BFC4", "#7AAF00", "#F8766D"))+
  geom_bar(position="stack", stat="identity")+
  geom_text(aes(label = paste0(round(column_3$proportion), "%")), 
            position = position_stack(vjust = 0.5), color = 'white', fontface = 'bold')+
  scale_y_discrete(limits=rev)+
  theme_classic()+
  labs(y = '', x = '')+
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=20,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"), 
        axis.text.y = element_blank(), 
        legend.text = element_text(colour = "black", size = 20))

plot3


List <- list(plot1,plot2,plot3)
#Plot
Plot <- wrap_plots(List,ncol = 3,nrow = 1)
Plot

############################################################################

# now do the same plots but for one month before each export event 

# first export column (one month before)
start_date <- '2021-09-08'
end_date <- '2021-10-08'

column_1 <- subset(morpho_conc_mean, date >= start_date & date <= end_date)
column_1 <- column_1 %>%
  group_by(layer, Cluster) %>%
  summarise(n = sum(n),
            conc = mean(conc, na.rm = TRUE))%>%
  mutate(total_count = sum(n),
         conc_total = sum(conc),
         proportion = conc*100/conc_total)


plot1 <- ggplot(column_1, aes(fill=Cluster, y=as.factor(layer), x=proportion)) + 
  scale_fill_manual('',values=c("#CD8BFF", "#00BFC4", "#7AAF00", "#F8766D"))+
  geom_bar(position="stack", stat="identity")+
  geom_text(aes(label = paste0(round(column_1$proportion), "%")), 
            position = position_stack(vjust = 0.5), color = 'white', fontface = 'bold')+
  scale_y_discrete(limits=rev)+
  theme_classic()+
  labs(y = 'Layer depth (m)', x = '')+
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=20,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"))

plot1

# second export column
start_date <- '2021-10-25'
end_date <- '2021-11-25'

column_2 <- subset(morpho_conc_mean, date >= start_date & date <= end_date)
column_2 <- column_2 %>%
  group_by(layer, Cluster) %>%
  summarise(n = sum(n),
            conc = mean(conc, na.rm = TRUE))%>%
  mutate(total_count = sum(n),
         conc_total = sum(conc),
         proportion = conc*100/conc_total)

# create the pie chart using ggplot
plot2 <- ggplot(column_2, aes(fill=Cluster, y=as.factor(layer), x=proportion)) + 
  scale_fill_manual('',values=c("#CD8BFF", "#00BFC4", "#7AAF00", "#F8766D"))+
  geom_bar(position="stack", stat="identity")+
  geom_text(aes(label = paste0(round(column_2$proportion), "%")), 
            position = position_stack(vjust = 0.5), color = 'white', fontface = 'bold')+
  scale_y_discrete(limits=rev)+
  theme_classic()+
  labs(y = '', x = 'Proportion (%)')+
  theme(legend.position = "top",
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=20,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"),
        axis.text.y = element_blank(), 
        legend.text = element_text(colour = "black", size = 20))+
  guides(color = guide_legend(nrow = 2))

plot2

# third export column
start_date <- '2022-02-10'
end_date <- '2022-03-10'

column_3 <- subset(morpho_conc_mean, date >= start_date & date <= end_date)
column_3 <- column_3 %>%
  group_by(layer, Cluster) %>%
  summarise(n = sum(n),
            conc = mean(conc, na.rm = TRUE))%>%
  mutate(total_count = sum(n),
         conc_total = sum(conc),
         proportion = conc*100/conc_total)

# create the pie chart using ggplot
plot3 <- ggplot(column_3, aes(fill=Cluster, y=as.factor(layer), x=proportion)) + 
  scale_fill_manual('',values=c("#CD8BFF", "#00BFC4", "#7AAF00", "#F8766D"))+
  geom_bar(position="stack", stat="identity")+
  geom_text(aes(label = paste0(round(column_3$proportion), "%")), 
            position = position_stack(vjust = 0.5), color = 'white', fontface = 'bold')+
  scale_y_discrete(limits=rev)+
  theme_classic()+
  labs(y = '', x = '')+
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=20,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"), 
        axis.text.y = element_blank(), 
        legend.text = element_text(colour = "black", size = 20))

plot3


List <- list(plot1,plot2,plot3)
#Plot
Plot <- wrap_plots(List,ncol = 3,nrow = 1)
Plot


# now try with profile after each export event 

# first export column
start_date <- '2021-10-23'
end_date <- '2021-11-23'

column_1 <- subset(morpho_conc_mean, date >= start_date & date <= end_date)
column_1 <- column_1 %>%
  group_by(layer, Cluster) %>%
  summarise(n = sum(n),
            conc = mean(conc, na.rm = TRUE))%>%
  mutate(total_count = sum(n),
         conc_total = sum(conc),
         proportion = conc*100/conc_total)


plot1 <- ggplot(column_1, aes(fill=Cluster, y=as.factor(layer), x=proportion)) + 
  scale_fill_manual('',values=c("#CD8BFF", "#00BFC4", "#7AAF00", "#F8766D"))+
  geom_bar(position="stack", stat="identity")+
  geom_text(aes(label = paste0(round(column_1$proportion), "%")), 
            position = position_stack(vjust = 0.5), color = 'white', fontface = 'bold')+
  scale_y_discrete(limits=rev)+
  theme_classic()+
  labs(y = 'Layer depth (m)', x = '')+
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=20,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"))

plot1

# second export column
start_date <- '2021-12-22'
end_date <- '2022-01-22'

column_2 <- subset(morpho_conc_mean, date >= start_date & date <= end_date)
column_2 <- column_2 %>%
  group_by(layer, Cluster) %>%
  summarise(n = sum(n),
            conc = mean(conc, na.rm = TRUE))%>%
  mutate(total_count = sum(n),
         conc_total = sum(conc),
         proportion = conc*100/conc_total)

# create the pie chart using ggplot
plot2 <- ggplot(column_2, aes(fill=Cluster, y=as.factor(layer), x=proportion)) + 
  scale_fill_manual('',values=c("#CD8BFF", "#00BFC4", "#7AAF00", "#F8766D"))+
  geom_bar(position="stack", stat="identity")+
  geom_text(aes(label = paste0(round(column_2$proportion), "%")), 
            position = position_stack(vjust = 0.5), color = 'white', fontface = 'bold')+
  scale_y_discrete(limits=rev)+
  theme_classic()+
  labs(y = '', x = 'Proportion (%)')+
  theme(legend.position = "top",
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=20,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"),
        axis.text.y = element_blank(), 
        legend.text = element_text(colour = "black", size = 20))+
  guides(color = guide_legend(nrow = 2))

plot2

# third export column
start_date <- '2022-04-15'
end_date <- '2022-05-15'

column_3 <- subset(morpho_conc_mean, date >= start_date & date <= end_date)
column_3 <- column_3 %>%
  group_by(layer, Cluster) %>%
  summarise(n = sum(n),
            conc = mean(conc, na.rm = TRUE))%>%
  mutate(total_count = sum(n),
         conc_total = sum(conc),
         proportion = conc*100/conc_total)

# create the pie chart using ggplot
plot3 <- ggplot(column_3, aes(fill=Cluster, y=as.factor(layer), x=proportion)) + 
  scale_fill_manual('',values=c("#CD8BFF", "#00BFC4", "#7AAF00", "#F8766D"))+
  geom_bar(position="stack", stat="identity")+
  geom_text(aes(label = paste0(round(column_3$proportion), "%")), 
            position = position_stack(vjust = 0.5), color = 'white', fontface = 'bold')+
  scale_y_discrete(limits=rev)+
  theme_classic()+
  labs(y = '', x = '')+
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=20,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"), 
        axis.text.y = element_blank(), 
        legend.text = element_text(colour = "black", size = 20))

plot3


List <- list(plot1,plot2,plot3)
#Plot
Plot <- wrap_plots(List,ncol = 3,nrow = 1)
Plot

##########################################
# now try to remove some profiles that are supected to not be in a physical column 
##########################################

# first export column
start_date <- '2021-10-08'
end_date <- '2021-10-23'

column_1 <- subset(morpho_conc_mean, date >= start_date & date <= end_date)
column_1 <- column_1 %>%
  group_by(layer, Cluster) %>%
  summarise(n = sum(n),
            conc = mean(conc, na.rm = TRUE))%>%
  mutate(total_count = sum(n),
         conc_total = sum(conc),
         proportion = conc*100/conc_total)

# plot proportion according to depth 

plot1 <- ggplot(column_1, aes(fill=Cluster, y=as.factor(layer), x=proportion)) + 
  scale_fill_manual('',values=c("#CD8BFF", "#00BFC4", "#7AAF00", "#F8766D"))+
  geom_bar(position="stack", stat="identity")+
  geom_text(aes(label = paste0(round(column_1$proportion), "%")), 
            position = position_stack(vjust = 0.5), color = 'white', fontface = 'bold')+
  scale_y_discrete(limits=rev)+
  theme_classic()+
  labs(y = 'Layer depth (m)', x = '')+
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=20,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"))

plot1

# second export column
start_date <- '2021-11-25'
end_date <- '2021-12-22'

column_2 <- subset(morpho_conc_mean, date >= start_date & date <= end_date) %>%
  filter(profile != "0081a_WMO6903095" &
         profile != "0082a_WMO6903095" &
         profile != "0083a_WMO6903095")

column_2 <- column_2 %>%
  group_by(layer, Cluster) %>%
  summarise(n = sum(n),
            conc = mean(conc, na.rm = TRUE))%>%
  mutate(total_count = sum(n),
         conc_total = sum(conc),
         proportion = conc*100/conc_total)

# plot proportion according to depth
plot2 <- ggplot(column_2, aes(fill=Cluster, y=as.factor(layer), x=proportion)) + 
  scale_fill_manual('',values=c("#CD8BFF", "#00BFC4", "#7AAF00", "#F8766D"))+
  geom_bar(position="stack", stat="identity")+
  geom_text(aes(label = paste0(round(column_2$proportion), "%")), 
            position = position_stack(vjust = 0.5), color = 'white', fontface = 'bold')+
  scale_y_discrete(limits=rev)+
  theme_classic()+
  labs(y = '', x = 'Proportion (%)')+
  theme(legend.position = "top",
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=20,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"),
        axis.text.y = element_blank(), 
        legend.text = element_text(colour = "black", size = 20))+
  guides(color = guide_legend(nrow = 2))

plot2

# third export column
start_date <- '2022-03-10'
end_date <- '2022-04-15'

column_3 <- subset(morpho_conc_mean, date >= start_date & date <= end_date) %>%
  filter(profile != "0116a_WMO6903095" &
           profile != "0117a_WMO6903095" &
           profile != "0121a_WMO6903095" &
           profile != "0122a_WMO6903095" &
           profile != "0123a_WMO6903095" &
           profile != "0124a_WMO6903095" &
           profile != "0127a_WMO6903095")


column_3 <- column_3 %>%
  group_by(layer, Cluster) %>%
  summarise(n = sum(n),
            conc = mean(conc, na.rm = TRUE))%>%
  mutate(total_count = sum(n),
         conc_total = sum(conc),
         proportion = conc*100/conc_total)

# plot proportion according to depth
plot3 <- ggplot(column_3, aes(fill=Cluster, y=as.factor(layer), x=proportion)) + 
  scale_fill_manual('',values=c("#CD8BFF", "#00BFC4", "#7AAF00", "#F8766D"))+
  geom_bar(position="stack", stat="identity")+
  geom_text(aes(label = paste0(round(column_3$proportion), "%")), 
            position = position_stack(vjust = 0.5), color = 'white', fontface = 'bold')+
  scale_y_discrete(limits=rev)+
  theme_classic()+
  labs(y = '', x = '')+
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=20,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"), 
        axis.text.y = element_blank(), 
        legend.text = element_text(colour = "black", size = 20))

plot3


List <- list(plot1,plot2,plot3)
#Plot
Plot <- wrap_plots(List,ncol = 3,nrow = 1)
Plot


  