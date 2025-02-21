# STARS method to detect shift in particles TS  

# libraries 
library(rshift)
library(ggplot2)
library(readr)
library(tidyverse)
library(dplyr)
library(scales)
library(patchwork)


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


# select only MaP first 

MaP <- ecopart %>% select(Profile, Date_Time, Depth, MaP_abun)

# try first to integrate all the water column for each profile 

MaP_all <- MaP %>%
  group_by(Profile)%>%
  summarise(Date_Time = mean(Date_Time, na.rm = TRUE),
            MaP_abun = mean(MaP_abun, na.rm = TRUE))

# STARS method 

MaP_STARS <- Rodionov(MaP_all, "MaP_abun", "Date_Time", 10, merge = TRUE)

# Filter the data for RSI > 0
shift_values <- MaP_STARS[MaP_STARS$RSI > 0,]

ggplot(MaP_STARS)+
  geom_path(aes(x = Date_Time, y = MaP_abun))+
  geom_vline(data = shift_values, aes(xintercept = as.numeric(Date_Time)), linetype = "dashed")


# try to integrate MaP abundance between 150 and 600m for each profile 

MaP_300_600 <- MaP %>%
  filter(Depth > 300, Depth < 600)%>%
  group_by(Profile)%>%
  summarise(Date_Time = mean(Date_Time, na.rm = TRUE),
            MaP_abun = mean(MaP_abun, na.rm = TRUE))

# STARS method 

MaP_STARS <- Rodionov(MaP_300_600, "MaP_abun", "Date_Time", 10, merge = TRUE)

mean_MaP <- mean(MaP_STARS$MaP_abun, na.rm = TRUE)

# Filter the data for RSI > 0
shift_values <- MaP_STARS[MaP_STARS$RSI > 0,]

#select only profiles that are in the periods 
first_period <- subset(MaP_STARS, Date_Time >= '2021-10-08' & Date_Time <= '2021-10-23')
second_period <- subset(MaP_STARS, Date_Time >= '2021-11-25' & Date_Time <= '2021-12-22')
third_period <- subset(MaP_STARS, Date_Time >= '2022-03-10' & Date_Time <= '2022-04-15')

#select profiles one month before each event  
before_first_period <- subset(MaP_STARS, Date_Time >= '2021-09-08' & Date_Time <= '2021-10-08')
before_second_period <- subset(MaP_STARS, Date_Time >= '2021-10-25' & Date_Time <= '2021-11-25')
before_third_period <- subset(MaP_STARS, Date_Time >= '2022-02-10' & Date_Time <= '2022-03-10')

#select profiles one month after each event  
after_first_period <- subset(MaP_STARS, Date_Time >= '2021-10-23' & Date_Time <= '2021-11-23')
after_second_period <- subset(MaP_STARS, Date_Time >= '2021-12-22' & Date_Time <= '2022-01-22')
after_third_period <- subset(MaP_STARS, Date_Time >= '2022-04-15' & Date_Time <= '2022-05-15')



profiles_in <- rbind(first_period, second_period, third_period)

# select profiles out 
profiles_out <- subset(MaP_STARS, Date_Time < '2021-10-08'| 
                                  Date_Time > '2021-10-23' & Date_Time < '2021-11-25'|
                                  Date_Time > '2021-12-22' & Date_Time < '2022-03-10'|
                                  Date_Time > '2022-04-15')

# save it 
write.csv(shift_values, "/Users/alexandreaccardo/Library/CloudStorage/OneDrive-Personnel/Stage M2/GIT/STARS_method/events_date.csv", row.names = FALSE)

df <- read_csv("/Users/alexandreaccardo/Library/CloudStorage/OneDrive-Personnel/Stage M2/GIT/STARS_method/events_date.csv")

#import FTLE values 

profiles_date <- select(MaP_all, Profile, Date_Time) %>%
  mutate(CYCLE_NUMBER = c(1:184))

FTLE <- read_csv('/Users/alexandreaccardo/Library/CloudStorage/OneDrive-Personnel/Stage M2/GIT/FTLE.csv') %>% 
  mutate(Profile = profile) %>%
  select(-profile) %>%
  merge(MaP_STARS, on = Profile) %>%
  rename(FTLE = Ftle_GlobEkmanDt_005daysBackward_mean_delta0ftle010)

# plot

ggplot(MaP_STARS)+
  geom_path(aes(x = Date_Time, y = MaP_abun))+
  geom_vline(data = shift_values, aes(xintercept = as.numeric(Date_Time)), linetype = "dashed")

plot1 <- ggplot()+
  geom_path(data = MaP_STARS, aes(x = Date_Time, y = MaP_abun), alpha = 0.5)+
  geom_point(data = profiles_out, aes(x = Date_Time, y = MaP_abun), color = 'black')+
  geom_point(data = first_period, aes(x = Date_Time, y = MaP_abun), color = 'blue')+
  geom_point(data = second_period, aes(x = Date_Time, y = MaP_abun), color = 'red')+
  geom_point(data = third_period, aes(x = Date_Time, y = MaP_abun), color = 'green')+
  geom_vline(data = shift_values, aes(xintercept = as.numeric(Date_Time)), linetype = "dashed")+
  #geom_hline(yintercept = mean_MaP)
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=20,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"))
plot1

plot2 <- ggplot()+
  geom_line(data = FTLE, aes(x = Date_Time, y = FTLE), color = 'black')+
  geom_point(data = profiles_out, aes(x = Date_Time, y = MaP_abun), color = 'grey')+
  geom_point(data = first_period, aes(x = Date_Time, y = MaP_abun), color = 'blue')+
  geom_point(data = second_period, aes(x = Date_Time, y = MaP_abun), color = 'red')+
  geom_point(data = third_period, aes(x = Date_Time, y = MaP_abun), color = 'green')+
  #geom_vline(data = shift_values, aes(xintercept = as.numeric(Date_Time)), linetype = "dashed")
  theme_bw()+
  labs(y = '', x = 'Date')+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=20,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"))
plot2


#import AOU values 

AOU <- read_csv("/Users/alexandreaccardo/Library/CloudStorage/OneDrive-Personnel/Stage M2/GIT/benguela/Data_created/AOU.csv")
AOU_300_600 <- AOU %>%
  merge(profiles_date, on = 'CYCLE_NUMBER')%>%
  filter(`Depth [m]` > 300, `Depth [m]` < 600)%>%
  group_by(CYCLE_NUMBER)%>%
  summarise(Date_Time = mean(Date_Time, na.rm = TRUE),
            AOU = mean(`AOU µmol/kg`, na.rm = TRUE))

#import Temperature and salinity values 

TS <- read_csv("/Users/alexandreaccardo/Library/CloudStorage/OneDrive-Personnel/Stage M2/GIT/raw_bgc_data/bgc_benguela_clean.csv")

TS_300_600 <- TS %>%
  merge(profiles_date, on = 'CYCLE_NUMBER')%>%
  filter(binned_depth > 300, binned_depth < 600)%>%
  group_by(CYCLE_NUMBER)%>%
  summarise(Date_Time = mean(Date_Time, na.rm = TRUE),
            Temp = mean(TEMP, na.rm = TRUE),
            Sal = mean(PSAL, na.rm = TRUE))

library(patchwork)
subplot <- plot2 / plot1
subplot



# now compute characterize each event 

# compute anomalies 

MaP_abun_mean <- mean(MaP_300_600$MaP_abun, na.rm = TRUE)
MaP_abun_sd <- sd(MaP_300_600$MaP_abun, na.rm = TRUE)

first_period <- first_period %>%
  mutate(event_mean = mean(MaP_abun, na.rm = TRUE),
         event_sd = sd(MaP_abun, na.rm = TRUE), 
         event_anomaly = (MaP_abun - MaP_abun_mean)/MaP_abun_sd)

before_first_period <- before_first_period %>%
  mutate(event_mean = mean(MaP_abun, na.rm = TRUE),
         event_sd = sd(MaP_abun, na.rm = TRUE), 
         event_anomaly = (MaP_abun - MaP_abun_mean)/MaP_abun_sd)

after_first_period <- after_first_period %>%
  mutate(event_mean = mean(MaP_abun, na.rm = TRUE),
         event_sd = sd(MaP_abun, na.rm = TRUE), 
         event_anomaly = (MaP_abun - MaP_abun_mean)/MaP_abun_sd)

second_period <- second_period %>%
  mutate(event_mean = mean(MaP_abun, na.rm = TRUE),
         event_sd = sd(MaP_abun, na.rm = TRUE), 
         event_anomaly = (MaP_abun - MaP_abun_mean)/MaP_abun_sd)

before_second_period <- before_second_period %>%
  mutate(event_mean = mean(MaP_abun, na.rm = TRUE),
         event_sd = sd(MaP_abun, na.rm = TRUE), 
         event_anomaly = (MaP_abun - MaP_abun_mean)/MaP_abun_sd)

after_second_period <- after_second_period %>%
  mutate(event_mean = mean(MaP_abun, na.rm = TRUE),
         event_sd = sd(MaP_abun, na.rm = TRUE), 
         event_anomaly = (MaP_abun - MaP_abun_mean)/MaP_abun_sd)

third_period <- third_period %>%
  mutate(event_mean = mean(MaP_abun, na.rm = TRUE),
         event_sd = sd(MaP_abun, na.rm = TRUE), 
         event_anomaly = (MaP_abun - MaP_abun_mean)/MaP_abun_sd)

before_third_period <- before_third_period %>%
  mutate(event_mean = mean(MaP_abun, na.rm = TRUE),
         event_sd = sd(MaP_abun, na.rm = TRUE), 
         event_anomaly = (MaP_abun - MaP_abun_mean)/MaP_abun_sd)

after_third_period <- after_third_period %>%
  mutate(event_mean = mean(MaP_abun, na.rm = TRUE),
         event_sd = sd(MaP_abun, na.rm = TRUE), 
         event_anomaly = (MaP_abun - MaP_abun_mean)/MaP_abun_sd)

p1 <- ggplot()+
  geom_path(data = MaP_STARS, aes(x = Date_Time, y = MaP_abun), alpha = 0.5)+
  geom_point(data = profiles_out, aes(x = Date_Time, y = MaP_abun), color = 'black')+
  geom_hline(yintercept = MaP_abun_mean)+
  geom_point(data = first_period, aes(x = Date_Time, y = MaP_abun), color = 'blue')+
  geom_segment(aes(x = as.Date("2021-10-08"), y = first_period$event_mean, xend = as.Date("2021-10-23"), yend = first_period$event_mean), color = "blue", size = 1.5)+
  geom_point(data = second_period, aes(x = Date_Time, y = MaP_abun), color = 'red')+
  geom_segment(aes(x = as.Date("2021-11-25"), y = second_period$event_mean, xend = as.Date("2021-12-22"), yend = second_period$event_mean), color = "red", size = 1.5)+
  geom_point(data = third_period, aes(x = Date_Time, y = MaP_abun), color = 'green')+
  geom_segment(aes(x = as.Date("2022-03-10"), y = third_period$event_mean, xend = as.Date("2022-04-15"), yend = third_period$event_mean), color = "green", size = 1.5)+
  #geom_vline(data = shift_values, aes(xintercept = as.numeric(Date_Time)), linetype = "dashed")+
  theme_bw()+
  labs(y = expression("MaP concentration (part.m"^-3*")"), x = 'Date')+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=20,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"))

ggsave(filename = "/Users/alexandreaccardo/Library/CloudStorage/OneDrive-Personnel/Stage M2/GIT/Figures/stars__events_detection.png", plot = p1, width = 10, height = 6, units = "in", dpi = 300)

# zoom on first period 

    # Map 
    
p1 <- ggplot()+
  geom_path(data = MaP_STARS, aes(x = Date_Time, y = MaP_abun), alpha = 0.5)+
  geom_point(data = before_first_period, aes(x = Date_Time, y = MaP_abun), color = 'black')+
  geom_segment(aes(x = as.Date("2021-09-08"), y = before_first_period$event_mean, xend = as.Date("2021-10-08"), yend = before_first_period$event_mean), color = "black", size = 1.5)+
  geom_point(data = first_period, aes(x = Date_Time, y = MaP_abun), color = 'blue')+
  geom_segment(aes(x = as.Date("2021-10-08"), y = first_period$event_mean, xend = as.Date("2021-10-23"), yend = first_period$event_mean), color = "blue", size = 1.5)+
  geom_point(data = after_first_period, aes(x = Date_Time, y = MaP_abun), color = 'black')+
  geom_segment(aes(x = as.Date("2021-10-23"), y = after_first_period$event_mean, xend = as.Date("2021-11-23"), yend = after_first_period$event_mean), color = "black", size = 1.5)+
  labs(y = 'MaP abundance (particles/m3)', x = 'Date')+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=10,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"))+
  scale_x_date(labels = date_format("%Y-%m"), limits = c(as.Date("2021-09-08"), as.Date("2021-11-23")))
  
p1
    # AOU 
p2 <- ggplot()+
  geom_path(data = AOU_300_600, aes(x = Date_Time, y = AOU), alpha = 0.5)+
  geom_point(data = AOU_300_600, aes(x = Date_Time, y = AOU), color = 'black')+
  labs(y = 'AOU (µmol/kg)', x = '')+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=10,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_blank())+
  scale_x_date(labels = date_format("%Y-%m"), limits = c(as.Date("2021-09-08"), as.Date("2021-11-23")))

    
    # FTLE 
p3 <- ggplot()+
  geom_path(data = FTLE, aes(x = Date_Time, y = FTLE), alpha = 0.5)+
  geom_point(data = FTLE, aes(x = Date_Time, y = FTLE), color = 'black')+
  labs(y = 'FTLE (days-1)', x = '')+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=10,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_blank())+
  scale_x_date(labels = date_format("%Y-%m"), limits = c(as.Date("2021-09-08"), as.Date("2021-11-23")))

 # Temp 
p4 <- ggplot()+
  geom_path(data = TS_300_600, aes(x = Date_Time, y = Temp), alpha = 0.5)+
  geom_point(data = TS_300_600, aes(x = Date_Time, y = Temp), color = 'black')+
  labs(y = 'Temperature (°C)', x = '')+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=10,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_blank())+
  scale_x_date(labels = date_format("%Y-%m"), limits = c(as.Date("2021-09-08"), as.Date("2021-11-23")))

    # Salinity 
p5 <- ggplot()+
  geom_path(data = TS_300_600, aes(x = Date_Time, y = Sal), alpha = 0.5)+
  geom_point(data = TS_300_600, aes(x = Date_Time, y = Sal), color = 'black')+
  labs(y = 'Salinity (psu)', x = '')+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=10,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_blank())+
  scale_x_date(labels = date_format("%Y-%m"), limits = c(as.Date("2021-09-08"), as.Date("2021-11-23")))


subplot <- p5 / p4 / p3 / p2 / p1
subplot

# zoom on second period 

    # Map

p1 <- ggplot()+
  geom_path(data = MaP_STARS, aes(x = Date_Time, y = MaP_abun), alpha = 0.5)+
  geom_point(data = before_second_period, aes(x = Date_Time, y = MaP_abun), color = 'black')+
  geom_segment(aes(x = as.Date("2021-10-25"), y = before_second_period$event_mean, xend = as.Date("2021-11-25"), yend = before_second_period$event_mean), color = "black", size = 1.5)+
  geom_point(data = second_period, aes(x = Date_Time, y = MaP_abun), color = 'blue')+
  geom_segment(aes(x = as.Date("2021-11-25"), y = second_period$event_mean, xend = as.Date("2021-12-22"), yend = second_period$event_mean), color = "blue", size = 1.5)+
  geom_point(data = after_second_period, aes(x = Date_Time, y = MaP_abun), color = 'black')+
  geom_segment(aes(x = as.Date("2021-12-22"), y = after_second_period$event_mean, xend = as.Date("2022-01-22"), yend = after_second_period$event_mean), color = "black", size = 1.5)+
  theme_bw()+
  labs(y = 'MaP abundance (particles/m3)', x = '')+
  xlim(as.Date("2021-10-25"), as.Date("2022-01-22"))+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=10,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"))+
  scale_x_date(labels = date_format("%Y-%m"), limits = c(as.Date("2021-10-25"), as.Date("2022-01-22")))
p1
    # AOU 
p2 <- ggplot()+
  geom_path(data = AOU_300_600, aes(x = Date_Time, y = AOU), alpha = 0.5)+
  geom_point(data = AOU_300_600, aes(x = Date_Time, y = AOU), color = 'black')+
  labs(y = 'AOU (µmol/kg)', x = '')+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=10,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_blank())+
  scale_x_date(labels = date_format("%Y-%m"), limits = c(as.Date("2021-10-25"), as.Date("2022-01-22")))

    # FTLE 
p3 <- ggplot()+
  geom_path(data = FTLE, aes(x = Date_Time, y = FTLE), alpha = 0.5)+
  geom_point(data = FTLE, aes(x = Date_Time, y = FTLE), color = 'black')+
  labs(y = 'FTLE (days-1)', x = '')+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=10,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_blank())+
  scale_x_date(labels = date_format("%Y-%m"), limits = c(as.Date("2021-10-25"), as.Date("2022-01-22")))

    # Temp 
p4 <- ggplot()+
  geom_path(data = TS_300_600, aes(x = Date_Time, y = Temp), alpha = 0.5)+
  geom_point(data = TS_300_600, aes(x = Date_Time, y = Temp), color = 'black')+
  labs(y = 'Temperature (°C)', x = '')+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=10,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_blank())+
  scale_x_date(labels = date_format("%Y-%m"), limits = c(as.Date("2021-10-25"), as.Date("2022-01-22")))

    # Salinity 
p5 <- ggplot()+
  geom_path(data = TS_300_600, aes(x = Date_Time, y = Sal), alpha = 0.5)+
  geom_point(data = TS_300_600, aes(x = Date_Time, y = Sal), color = 'black')+
  labs(y = 'Salinity (psu)', x = '')+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=10,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_blank())+
  scale_x_date(labels = date_format("%Y-%m"), limits = c(as.Date("2021-10-25"), as.Date("2022-01-22")))


subplot <- p5 / p4 / p3 / p2 / p1
subplot

# zoom on third period 

p1 <- ggplot()+
  geom_path(data = MaP_STARS, aes(x = Date_Time, y = MaP_abun), alpha = 0.5)+
  geom_point(data = before_third_period, aes(x = Date_Time, y = MaP_abun), color = 'black')+
  geom_segment(aes(x = as.Date("2022-02-10"), y = before_third_period$event_mean, xend = as.Date("2022-03-10"), yend = before_third_period$event_mean), color = "black", size = 1.5)+
  geom_point(data = third_period, aes(x = Date_Time, y = MaP_abun), color = 'blue')+
  geom_segment(aes(x = as.Date("2022-03-10"), y = third_period$event_mean, xend = as.Date("2022-04-15"), yend = third_period$event_mean), color = "blue", size = 1.5)+
  geom_point(data = after_third_period, aes(x = Date_Time, y = MaP_abun), color = 'black')+
  geom_segment(aes(x = as.Date("2022-04-15"), y = after_third_period$event_mean, xend = as.Date("2022-05-15"), yend = after_third_period$event_mean), color = "black", size = 1.5)+
  theme_bw()+
  labs(y = 'MaP abundance (particles/m3)', x = 'Date')+
  xlim(as.Date("2022-02-10"), as.Date("2022-05-15"))+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=10,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"))+
  scale_x_date(labels = date_format("%Y-%m"), limits = c(as.Date("2022-02-10"), as.Date("2022-05-15")))

# AOU 
p2 <- ggplot()+
  geom_path(data = AOU_300_600, aes(x = Date_Time, y = AOU), alpha = 0.5)+
  geom_point(data = AOU_300_600, aes(x = Date_Time, y = AOU), color = 'black')+
  labs(y = 'AOU (µmol/kg)', x = '')+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=10,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_blank())+
  scale_x_date(labels = date_format("%Y-%m"), limits = c(as.Date("2022-02-10"), as.Date("2022-05-15")))

# FTLE 
p3 <- ggplot()+
  geom_path(data = FTLE, aes(x = Date_Time, y = FTLE), alpha = 0.5)+
  geom_point(data = FTLE, aes(x = Date_Time, y = FTLE), color = 'black')+
  labs(y = 'FTLE (days-1)', x = '')+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=10,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_blank())+
  scale_x_date(labels = date_format("%Y-%m"), limits = c(as.Date("2022-02-10"), as.Date("2022-05-15")))

  # Temp 
p4 <- ggplot()+
  geom_path(data = TS_300_600, aes(x = Date_Time, y = Temp), alpha = 0.5)+
  geom_point(data = TS_300_600, aes(x = Date_Time, y = Temp), color = 'black')+
  labs(y = 'Temperature (°C)', x = '')+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=10,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_blank())+
  scale_x_date(labels = date_format("%Y-%m"), limits = c(as.Date("2022-02-10"), as.Date("2022-05-15")))

  # Salinity 
p5 <- ggplot()+
  geom_path(data = TS_300_600, aes(x = Date_Time, y = Sal), alpha = 0.5)+
  geom_point(data = TS_300_600, aes(x = Date_Time, y = Sal), color = 'black')+
  labs(y = 'Salinity (psu)', x = '')+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.title=element_text(size=10,face="bold"),
        axis.text = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_blank())+
  scale_x_date(labels = date_format("%Y-%m"), limits = c(as.Date("2022-02-10"), as.Date("2022-05-15")))

subplot <- p5 / p4 / p3 / p2 / p1
subplot
