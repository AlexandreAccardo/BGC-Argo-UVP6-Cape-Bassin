# Script to detect spikes as big particles 
# Method developed by Briggs et al., 2020 (nature)


library(zoo)
library(tidyverse)
library(readr)
library(ggplot2)
library(viridis)

Briggs_2020_BBP <- function(BBP700, PRES, cycleNumber){
  
  final <- data.frame(BBP700, PRES, cycleNumber)
  
  final <- final %>%
    group_by(cycleNumber) %>%
    #removes spikes : filtered = bbsr
    mutate(filtered = rollapply(BBP700, 11, min, fill = "extend")) %>%
    mutate(filtered = rollmax(filtered, 11, fill = "extend"),
           residual_spikes = BBP700 - filtered) %>% # residual spikes = bbl + noise
    ungroup()
  
  # determination du bruit de chaque flotteur
  NOISE <-  final %>%
    filter(PRES>300) %>%
    mutate(mediane = median(residual_spikes)) %>% # mediane de tout les profils d'un flotteur sous 300m
    # filter((bbl_noise<0)==F) %>%
    mutate (bin = NA)
  
  # rattribuer des bins par tranche de 50m sous 300m
  for(j in 1:length(NOISE$BBP700)){
    for(l in 1:14){
      if(NOISE$PRES[j]<(300+l*50) & NOISE$PRES[j]>=(300+(l-1)*50)){
        NOISE$bin[j] <- l
      }
    }
  }
  
  NOISE <- NOISE %>%
    dplyr::group_by(cycleNumber, bin) %>%
    dplyr::summarise(residual_spikes = residual_spikes, mediane = mediane, moyenne = mean(residual_spikes), test = moyenne < (2*mediane)) %>%
    ungroup() %>%
    filter(test = TRUE) %>% # retirer si moyenne superieure a 2 fois la mediane
    dplyr::summarise(mediane_final = median(residual_spikes)) %>% # calcul de nouveu la mediane sans les spikes = bruit final du flotteur
    as.numeric()
  
  # bbrs et bbl ajoutés au fichier final
  final <- final %>%
    mutate(noise_float = NOISE,
           bbsr = filtered + NOISE,
           bbl = ifelse(BBP700-bbsr>0, BBP700-bbsr,0))
  
  # bbr = blanc profond determiné à partir du bbsr de fond
  bbr <- final %>%
    filter(PRES>850 & PRES<900) %>%
    dplyr::summarise(quantile(bbsr,0.25, na.rm=T))%>%
    as.numeric()
  
  #bbs ajouté au fichier final
  final <- final %>%
    mutate(bbr = bbr,
           bbs = bbsr - bbr)
  
  filtered <- final %>% pull(filtered)
  residual_spikes <- final %>% pull(residual_spikes)
  noise_float <- final %>% pull(noise_float)
  bbsr <- final %>% pull(bbsr)
  bbr <- final %>% pull(bbr)
  bbl <- final %>% pull(bbl)
  bbs <- final %>% pull(bbs)
  
  return(list(filtered_bbp = filtered, residual_spikes_bbp = residual_spikes, noise_bbp = noise_float, bbsr= bbsr, bbr= bbr, bbs=bbs, bbl=bbl))
}

bgc_benguela <- read_csv("C:/Users/lemoi/Desktop/GIT/raw_bgc_data/bgc_benguela.csv") %>% 
  select(CYCLE_NUMBER, PRES, BBP700) %>%
  filter(BBP700 != 'NA')

bgc_benguela <- bgc_benguela %>% 
  mutate(bbsr = Briggs_2020_BBP(bgc_benguela$BBP700, 
                                bgc_benguela$PRES, 
                                cycleNumber = bgc_benguela$CYCLE_NUMBER)$bbsr,
         bbl = Briggs_2020_BBP(bgc_benguela$BBP700, 
                               bgc_benguela$PRES, 
                               cycleNumber = bgc_benguela$CYCLE_NUMBER)$bbl)

# try to represent selected profiles to see how it's look like 
bbp <- bgc_benguela %>%
  filter(CYCLE_NUMBER ==64)

bbp %>%
  ggplot(aes(x = PRES, y = BBP700))+
  geom_point()+
  geom_line(aes(x = PRES, y = bbsr), color = 'green')+
  geom_line(aes(x = PRES, y = bbl), color = 'red')+ 
  scale_x_reverse() +
  scale_y_continuous() +
  coord_flip()+
  theme_bw()


# time series of bbp

# import csv file with date of each profile and merge it to bbp data 

profiles <-  read_csv("C:/Users/lemoi/Desktop/GIT/benguela/benguela_profiles.csv")

profiles <- profiles[grep("a", profiles$Profile), ] %>%
  select(Profile, date) %>% 
  mutate(profile = Profile, CYCLE_NUMBER = c(1:184)) %>%
  select(-Profile)

# merge it and select only depth > 1000 m and log transfrom
bgc_benguela <- bgc_benguela %>%
  merge(profiles, by = 'CYCLE_NUMBER')%>%
  select(-profile) %>%
  filter(PRES < 1000) %>%
  mutate(log_bbl = log10(bbl+1), 
         log_bbsr = log10(bbsr))

write.csv(bgc_benguela, "C:/Users/lemoi/Desktop/GIT/raw_bgc_data/bbp_decomposed.csv", row.names = FALSE)

# represent time series

ggplot(bgc_benguela, aes(x = date, y = -PRES, color = log_bbsr+10)) +
  geom_point(size = 3) +
  scale_color_viridis() +#limits = c(min(bgc_benguela$log_bbsr+10), 6.7), na.value = "yellow")+
  # set up the upper limit of the color bar to see spikes that are less intense
  #scale_color_viridis(limits = c(0, quantile(bgc_benguela$bbl, 0.95)), na.value = "yellow") + 
  
  #represent period of column 1 
  geom_vline(xintercept = as.numeric(as.POSIXct("2021-10-01")), 
             linetype = "dashed", color = "white", size = 1.5)+
  geom_vline(xintercept = as.numeric(as.POSIXct("2021-10-17")), 
             linetype = "dashed", color = "white", size = 1.5)+
  annotate("text", label = "MaP Column 1", x = as.POSIXct("2021-10-10"), 
           y = 20, size = 3.5, colour = "black", fontface = "bold.italic") +
  
  #represent period of column 2
  geom_vline(xintercept = as.numeric(as.POSIXct("2021-12-01")), 
             linetype = "dashed", color = "white", size = 1.5)+
  geom_vline(xintercept = as.numeric(as.POSIXct("2021-12-19")), 
             linetype = "dashed", color = "white", size = 1.5)+
  annotate("text", label = "MaP Column 2", x = as.POSIXct("2021-12-11"), 
           y = 20, size = 3.5, colour = "black", fontface = "bold.italic") +
  
  #represent period of column 3
  geom_vline(xintercept = as.numeric(as.POSIXct("2022-01-25")), 
             linetype = "dashed", color = "white", size = 1.5)+
  geom_vline(xintercept = as.numeric(as.POSIXct("2022-04-25")), 
             linetype = "dashed", color = "white", size = 1.5)+
  annotate("text", label = "MaP Column 3", x = as.POSIXct("2022-03-17"), 
           y = 20, size = 3.5, colour = "black", fontface = "bold.italic") +
  
  # set the limit
  ylim(-600, 20) + 
  theme_bw()

# look at the frequency histogram 
bgc_benguela |>
  #pivot_longer(everything()) |>
  ggplot() +
  geom_histogram(aes(x=bbl))+
  xlim(c(0,4*10^-5))


# Try to compute bins that are different according to depth because the resolution of measurement is not the same depending on the sampling depth

bbp_binned <- bgc_benguela %>% 
  mutate(PRES = ifelse(PRES < 300, 
                               round(PRES/5)*5 + 2.5, 
                               round(PRES/20)*20 + 10)) %>% 
  group_by(CYCLE_NUMBER, PRES) %>%
  summarise(bbsr = mean(bbsr, na.rm = TRUE), 
            bbl = mean(bbl, na.rm = TRUE),
            date= mean(date, na.rm = TRUE))

# represent time series

ggplot(bbp_binned, aes(x = date, y = -PRES, color = bbl)) +
  geom_point(size = 3) +
  # set up the upper limit of the color bar to see spikes that are less intense
  scale_color_viridis(limits = c(0, 0.00012), na.value = "yellow") + 
  
  #represent period of column 1 
  geom_vline(xintercept = as.numeric(as.POSIXct("2021-10-01")), 
             linetype = "dashed", color = "white", size = 1.5)+
  geom_vline(xintercept = as.numeric(as.POSIXct("2021-10-17")), 
             linetype = "dashed", color = "white", size = 1.5)+
  annotate("text", label = "MaP Column 1", x = as.POSIXct("2021-10-10"), 
           y = 20, size = 3.5, colour = "black", fontface = "bold.italic") +
  
  #represent period of column 2
  geom_vline(xintercept = as.numeric(as.POSIXct("2021-12-01")), 
             linetype = "dashed", color = "white", size = 1.5)+
  geom_vline(xintercept = as.numeric(as.POSIXct("2021-12-19")), 
             linetype = "dashed", color = "white", size = 1.5)+
  annotate("text", label = "MaP Column 2", x = as.POSIXct("2021-12-11"), 
           y = 20, size = 3.5, colour = "black", fontface = "bold.italic") +
  
  #represent period of column 3
  geom_vline(xintercept = as.numeric(as.POSIXct("2022-01-25")), 
             linetype = "dashed", color = "white", size = 1.5)+
  geom_vline(xintercept = as.numeric(as.POSIXct("2022-04-25")), 
             linetype = "dashed", color = "white", size = 1.5)+
  annotate("text", label = "MaP Column 3", x = as.POSIXct("2022-03-17"), 
           y = 20, size = 3.5, colour = "black", fontface = "bold.italic") +
  
  # set the limit
  ylim(-600, 20) + 
  theme_bw()

write.csv(bbp_binned, "C:/Users/lemoi/Desktop/GIT/raw_bgc_data/bbp_decomposed_binned.csv", row.names = FALSE)





