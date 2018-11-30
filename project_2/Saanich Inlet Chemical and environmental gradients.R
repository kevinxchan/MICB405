raw_dat <- read_csv(file = "Saanich_Data.csv", col_names = T)

Chemical_gradients <-  raw_dat %>%
  dplyr::select(Cruise, Date, Depth, WS_O2, WS_PO4, SI, WS_NO3, Mean_NH4, Mean_NO2,
                WS_H2S, Cells.ml, Mean_N2, Mean_O2, Mean_co2, Mean_N2O, Mean_CH4) %>%
  mutate(Depth_m = Depth * 1000) %>%
  gather(key = "Chemical", value = "Concentrations",
         -Cruise, -Date, -Depth_m, -Depth, -Cells.ml) %>%
  mutate(Chemical = gsub("Mean_co2", "Mean_CO2", Chemical),
         Chemical = gsub("_", " ", Chemical))

ggplot(Chemical_gradients, aes(y = Depth_m, x = Concentrations, colour = Chemical)) +
  geom_point(alpha = 0.25) +
  geom_line(aes(group = Cruise), alpha = 0.25) +
  facet_wrap(~ Chemical, scales = "free_x") +
  scale_y_reverse() +
  theme_bw() +
  scale_colour_discrete(guide = F) +
  labs(y = "Depth (meters)", x = expression("Concentration ("*mu*"M)")) +
  theme(
    text = element_text(size = 14),
    panel.background = element_blank(),
    panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
    panel.grid.minor = element_blank()
  )

Chemical_gradients %>%
  filter(Cruise == 42 |
           Cruise == 48 | 
           Cruise == 72 | 
           Cruise == 73 | 
           Cruise == 74 | 
           Cruise == 75) %>%
  filter(Chemical != "Mean CO2",
         Chemical != "Mean N2",
         Chemical != "Mean O2") %>%
  mutate(Cruise = gsub("42", "SI042", Cruise),
         Cruise = gsub("48", "SI048", Cruise),
         Cruise = gsub("72", "SI072", Cruise),
         Cruise = gsub("73", "SI073", Cruise),
         Cruise = gsub("74", "SI074", Cruise),
         Cruise = gsub("75", "SI075", Cruise)) %>%
  ggplot(aes(y = Depth_m, x = Concentrations)) +
  geom_point(alpha = 0.75, size = 2, aes(pch = factor(Cruise), colour = factor(Date))) +
  geom_line(aes(group = Cruise, colour = factor(Date)), alpha = 0.5) +
  facet_wrap(~ Chemical, scales = "free_x") +
  scale_y_reverse() +
  scale_shape_manual(values = c(16, 15, 17, 18, 1, 6), name = "Cruise") +
  theme_bw() +
  scale_colour_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
                                   "#E6AB02"),
                        name = "Date of sampling") +
  labs(y = "Depth (meters)", x = expression("Concentration ("*mu*"M)")) +
  theme(
    text = element_text(size = 14),
    panel.background = element_blank(),
    panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
    panel.grid.minor = element_blank()
  )

# "Note: different scales; lines represent measurements taken from the same cruise"

# Seeing if time may be a possible factor for chemical concentration gradients
## Ideally, you'd expect some sort of cycling of nutrients by seasons (maybe at some of the depths)
ggplot(Chemical_gradients, aes(x = Date, y = Concentrations, colour = Chemical)) +
  geom_point(alpha = 0.25) +
  geom_line(aes(group = Cruise), alpha = 0.25) +
  facet_wrap(~ Chemical, scales = "free_y") +
  theme_bw() +
  stat_smooth(se = F, method = "loess", colour = "black", size = 0.5) 

# quick correlation to see what's up (in terms of correlations)
library(PerformanceAnalytics)

Chemical_gradients_corr <- Chemical_gradients %>%
  spread(key = Chemical, value = Concentrations) %>%
  dplyr::select(-Depth, -Cells.ml) 

chart.Correlation(data.matrix(Chemical_gradients_corr)) 

# Looking at other aspects of the samples as a function of depth

depth_characteristics <-  raw_dat %>%
  dplyr::select(Cruise, Date, Depth, Longitude, Latitude, Cells.ml, Temperature,
                Salinity, Density) %>%
  mutate(Depth_m = Depth * 1000,
         log_Cells.ml = log10(Cells.ml)) %>%
  dplyr::select(-Depth, -Cells.ml) %>%
  gather(key = "Characteristics", value = "Measures",
         -Cruise, -Date, -Depth_m, -Longitude, -Latitude) %>%
  mutate(Characteristics = gsub("log_Cells.ml", "Log Cells per ml", Characteristics),
         Measures = gsub("-Inf", NA, Measures))

depth_characteristics$Measures <- as.numeric(depth_characteristics$Measures)

## Cruise 93 is -Inf, so must filter
tibble <- depth_characteristics %>%
  filter(Characteristics == "Log Cells per ml") %>%
  group_by(Cruise) %>%
  filter(!is.na(Measures)) %>%
  summarize(mean = mean(Measures),
            sum = sum(Measures))

Cruise_93 <- depth_characteristics %>%
  filter(Cruise == 93)

depth_characteristics %>%
  filter(Cruise == 42 |
           Cruise == 48 | 
           Cruise == 72 | 
           Cruise == 73 | 
           Cruise == 74 | 
           Cruise == 75) %>%
  mutate(Cruise = gsub("42", "SI042", Cruise),
         Cruise = gsub("48", "SI048", Cruise),
         Cruise = gsub("72", "SI072", Cruise),
         Cruise = gsub("73", "SI073", Cruise),
         Cruise = gsub("74", "SI074", Cruise),
         Cruise = gsub("75", "SI075", Cruise)) %>%
ggplot(aes(y = Depth_m, x = Measures, colour = Characteristics)) +
  geom_point(alpha = 0.75, size = 2, aes(pch = factor(Cruise), colour = factor(Date))) +
  geom_line(aes(group = Cruise, colour = factor(Date)), alpha = 0.5) +
  facet_wrap(~ Characteristics, scales = "free_x") +
  scale_y_reverse() +
  scale_shape_manual(values = c(16, 15, 17, 18, 1, 6), name = "Cruise") +
  theme_bw() +
  scale_colour_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
                                 "#E6AB02"),
                      name = "Date of sampling") +
  labs(subtitle = expression("Units: Density ("*theta*"); Cells (per ml); Salinity (psu); Temperature (Celsius)"),
       y = "Depth (meters)", x = "Evironmental parameters") +
  theme(
    text = element_text(size = 14),
    panel.background = element_blank(),
    panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
    panel.grid.minor = element_blank()
  )
  
  
  
  geom_point(alpha = 0.25) +
  geom_line(aes(group = Cruise), alpha = 0.25) +
  facet_wrap(~ Characteristics, scales = "free_x") +
  scale_y_reverse() +
  theme_bw() +
  scale_colour_discrete(guide = F) +
  labs(subtitle = expression("Units: Density ("*theta*"); Cells (per ml); Salinity (psu); Temperature (Celsius)"),
       y = "Depth (meters)", x = "Evironmental parameters") +
  theme(
        text = element_text(size = 14),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank()
  )



  
  
  
  
  
  
  
  
  
  raw_dat %>%
  dplyr::select(Cruise, Date, Depth, Temperature,
                WS_O2, WS_NO3, WS_H2S, SI, WS_PO4, Salinity, Density) %>%
  dplyr::rename(O2_uM=WS_O2, NO3_uM=WS_NO3, H2S_uM=WS_H2S, PO4_uM = WS_PO4) %>%
  mutate(Depth_m=Depth*1000) %>%
  dplyr::select(-Depth) %>%
  gather(key = "Chemical", value = "Concentration", -Cruise, -Date, -Temperature, -SI,
         -Salinity, -Density, -Depth_m)

dat_project_2 %>%
  dplyr::select(Depth_m, Temperature, Salinity, Density) %>%
  melt(id = c("Depth_m")) %>%
  ggplot(aes(x = Depth_m, y = value, colour = variable)) +
  geom_point(alpha = 0.10) +
  facet_wrap(~ variable, scales = "free_y") +
  theme_bw() +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", 
               size = 0.75, width = 10, alpha = 0.75) +
  stat_smooth(method = "loess", se = T) +
  labs(x = "Depth (m)", y = "Measure", subtitle = "Temperature (Celsius); Salinity ; Density ")


dat_class %>%
  melt(id = c("Cruise", "Date", "Temperature", "Depth_m", "Depth")) %>%
  ggplot(aes(x = Depth_m, y = value, colour = variable)) +
  geom_point(alpha = 0.15) +
  facet_wrap(~ variable, scales = "free_y") +
  theme_bw() +
  stat_smooth(method = "loess", se = T) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar", 
               size = 0.75, width = 10, alpha = 0.75) +
  geom_abline(slope = 0, intercept = 0, lty = 3, alpha = 0.75) +
  labs(x = "Depth (m)", y = "Concentration (uM)")

Chemical_gradients %>%
  filter(Cruise == 42 |
           Cruise == 48 | 
           Cruise == 72 | 
           Cruise == 73 | 
           Cruise == 74 | 
           Cruise == 75) %>%
  ggplot(aes(x = Date, y = Concentration, colour = Chemical)) +
  geom_point() +
  facet_grid(~ Cruise, scales = "free_y")

