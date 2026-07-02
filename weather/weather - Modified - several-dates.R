#This code gives mean_max T, mean_min T, n of days with T above 32,
#n of days with T above 35, cum_rain in mm, and cumulative radiation
## Script name: ====

## Script objective: ====

## Author: Elvis F. Elli ====

## Script was created on: 2024-09-13====

## Cleaning up environment ====
rm(list=ls())

## Libraries ====
library(rstudioapi)
library(tidyverse)
library(readxl)
library(apsimx)
library(data.table)
library(nasapower)

## Set working directory ====
setwd(dirname(getActiveDocumentContext()$path))

locations <- read_excel("locations.june.july.xlsx")
all_weather_data <- list()

### Weather data
for (i in 1:nrow(locations)) {
  tryCatch({
    
    location_name <- locations$location[i]
    lonlat <- unlist(locations[i, c("longitude", "latitude")])
    start_doy <- as.numeric(locations$start_date[i])
    end_doy <- as.numeric(locations$end_date[i])
    year <- as.numeric(locations$year[i])
    
    if (is.na(year) | is.na(start_doy) | is.na(end_doy)) {
      print(paste("Skipping due to missing values for", location_name, "in year", year))
      next}
    
    ## Filter locations to keep only specific ones
    if(!(location_name %in% c("rohwer", "keiser_nere", "marianna", "stuttgart", "des_arc", "marianna_late", 
                              "rohwer_late", "stuttgart_late", "harrisburg_neerec"))) {
      print(paste("Skipping data for", location_name))  
      next}
    
    ## Get weather data for each year/location and date range
    for (y in unique(year)) {  
      pwr <- get_power_apsim_met(lonlat = lonlat,
                                 dates = c(paste0(y, "-01-01"), paste0(y, "-12-31")))  
      
      ## Filtering data within the specified date range (DOY)
      pwr_filtered <- pwr %>%
        filter(year == y & day >= start_doy & day <= end_doy) %>%
        mutate(location = location_name,
               tav = (maxt + mint)/2,             
               SVP = 0.61078 * exp((17.27 * tav)/(tav + 237.3)),  
               AVP = SVP * (rh/100),              
               VPD = SVP - AVP) %>% 
        arrange(year, location, day)
      
      ## Calculating the mean of the variables
      pwr_mean <- pwr_filtered %>%
        group_by(year, location) %>%
        summarise(mean_maxt = mean(maxt, na.rm = TRUE),
                  mean_mint = mean(mint, na.rm = TRUE),
                  cum_radn = sum(radn, na.rm = TRUE),
                  cum_rain = sum(rain, na.rm = TRUE),
                  avg_t = mean((maxt + mint) / 2, na.rm = TRUE),
                  above_32 = sum(maxt > 32, na.rm = TRUE),
                  above_35 = sum(maxt > 35, na.rm = TRUE),
                  mean_vpd = mean(VPD, na.rm = TRUE),
                  mean_rh = mean(rh))
      
      ## Adding location name as new column
      all_weather_data[[length(all_weather_data) + 1]] <- pwr_mean
      print(paste("Weather data processed for", location_name, "in year", y))}
    
  }, error = function(e){
    print(paste("Failed to process data for", locations$location[i], ":", e$message))
  })
}

# Combining all data into a single dataframe
combined_weather_data <- bind_rows(all_weather_data)

# Save the combined data as a CSV file
write_csv(combined_weather_data, "outputs/test3.csv")
#############################################################################

locations <- read_excel("locations.months.xlsx",sheet = "historical")

all_weather_data <- list()

for (i in 1:nrow(locations)) {
  tryCatch({
    
    loc <- locations$location[i]
    lonlat <- unlist(locations[i, c("longitude", "latitude")])
    y <- locations$year[i]
    start_doy <- locations$start_date[i]
    end_doy <- locations$end_date[i]
    mes <- locations$month[i]
    
    # Converting DOY to date
    start_date <- as.Date(paste0(y, "-01-01")) + start_doy - 1
    end_date   <- as.Date(paste0(y, "-01-01")) + end_doy - 1
    
    # Obtain only days of the month
    pwr <- get_power_apsim_met(lonlat = lonlat,
                               dates = c(as.character(start_date), as.character(end_date)))
    
    # Calculate variables
    pwr_summary <- pwr %>%
      mutate(tav = (maxt + mint)/2,
             SVP = 0.61078 * exp((17.27 * tav)/(tav + 237.3)),
             AVP = SVP * (rh/100),
             VPD = SVP - AVP) %>%
      summarise(location = loc,
                year = y,
                month = mes,
                min_mint = min(mint, na.rm = TRUE),
                mean_mint = mean(mint, na.rm = TRUE),
                mean_maxt = mean(maxt, na.rm = TRUE),
                avg_t = mean(tav, na.rm = TRUE),
                cum_radn = sum(radn, na.rm = TRUE),
                cum_rain = sum(rain, na.rm = TRUE),
                above_32 = sum(maxt > 32, na.rm = TRUE),
                above_35 = sum(maxt > 35, na.rm = TRUE),
                mean_vpd = mean(VPD, na.rm = TRUE),
                mean_rh = mean(rh, na.rm = TRUE))
    
    all_weather_data[[length(all_weather_data) + 1]] <- pwr_summary
    
    message("Procesado: ", loc, " ", y, " ", mes)
    
  }, error = function(e) {
    message("Error en ", loc, " año ", y, ": ", e$message)
  })
}

combined_weather_data <- bind_rows(all_weather_data)
write_csv(combined_weather_data, "outputs/historical_weather.csv")

########################################################################
library(readxl)
library(dplyr)
library(purrr)
library(apsimx)
library(ggplot2)
library(patchwork)

loc <- read_excel("locations.2024.2025.xlsx") 
  

# --- 1. Pull 2024-2025 data ---
all_weather <- map_dfr(1:nrow(loc), function(i) {
  
  lonlat     <- unlist(loc[i, c("longitude", "latitude")])
  y          <- loc$year[i]
  location   <- loc$location[i]
  start_date <- as.Date(paste0(y, "-01-01")) + loc$start_date[i] - 1
  end_date   <- as.Date(paste0(y, "-01-01")) + loc$end_date[i] - 1
  
  pwr <- get_power_apsim_met(
    lonlat = lonlat,
    dates  = c(as.character(start_date), as.character(end_date)))
  
  pwr %>%
    mutate(
      date      = as.Date(paste(year, day, sep = "-"), format = "%Y-%j"),
      mean_temp = (maxt + mint) / 2,
      month     = as.integer(format(date, "%m")),
      location  = location,
      year      = y) %>%
    group_by(location, year, month) %>%
    summarise(
      mean_temp  = mean(mean_temp, na.rm = TRUE),
      mean_maxt  = mean(maxt,      na.rm = TRUE),
      mean_mint  = mean(mint,      na.rm = TRUE),
      total_rain = sum(rain,       na.rm = TRUE),
      mean_radn  = mean(radn,      na.rm = TRUE),
      min_radn   = min(radn,       na.rm = TRUE),
      max_radn   = max(radn,       na.rm = TRUE),
      mean_rh    = mean(rh,        na.rm = TRUE),
      .groups = "drop")
})

write_csv(all_weather, "weather_monthly_summary.csv")


# --- 2. Pull historical data (1995-2025) ---
loc_unique <- loc %>% distinct(location, longitude, latitude)

historical <- map_dfr(1:nrow(loc_unique), function(i) {
  
  lonlat   <- unlist(loc_unique[i, c("longitude", "latitude")])
  location <- loc_unique$location[i]
  
  pwr <- get_power_apsim_met(
    lonlat = lonlat,
    dates  = c("1995-01-01", "2025-12-31"))
  
  pwr %>%
    mutate(
      date      = as.Date(paste(year, day, sep = "-"), format = "%Y-%j"),
      mean_temp = (maxt + mint) / 2,
      month     = as.integer(format(date, "%m")),
      location  = location) %>%
    filter(month %in% 3:10) %>%
    group_by(location, month) %>%
    summarise(
      mean_temp  = mean(mean_temp, na.rm = TRUE),
      mean_maxt  = mean(maxt,      na.rm = TRUE),
      mean_mint  = mean(mint,      na.rm = TRUE),
      total_rain = sum(rain,       na.rm = TRUE) / n_distinct(year),
      mean_radn  = mean(radn,      na.rm = TRUE),
      .groups = "drop")
}) %>%
  mutate(month = factor(month, levels = 3:10,
                        labels = c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct")),
         location = recode(location,
                           "fayetteville_sarec" = "Fayetteville",
                           "keiser_nere"        = "Keiser",
                           "marianna"           = "Marianna",
                           "rohwer"             = "Rohwer",
                           "stuttgart"          = "Stuttgart"))

# --- 3. Prep plot data ---
year_colors <- c("2024" = "#3A7BD5", "2025" = "#E8A838", "Historical" = "black")

plot_data <- all_weather %>%
  filter(month %in% 3:10) %>%
  mutate(year = factor(year),
         month = factor(month, levels = 3:10,
                        labels = c("Mar","Apr","May","Jun","Jul","Aug","Sep","Oct")),
         location = recode(location,
                           "fayetteville_sarec" = "Fayetteville",
                           "keiser_nere"        = "Keiser",
                           "marianna"           = "Marianna",
                           "rohwer"             = "Rohwer",
                           "stuttgart"          = "Stuttgart"))

hist_rain <- historical %>%
  mutate(year = factor("Historical"))

plot_data_rain <- plot_data %>%
  bind_rows(hist_rain %>% select(location, month, total_rain, year))

# --- 4. Plots ---
p_rain <- ggplot(plot_data_rain, aes(x = month, y = total_rain, fill = year)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  facet_wrap(~location, nrow = 1) +
  scale_fill_manual(values = year_colors) +
  labs(y = "Precipitation (mm)", x = NULL, fill = NULL) +
  theme_bw() +
  theme(strip.text = element_text(size = 9),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top")

p_radn <- ggplot(plot_data, aes(x = month, y = mean_radn, color = year, group = year)) +
  # geom_ribbon(aes(ymin = min_radn, ymax = max_radn, fill = year),
  #             alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_line(data = historical, aes(x = month, y = mean_radn, group = 1, color = "Historical"),
            inherit.aes = FALSE, linetype = "dashed", linewidth = 0.7) +
  geom_point(data = historical, aes(x = month, y = mean_radn, color = "Historical"),
             inherit.aes = FALSE, shape = 21, fill = "white", size = 2) +
  facet_wrap(~location, nrow = 1) +
  scale_color_manual(values = year_colors) +
  scale_fill_manual(values = year_colors) +
  labs(y = expression("Solar radiation (MJ m"^-2~"day"^-1~")"), x = NULL, color = NULL, fill = NULL) +
  theme_bw() +
  theme(strip.text = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

p_temp <- ggplot(plot_data, aes(x = month, group = year, color = year)) +
  # geom_ribbon(aes(ymin = mean_mint, ymax = mean_maxt, fill = year),
  #             alpha = 0.15, color = NA) +
  geom_line(aes(y = mean_temp), linewidth = 0.8) +
  geom_point(aes(y = mean_temp), size = 2) +
  geom_line(data = historical, aes(x = month, y = mean_temp, group = 1, color = "Historical"),
            inherit.aes = FALSE, linetype = "dashed", linewidth = 0.7) +
  facet_wrap(~location, nrow = 1) +
  scale_color_manual(values = year_colors) +
  scale_fill_manual(values = year_colors) +
  labs(y = "Temperature (°C)", x = "Month", color = NULL, fill = NULL) +
  theme_bw() +
  theme(strip.text = element_blank(),
        legend.position = "none")

# --- 5. Save ---
ggsave(
  plot = (p_rain / p_radn / p_temp) + plot_layout(heights = c(1, 1, 1)),
  "environmental_characterization.tiff",
  width = 12, height = 8, units = "in", dpi = 600, bg = "white")
######################################################################################################

#This code gives as an output mean_t and rain for several months for several selected months

### weather per month
## Script name: ====

## Script objective: Monthly weather summaries from NASA POWER data for specific locations

## Author: Elvis F. Elli ====

## Script was created on: 2024-09-13====

## Cleaning up environment ====
rm(list = ls())

## Libraries ====
library(rstudioapi)
library(tidyverse)
library(readxl)
library(apsimx)
library(data.table)
library(nasapower)
library(lubridate)

## Set working directory ====
setwd(dirname(getActiveDocumentContext()$path))

## Read location info ====
locations <- read_excel("locations.3.xlsx")
all_weather_data <- list()

### Weather data processing ====
for (i in 1:nrow(locations)) {
  tryCatch({
    
    location_name <- locations$location[i]
    lonlat <- unlist(locations[i, c("longitude", "latitude")])
    start_doy <- as.numeric(locations$start_date[i])
    end_doy <- as.numeric(locations$end_date[i])
    year <- as.numeric(locations$year[i])
    
    if (is.na(year) | is.na(start_doy) | is.na(end_doy)) {
      print(paste("Skipping due to missing values for", location_name, "in year", year))
      next
    }
    
    ## Filter locations to keep only specific ones
    if (!(location_name %in% c("rohwer", "keiser_nere", "marianna", "stuttgart", "bell_farm", 
                               "marianna_late", "rohwer_late", "stuttgart_late", "harrisburg_neerec"))) {
      print(paste("Skipping data for", location_name))
      next
    }
    
    ## Get weather data for the specified year
    pwr <- get_power_apsim_met(
      lonlat = lonlat,
      dates = c(paste0(year, "-01-01"), paste0(year, "-12-31"))
    )
    
    ## Filter to planting-harvest window
    pwr_filtered <- pwr %>%
      filter(year == year & day >= start_doy & day <= end_doy) %>%
      mutate(
        location = location_name,
        date = as.Date(paste0(year, "-01-01")) + (day - 1),
        month = lubridate::month(date)
      ) %>%
      arrange(year, location, day)
    
    ## Monthly summaries
    pwr_monthly <- pwr_filtered %>%
      group_by(year, month, location) %>%
      summarise(
        mean_tavg = mean((maxt + mint) / 2, na.rm = TRUE),
        cum_rain = sum(rain, na.rm = TRUE),
        .groups = "drop"
      )
    
    ## Store
    all_weather_data[[length(all_weather_data) + 1]] <- pwr_monthly
    print(paste("Monthly weather summary completed for", location_name, "in year", year))
    
  }, error = function(e) {
    print(paste("Failed to process data for", locations$location[i], ":", e$message))
  })
}

## Combine all data
combined_weather_data <- bind_rows(all_weather_data)

## Save to CSV
write_csv(combined_weather_data, "outputs/test1.csv")
