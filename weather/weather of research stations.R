#load packages
library(tidyverse)
library(readxl)
library(lubridate)

## Read location and sowing information ====
locations <- read_excel("locations.3.xlsx")

# Filter for rohwer only
rohwer_info <- locations %>%
  filter(location == "rohwer") %>%
  mutate(start_doy = as.numeric(start_date),
         end_doy = as.numeric(end_date),
         year = as.numeric(year)) %>%
  drop_na(year, start_doy, end_doy)

## Read weather data from research station
weather_raw <- read_delim("rohwer.weather.txt", delim = "\t", show_col_types = FALSE)

# Clean and prepare weather data
weather_clean <- weather_raw %>%
  rename(
    date = Date,
    maxt = `Air Temp Max(ºC)`,
    mint = `Air Temp Min(ºC)`,
    rain = `Rainfall (cm)`
  ) %>%
  mutate(
    date = mdy(date),
    year = year(date),
    doy = yday(date),
    tavg = (maxt + mint) / 2,
    above_32 = if_else(maxt > 32, 1, 0)
  ) %>%
  filter(!is.na(year) & !is.na(doy))

## Process for each year
all_weather_data <- list()

for (i in 1:nrow(rohwer_info)) {
  y <- rohwer_info$year[i]
  s_doy <- rohwer_info$start_doy[i]
  e_doy <- rohwer_info$end_doy[i]
  
  filtered <- weather_clean %>%
    filter(year == y, doy >= s_doy, doy <= e_doy)
  
  summary_stats <- filtered %>%
    summarise(
      year = y,
      location = "rohwer",
      mean_maxt = mean(maxt, na.rm = TRUE),
      mean_mint = mean(mint, na.rm = TRUE),
      cum_rain = sum(rain, na.rm = TRUE),
      avg_t = mean(tavg, na.rm = TRUE),
      above_32 = sum(above_32, na.rm = TRUE)
    )
  
  all_weather_data[[i]] <- summary_stats
  print(paste("Processed year", y, "for rohwer"))
}

# Combine all and write to CSV
combined_weather_data <- bind_rows(all_weather_data)
write_csv(combined_weather_data, "outputs/rohwer_station_weather_summary.csv")
