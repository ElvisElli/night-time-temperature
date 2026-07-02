
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

## Set working directory ====
setwd(dirname(getActiveDocumentContext()$path))

locations <- read_excel("locations.xlsx")


### Weather data
for (i in 1:dim(locations)[1]) {
  
  weather_file_path <- file.path('outputs', paste0(locations$location[i], '.met'))
  
  pwr <- get_power_apsim_met(lonlat = unlist(locations[i, 2:3]),
                             dates =  c('1985-01-01', '2024-10-17'))
  
  pwr_imptd <- impute_apsim_met(pwr,verbose = TRUE)#using linear interpolation
  
  write_apsim_met(pwr_imptd,
                  wrt.dir = dirname(weather_file_path), 
                  filename = basename(weather_file_path))
  
}


##how to read it

marianna <- read_table("outputs/marianna.met",skip = 9,col_names = F)

write_csv(marianna,"marianna.csv")
