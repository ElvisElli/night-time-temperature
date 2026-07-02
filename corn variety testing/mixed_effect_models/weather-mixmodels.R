#loading libraries: 
library(tidyverse)
library(dplyr)
library(nlme)
library(MuMIn)
library(dplyr)
library(nlraa)
library(car)
library(corrplot)

#loading data:

df <- read.csv("input/Corn data organized2.csv", stringsAsFactors = FALSE)
data_model <- df %>%
  select(yieldkg.ha, year, location, hybrid, mean_mint, above_35, cum_radn.j,vpd, 
         april, may, june, july,august, sept, precipitation) %>% 
  mutate(yieldkg.ha = as.numeric(yieldkg.ha),
         above_35 = as.numeric(above_35))%>%
  drop_na() 

data_model_z <- data_model %>%
  mutate(
    z_mean_mint = as.numeric(scale(mean_mint)),
    z_above_35  = as.numeric(scale(above_35)),
    z_precipitation  = as.numeric(scale(precipitation)),
    z_cum_radn  = as.numeric(scale(cum_radn.j)),
    z_vpd       = as.numeric(scale(vpd)),
    z_april     = as.numeric(scale(april)),
    z_may       = as.numeric(scale(may)),
    z_june      = as.numeric(scale(june)),
    z_july      = as.numeric(scale(july)),
    z_august    = as.numeric(scale(august)),
    z_sept      = as.numeric(scale(sept)))

vars <- data_model_z %>% 
  select(z_mean_mint, z_above_35, z_cum_radn, z_vpd,
         z_april, z_may, z_july, z_august, z_sept, year)

cor_mat <- cor(vars, use = "pairwise.complete.obs")
cor_mat

#Including year as a random effect separates variation between years from variation within a year. 
#The fixed effect of radiation then represents only the within-year variation, which can be negative if, within the same year, 
#locations with more radiation do not consistently yield more. This does not mean radiation reduces yield overall; 
#it simply reflects the residual effect after accounting for year-to-year differences.

m_climate_z <- lme(
  fixed  = yieldkg.ha ~ z_mean_mint + z_cum_radn + z_vpd + z_april +
    z_may  + z_july + z_august + z_sept,
  random = ~ 1 | location/hybrid,
  weights = varIdent(form = ~1 | location),  
  data   = data_model_z)
summary(m_climate_z)
VarCorr(m_climate_z)

vc_mat <- VarCorr(m_climate_z)

# Checking for VIF (variance inflation factor)
fixed_vars <- data_model_z %>%
  select(yieldkg.ha,z_mean_mint, z_cum_radn, z_vpd, z_april, z_may, z_july, z_august, z_sept)

# Ajustar un lm temporal
lm_fixed <- lm(yieldkg.ha ~ ., data = fixed_vars)

# Calcular VIF
vif(lm_fixed)

var_loc <- as.numeric(vc_mat[2, "Variance"])
var_hyb <- as.numeric(vc_mat[4, "Variance"])
var_res <- as.numeric(vc_mat[5, "Variance"])

var_total <- var_loc + var_hyb + var_res

variance_partition <- data.frame(
  component = c("Location", "Hybrid", "Residual"),
  variance  = c(var_loc, var_hyb, var_res),
  percent   = c(var_loc, var_hyb, var_res) / var_total * 100)

variance_partition
### 

#plot residuals

qqnorm(residuals(m_climate_z))
qqline(residuals(m_climate_z), col = "red", lwd = 2)

# rank variables
coef_rank_z <- summary(m_climate_z)$tTable %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  filter(variable != "(Intercept)") %>%
  mutate(abs_effect = abs(Value)) %>%
  arrange(desc(abs_effect)) %>%
  select(variable, Value, Std.Error, `t-value`, `p-value`, abs_effect)

print(coef_rank_z)

#Rank variables:

coef_rank_z <- coef_rank_z %>%
  mutate(variable = str_replace(variable, "z_mean_mint", "June-July Min Temperature"),
         variable = str_replace(variable, "z_above_35", "Days T > 35°C"),
         variable = str_replace(variable, "z_cum_radn", "Cumulative radiation"),
         variable = str_replace(variable, "z_april", "Apr rain (mm)"),
         variable = str_replace(variable, "z_may", "May rain (mm)"),
         variable = str_replace(variable, "z_june", "Jun rain (mm)"),
         variable = str_replace(variable, "z_july", "July rain (mm)"),
         variable = str_replace(variable, "z_august", "Aug rain (mm)"),
         variable = str_replace(variable, "z_sept", "Sep rain (mm)"),
         variable = str_replace (variable, "z_vpd", "Jun-July VPD (kPa)"))

# final plot

weather2 <- ggplot(coef_rank_z, aes(x = reorder(variable, Value), y = Value, fill = Value)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.0f", Value)), hjust = -0.2, size = 4.5, family = "sans") +
  scale_fill_gradient(low = "tomato", high = "steelblue") +
  coord_flip() +
  labs(
    x = NULL,
    y = expression("Standardized effect size (|"*beta*"|)")) +
  theme_classic(base_size = 16) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "gray40"),
    panel.grid.major.x = element_line(color = "gray85", linetype = "dotted"),
    legend.position = "none",
    plot.margin = margin(10, 30, 10, 10)) +
  expand_limits(y = max(coef_rank_z$abs_effect) * 1.15)
ggsave(plot = weather2, "output/Figure1.tiff",
       width = 8, height = 7, unit = "in", dpi = 600, bg = "white")

##################################################################

sd_mean_mint <- sd(data_model$mean_mint, na.rm = TRUE)
sd_mean_mint

data_model %>% 
  ggplot(aes(x= mean_mint, y= yieldkg.ha))+
  geom_point()+ 
  geom_smooth(method = "loess")


data_mean <- data_model %>% 
  group_by(mean_mint) %>% 
  summarize(mean_yieldkg.ha = mean(yieldkg.ha), .groups = "drop")

mod2 <- lm(mean_yieldkg.ha ~ mean_mint, data = data_mean)

# Extract slope, intercept, p-value
slope_value <- coef(mod2)[2]
intercept_value <- coef(mod2)[1]
p_value <- summary(mod2)$coefficients["mean_mint", "Pr(>|t|)"]

# Create label
label <- paste0(
  "Slope = ", round(slope_value), " kg/ha/°C\n",
  "p = ", formatC(p_value, format = "e", digits = 2))

p2 <- data_mean %>%
  mutate(Y_pred = intercept_value + slope_value * mean_mint) %>%
  ggplot(aes(x = mean_mint, y = mean_yieldkg.ha)) +
  geom_line(aes(y = Y_pred), linewidth = 1) +
  geom_point(fill = "green4", shape = 21, size = 4, alpha = 0.8) +
  scale_y_continuous(breaks = seq(10000, 18000, by = 2000),
    limits = c(10000, 18000)) +
  labs(
    x = "Night-time temperature (°C)",
    y = expression("Mean yield (kg ha"^{-1}*")")) +
  geom_text(x = 23.5, 
            y = 17000, 
            label = label, size = 6, fontface = "italic") +
  theme_classic(base_size = 25)
ggsave(plot = p2, "output/night-time.tiff", width = 18, height = 14 , unit = "cm", dpi = 600, bg= "white")

############### increasing night temperature

df <- read.csv("input/Corn data organized2.csv", stringsAsFactors = FALSE)
data_model <- df %>%
  select(yield, yieldkg.ha, year, location, hybrid, mean_mint, above_35, cum_radn.j,vpd, 
         april, may, june, july,august, sept, precipitation) %>% 
  mutate(yield= as.numeric(yield),
         yieldkg.ha= as.numeric(yieldkg.ha),
         above_35 = as.numeric(above_35))%>%
  drop_na() 

## scaling variables
data_model_z <- data_model %>%
  mutate(
    z_mean_mint = as.numeric(scale(mean_mint)),
    z_above_35  = as.numeric(scale(above_35)),
    z_precipitation  = as.numeric(scale(precipitation)),
    z_cum_radn  = as.numeric(scale(cum_radn.j)),
    z_vpd       = as.numeric(scale(vpd)),
    z_april     = as.numeric(scale(april)),
    z_may       = as.numeric(scale(may)),
    z_june      = as.numeric(scale(june)),
    z_july      = as.numeric(scale(july)),
    z_august    = as.numeric(scale(august)),
    z_sept      = as.numeric(scale(sept)))

m_climate_z <- lme(
  fixed  = yieldkg.ha ~ z_mean_mint + z_cum_radn + z_vpd + z_april +
    z_may  + z_july + z_august + z_sept,
  random = ~ 1 | location/hybrid,
  weights = varIdent(form = ~1 | location),  
  data   = data_model_z)
summary(m_climate_z)

qqnorm(residuals(m_climate_z))
qqline(residuals(m_climate_z), col = "red")

## add scenarios of increasing T
data_model_z <- data_model_z %>%
  mutate(
    plus1 = mean_mint + 1,
    plus2 = mean_mint + 2,
    plus3 = mean_mint + 3,
    plus4 = mean_mint + 4)

## pivot longer for data analysis
data_long <- data_model_z %>%
  pivot_longer(
    cols = c(mean_mint, plus1, plus2, plus3, plus4),
    names_to = "scenario",
    values_to = "mean_mint_temp") %>%
  mutate(
    z_mean_mint = (mean_mint_temp - mean(data_model$mean_mint, na.rm = TRUE)) /
      sd(data_model$mean_mint, na.rm = TRUE))

# --- NUEVO: calcular baseline para escenario 0 ---
baseline_yields <- data_long %>%
  filter(scenario == "mean_mint") %>%
  mutate(original_yield = predict(m_climate_z, newdata = cur_data())) %>%
  select(year, location, hybrid, original_yield)

# --- Predecir rendimientos para todos los escenarios ---
data_long <- data_long %>%
  left_join(baseline_yields, by = c("year", "location", "hybrid")) %>%
  rowwise() %>%
  mutate(
    predicted_yield = if_else(
      scenario == "mean_mint",
      original_yield,  
      predict(m_climate_z, newdata = cur_data())),
    pct_change = 100 * (predicted_yield - original_yield) / original_yield) %>%
  ungroup()

# --- resumen por escenario ---
yield_change_by_scenario <- data_long %>%
  group_by(scenario) %>%
  summarise(
    mean_yield_change = mean(predicted_yield - original_yield, na.rm = TRUE),
    mean_pct_change   = mean(pct_change, na.rm = TRUE))

yield_change_by_scenario

#### graphs
p3 <- data_long %>%
  mutate(
    temp_increase = case_when(
      scenario == "mean_mint" ~ 0,
      scenario == "plus1" ~ 1,
      scenario == "plus2" ~ 2,
      scenario == "plus3" ~ 3,
      scenario == "plus4" ~ 4)) %>%
  ggplot(aes(x = temp_increase, y = pct_change)) +
  geom_line(aes(group = interaction(location, year)), alpha = 0.3) +
  geom_point(size = 3, alpha = 0.6, color = "green4") +
  geom_smooth(se = FALSE, linewidth = 1.2, color = "black") +
  labs(
    x = "Night-time temperature increase (°C)",
    y = "Yield change (%)") +
  theme_classic(base_size = 25)

ggsave(plot = p3, "output/NT.predictedyield.tiff", width = 18, height = 14 , unit = "cm", dpi = 600, bg= "white")

data_long %>%
  mutate(
    temp_increase = case_when(
      scenario == "mean_mint" ~ 0,
      scenario == "plus1" ~ 1,
      scenario == "plus2" ~ 2,
      scenario == "plus3" ~ 3,
      scenario == "plus4" ~ 4)) %>%
  ggplot(aes(x = temp_increase, y = pct_change)) +
  # Líneas individuales por Location-Year y escenario
  geom_line(aes(group = interaction(location, year, scenario)), alpha = 0.3) +
  # Individual points
  geom_point(size = 3, alpha = 0.6, color = "green4") +
  # T Overall slope by scenario
  geom_smooth(aes(group = scenario), se = FALSE, linewidth = 1.2, color = "black") +
  labs(
    x = "Night-time temperature increase (°C)",
    y = "Yield change (%)") +
  theme_classic(base_size = 20)


p4 <- data_long %>% 
  mutate(temp_increase = case_when(
    scenario == "mean_mint" ~ 0,
    scenario == "plus1"      ~ 1,
    scenario == "plus2"      ~ 2,
    scenario == "plus3"      ~ 3,
    scenario == "plus4"      ~ 4)) %>%
  group_by(temp_increase) %>%
  summarise(mean_pct = mean(pct_change, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(temp_increase %in% 1:4) %>%
  ggplot(aes(x = factor(temp_increase), y = mean_pct)) +  
  geom_col(fill = "#CC6666", alpha = 0.8) +
  scale_y_continuous(breaks = seq(-15,0, by=5), limits = c(-16,0)) +
  labs(
    x = "Night-time temperature increase (°C)",
    y = "Mean yield change (%)") +
  theme_minimal(base_size = 25)
ggsave(plot = p4, "output/NTyield.paper1.tiff", width = 8, height = 7 , unit = "in", dpi = 600, bg= "white")
#####
#mean of observed data
baseline_yield <- original_yields %>%
  summarise(mean_baseline = mean(original_yield, na.rm = TRUE))

#mean of predicted data in mean_mint scenario
# Average predicted yield for the baseline scenario (mean_mint)
mean_predicted_baseline <- data_long %>%
  filter(scenario == "mean_mint") %>%
  summarise(avg_predicted = mean(predicted_yield, na.rm = TRUE))

mean_predicted_baseline
### RMSE
obs <- data_model_z$yield
pred <- fitted(m_climate_z)

rmse <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
rmse

#### RMSE CV
cv_rmse <- data_model_z %>%
  group_split(location) %>%
  lapply(function(test) {
    
    train <- anti_join(data_model_z, test, by = c("year","location","hybrid"))
    
    m <- lme(
      yield ~ z_mean_mint + z_cum_radn + z_vpd +
        z_april + z_may + z_july + z_august + z_sept,
      random = ~1 | location/hybrid,
      weights = varIdent(form = ~1 | location),
      data = train)
    
    data.frame(
      obs  = test$yield,
      pred = predict(m, test, level = 0))
  }) %>%
  bind_rows() %>%
  summarise(rmse = sqrt(mean((obs - pred)^2))) %>%
  pull(rmse)

cv_rmse

### June july T

df <- read.csv("input/Corn data organized2.csv", stringsAsFactors = FALSE)
data_model <- df %>%
  select(yield, yieldkg.ha, year, location, hybrid, mean_mint, above_35, cum_radn.j,vpd, 
         april, may, june, july,august, sept, precipitation) %>% 
  mutate(yield= as.numeric(yield),
         yieldkg.ha= as.numeric(yieldkg.ha),
         above_35 = as.numeric(above_35))%>%
  drop_na() 

ggplot(df, aes(x = year, y = mean_mint)) +
  geom_line(aes(color = location),alpha = 0.7, linewidth = 0.8) +
  geom_point(aes(color = location), size = 2) +
  scale_x_continuous(breaks = seq(2010,2024,by=2))+
  scale_y_continuous(breaks = seq(20,25,by=1), limits = c(20,25))+
  labs(
    x = "Year",
    y = "Night-time temperature (°C)",
    title = "Observed Night-time Temperature by Year and Location") +
  theme_classic(base_size = 16)

ggplot(df, aes(x = year, y = avg_t)) +
  geom_line(aes(color = location),alpha = 0.7, linewidth = 0.8) +
  geom_point(aes(color = location), size = 2) +
  scale_x_continuous(breaks = seq(2010,2024,by=2))+
  scale_y_continuous(breaks = seq(24,32,by=1), limits = c(24,32))+
  labs(
    x = "Year",
    y = "Average June-July T (°C)",
    title = "Average T") +
  theme_classic(base_size = 16)


ggplot(df, aes(x = year, y = mean_maxt)) +
  geom_line(aes(color = location),alpha = 0.7, linewidth = 0.8) +
  geom_point(aes(color = location), size = 2) +
  scale_x_continuous(breaks = seq(2010,2024,by=2))+
  scale_y_continuous(breaks = seq(30,38,by=1), limits = c(30,38))+
  labs(
    x = "Year",
    y = "Maximum June-July T (°C)",
    title = "Maximum T") +
  theme_classic(base_size = 16)


df_long <- df %>%
  pivot_longer(cols = c(mean_maxt, mean_mint, avg_t), names_to = "variable", values_to = "temperature") %>%
  mutate(variable = dplyr::recode(variable, "mean_maxt" = "June–July Maximum T", 
                                  "mean_mint" = "June–July Minimum T", 
                                  "avg_t" = "June–July Mean T"))

temperatures <- ggplot(df_long, aes(x = year, y = temperature, color = variable, shape = variable)) +
  geom_line(aes(group = interaction(location, variable)), alpha = 0.7, linewidth = 1) +
  geom_point(size = 2) +
  geom_smooth(aes(group = variable), method = "lm", se = FALSE, linewidth = 1.2) +
  scale_x_continuous(breaks = seq(2010, 2024, by = 2)) +
  scale_y_continuous(breaks = seq(20, 38, by = 2)) +
  labs(x = "Year", y = "Temperature (°C)", color = NULL, shape = NULL) +
  theme_classic(base_size = 25)

ggsave(plot = temperatures, "output/temperatures.paper1.tiff", width = 8, height = 7 , unit = "in", dpi = 600, bg= "white")

df_avg <- df_long %>%
  group_by(year, variable) %>%
  summarise(temperature = mean(temperature, na.rm = TRUE), .groups = "drop")

# Gráfico
temperatures.avg <- ggplot(df_avg, aes(x = year, y = temperature, color = variable, shape = variable)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_smooth(aes(group = variable), method = "lm", se = FALSE, linewidth = 1.2) +
  scale_x_continuous(breaks = seq(2010, 2024, by = 2)) +
  scale_y_continuous(breaks = seq(20, 38, by = 2)) +
  labs(x = "Year", y = "Temperature (°C)", color = NULL, shape = NULL) +
  theme_classic(base_size = 25) +
  theme(legend.text = element_text(size = 12),
        legend.position = "top")

ggsave(plot = temperatures.avg, "output/temperatures.paper1.tiff", width = 8, height = 7 , unit = "in", dpi = 600, bg= "white")

yield_table <- data_model %>%
  group_by(location) %>%
  summarise(
    Years = n_distinct(year),
    Mean  = mean(yieldkg.ha, na.rm = TRUE),
    Min   = min(yieldkg.ha, na.rm = TRUE),
    Max   = max(yieldkg.ha, na.rm = TRUE),
    Range = Max - Min,
    .groups = "drop") %>%
  arrange(desc(Mean))
###########
library(broom)

df <- read_csv("input/historical_weather.csv")
df1 <- df %>%
  group_by(year) %>%
  summarise(
    mean_maxt = mean(mean_maxt, na.rm = TRUE),
    mean_mint = mean(mean_mint, na.rm = TRUE),
    avg_t     = mean(avg_t, na.rm = TRUE),
    .groups = "drop") %>%
  pivot_longer(-year, names_to = "variable", values_to = "temperature") %>%
  mutate(variable = dplyr::recode(variable,
                           mean_maxt = "Maximum T",
                           mean_mint = "Minimum T",
                           avg_t     = "Mean T"))
slopes <- df1 %>%
  group_by(variable) %>%
  summarise(
    slope = coef(lm(temperature ~ year))[2],
    .groups = "drop") %>%
  mutate(label = paste0("Slope = ", round(slope*10, 2), " °C/decade"))

label_pos <- df1 %>%
  group_by(variable) %>%
  filter(year == max(year)) %>%
  left_join(slopes, by = "variable")

 T1 <- ggplot(df1, aes(year, temperature,
                    color = variable, shape = variable)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
  scale_y_continuous(breaks = seq(18, 38, by = 2)) +
  scale_x_continuous(breaks = seq(1980, 2025, by = 5)) +
   geom_text(data = label_pos,
             aes(x = year - 5,y = temperature + 0.7, label = label),
             show.legend = FALSE,
             size = 6) +
  labs(x = "Year", y = "Temperature (°C)",
       color = NULL, shape = NULL) +
  theme_classic(base_size = 25) +
  expand_limits(x = max(df1$year) + 2)

ggsave(plot = T1, "output/historical.ES.png", width = 12, height = 7 , unit = "in", dpi = 600, bg= "white")

####### June-July 
df <- read_csv("input/historical_weather.csv")
df_jj <- df %>%
  filter(month %in% c("june", "july")) %>%
  group_by(year) %>%
  summarise(mean_maxt = mean(mean_maxt, na.rm = TRUE),
            mean_mint = mean(mean_mint, na.rm = TRUE),
            avg_t = mean(avg_t, na.rm = TRUE),
            .groups = "drop") %>%
  pivot_longer(cols = -year, names_to = "variable", values_to = "temperature") %>%
  mutate(variable = dplyr::recode(variable,
                           mean_maxt = "June–July Maximum T",
                           mean_mint = "June–July Minimum T",
                           avg_t = "June–July Mean T"))

slopes_jj <- df_jj %>%
  group_by(variable) %>%
  summarise(slope = coef(lm(temperature ~ year, data = cur_data()))[2],
            .groups = "drop") %>%
  mutate(label = paste0("+", round(slope*10, 2), " °C/decade"))
label_pos_jj <- df_jj %>%
  group_by(variable) %>%
  filter(year == max(year)) %>%
  summarise(temperature = max(temperature), .groups = "drop") %>%
  left_join(slopes_jj, by = "variable") %>%
  arrange(desc(temperature)) %>%
  mutate(x_pos = max(df_jj$year) - 1,
    y_pos = temperature + seq(0.5, 1.5, length.out = n()))


 T2 <- ggplot(df_jj, aes(year, temperature, color = variable, shape = variable)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
  geom_text(data = tibble(
    variable = c("June–July Maximum T", "June–July Mean T", "June–July Minimum T"),
    label = slopes_jj$label,  
    x_pos = c(2018, 2018, 2018),
    y_pos = c(35, 29, 24)),
  aes(x = x_pos, y = y_pos, label = label, color = variable),
  show.legend = FALSE,
  size = 6,
  hjust = 0) +
  geom_vline(xintercept = c(2010, 2024), linetype = "dashed", linewidth = 0.7, color = "black") +
  scale_y_continuous(breaks = seq(18, 38, by = 2)) +
  scale_x_continuous(breaks = seq(1980, 2025, by = 5)) +
  labs(x = "Year", y = "Temperature (°C)", color = NULL, shape = NULL) +
  theme_classic(base_size = 25) +
  theme(legend.position = "top") +
  guides(color = guide_legend(override.aes = list(linewidth = 2.5))) +
  expand_limits(x = max(df_jj$year) + 3) +
  coord_cartesian(clip = "off")


####################################
 # Figure S6
df <- read_csv("input/historical_weather.csv")

df_jj <- df %>%
  filter(month %in% c("june", "july")) %>%
  group_by(year) %>%
  summarise(
    sd_maxt   = sd(mean_maxt, na.rm = TRUE),
    sd_mint   = sd(mean_mint, na.rm = TRUE),
    sd_avgt   = sd(avg_t, na.rm = TRUE),
    mean_maxt = mean(mean_maxt, na.rm = TRUE),
    mean_mint = mean(mean_mint, na.rm = TRUE),
    avg_t     = mean(avg_t, na.rm = TRUE),
    .groups = "drop")

df_jj_long <- df_jj %>%
  pivot_longer(
    cols = c(mean_maxt, mean_mint, avg_t, sd_maxt, sd_mint, sd_avgt),
    names_to = "variable", values_to = "value") %>%
  mutate(
    type = ifelse(grepl("^sd_", variable), "sd", "mean"),
    var  = case_when(
      variable %in% c("mean_maxt", "sd_maxt") ~ "June–July Maximum T",
      variable %in% c("mean_mint", "sd_mint") ~ "June–July Minimum T",
      variable %in% c("avg_t",     "sd_avgt") ~ "June–July Mean T")) %>%
  select(-variable) %>%
  pivot_wider(names_from = type, values_from = value) %>%
  rename(variable = var, temperature = mean)

slopes_jj <- df_jj_long %>%
  group_by(variable) %>%
  summarise(slope = coef(lm(temperature ~ year, data = cur_data()))[2],
            .groups = "drop") %>%
  mutate(label = paste0("+", round(slope*10, 2), " °C/decade"))

label_pos_jj <- df_jj_long %>%
  group_by(variable) %>%
  filter(year == max(year)) %>%
  summarise(temperature = max(temperature), .groups = "drop") %>%
  left_join(slopes_jj, by = "variable") %>%
  arrange(desc(temperature)) %>%
  mutate(x_pos = max(df_jj_long$year) - 1,
         y_pos = temperature + seq(0.5, 1.5, length.out = n()))

S6 <- ggplot(df_jj_long, aes(year, temperature, color = variable, shape = variable, fill = variable)) +
  geom_ribbon(aes(ymin = temperature - sd, ymax = temperature + sd), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
  geom_text(data = tibble(
    variable = c("June–July Maximum T", "June–July Mean T", "June–July Minimum T"),
    label = slopes_jj$label,
    x_pos = c(2018, 2018, 2018),
    y_pos = c(35, 29, 24)),
    aes(x = x_pos, y = y_pos, label = label, color = variable),
    show.legend = FALSE, size = 6, hjust = 0) +
  geom_vline(xintercept = c(2010, 2024), linetype = "dashed", linewidth = 0.7, color = "black") +
  scale_y_continuous(breaks = seq(18, 38, by = 2)) +
  scale_x_continuous(breaks = seq(1980, 2025, by = 5)) +
  labs(x = "Year", y = "Temperature (°C)", color = NULL, shape = NULL, fill = NULL) +
  theme_classic(base_size = 25) +
  theme(legend.position = "top") +
  guides(color = guide_legend(override.aes = list(linewidth = 2.5))) +
  expand_limits(x = max(df_jj_long$year) + 3) +
  coord_cartesian(clip = "off")

ggsave("output/FigureS6.png",
       plot = S6,
       width = 12, height = 7 ,
       units = "in", dpi = 600, bg = "white")


####### N of days with T above 35
#avg number of days with T >35 in each month
df_month_avg <- df %>%
  filter(month %in% c("april","may","june","july","august","september")) %>%
  group_by(year, month) %>%
  summarise(
    avg_above_35 = mean(above_35, na.rm = TRUE),
    .groups = "drop")
# number of days with T>35 in each season
df_season <- df_month_avg %>%
  group_by(year) %>%
  summarise(
    season_above_35 = sum(avg_above_35),
    .groups = "drop")

slope_days <- coef(lm(season_above_35 ~ year, data = df_season))[2]

label_days <- data.frame(
  year = max(df_season$year),
  season_above_35 = df_season$season_above_35[df_season$year == max(df_season$year)],
  label = paste0("Slope = ", round(slope_days*10, 1), " days/decade"))

 T3 <- ggplot(df_season, aes(x = year, y = season_above_35)) +
  geom_line(color = "firebrick", linewidth = 1) +
  geom_point(color = "firebrick", size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1.2) +
  geom_text(data = label_days,
            aes(x = year - 5,
                y = 65,
                label = label),
            size = 6,
            color = "black") +
  scale_y_continuous(breaks = seq(0, max(df_season$season_above_35)+5, by = 5)) +
  scale_x_continuous(breaks = seq(1980, 2025, by = 5)) +
  labs(x = "Year", y = "N Days > 35°C ") +
  theme_classic(base_size = 25) +
  expand_limits(x = max(df_season$year) + 2)

ggsave("output/days_above_35_seasonal.png",
       plot = T3,
       width = 8, height = 7,
       units = "in", dpi = 600, bg = "white")


#######################################################
df <- read.csv("input/Corn data organized2.csv", stringsAsFactors = FALSE)
data_model <- df %>%
  select(yield, year, location, hybrid, mean_mint) %>% 
  mutate(yield= as.numeric(yield)) %>% 
  drop_na() 
env_summary <- data_model %>%
  group_by(year, location) %>%
  summarise(
    yield     = mean(yield, na.rm = TRUE),
    mean_mint = mean(mean_mint, na.rm = TRUE),
    .groups   = "drop")

p1 <- ggplot(env_summary, aes(x = mean_mint, y = yield)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "loess", se = TRUE)+
  labs(x = "Mean night temperature (°C, June–July)",
       y = "Yield (t/ha)") +
  theme_classic(base_size = 16)
ggsave(plot = p1, "output/nt.against.yield.png",
       width = 12, height = 7 , unit = "in", dpi = 600, bg= "white")

p2 <- ggplot(env_summary, aes(x = mean_mint, y = yield)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE)+
  labs(x = "Mean night temperature (°C, June–July)",
       y = "Yield (t/ha)") +
  theme_classic(base_size = 16)
ggsave(plot = p2, "output/nt.against.yield2.png",
       width = 12, height = 7 , unit = "in", dpi = 600, bg= "white")
