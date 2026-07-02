library(dplyr)
library(tidyverse)
library(nlme)
library(glmnet)
data <- read.csv("input/Corn data organized2.csv")
data$yield <- as.numeric(data$yield)
data <- data %>%
  select(location,year,sowing.doy, crm, yield, mean_maxt,mean_mint, cum_radn.j, above_35,vpd,april, may,june,july,august,sept) %>% 
  na.omit()

Y <- data %>% select(yield) %>% as.matrix()
X <- model.matrix(
  yield ~ crm + mean_maxt + mean_mint +
    cum_radn.j + above_35 + vpd + april + may + june + july +
    august + sept,
  data = data)[,-1]
lambdas <- 10^seq(-3, 3, length.out = 100)
lasso.fit <- glmnet(X, Y, alpha = 1, lambda = lambdas)
 plot(lasso.fit, xvar = "lambda")

 
lasso.cv <- cv.glmnet(X, Y, alpha = 1, lambda = lambdas,
                       nfolds = 10)
plot(lasso.cv) 

lasso.best <- glmnet(X, Y, alpha = 1,
                     lambda = lasso.cv$lambda.min)
coef(lasso.best)
lasso.coef <- coef(lasso.best)
coef_df <- data.frame(
  variable = rownames(lasso.coef),
  coef = as.numeric(lasso.coef))

coef_df$abs_coef <- abs(coef_df$coef)
#eliminate intercept and rank them
coef_df <- coef_df %>%
  filter(variable != "(Intercept)") %>%
  arrange(desc(abs_coef))
#plot
ggplot(coef_df, aes(x = reorder(variable, abs_coef), y = abs_coef, fill = abs_coef)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "lightblue", high = "steelblue") +
  labs(x = "Variable", y = "Importancia (valor absoluto del coeficiente)",
       title = "Ranking de variables según importancia LASSO") +
  theme_minimal()


###new approach scaling weather variables


clim_vars <- c("mean_maxt","mean_mint","cum_radn.j","above_35","vpd",
               "april","may","june","july","august","sept")


data[clim_vars] <- scale(data[clim_vars])

# Matriz for glmnet
X <- model.matrix(
  yield ~  mean_mint +
    cum_radn.j + vpd + april + may + june + july +
    august + sept,
  data = data)[,-1]

# Response
Y <- data %>% select(yield) %>% as.matrix()

# Lasso con CV
lambdas <- 10^seq(-3, 3, length.out = 100)
lasso.cv <- cv.glmnet(X, Y, alpha = 1, lambda = lambdas, nfolds = 10)
lasso.best <- glmnet(X, Y, alpha = 1, lambda = lasso.cv$lambda.min)

# Coefficients
lasso.coef <- coef(lasso.best)
coef_df <- data.frame(
  variable = rownames(lasso.coef),
  coef = as.numeric(lasso.coef))
coef_df$abs_coef <- abs(coef_df$coef)
coef_df <- coef_df %>% filter(variable != "(Intercept)") %>% arrange(desc(abs_coef))


var_labels <- c(
 # mean_maxt   = "Avg max T",
  mean_mint   = "Avg min T",
  cum_radn.j  = "Cumulative Radiation",
 # above_35    = "Days >35°C",
  april       = "Apr (mm)",
  may         = "May (mm)",
  june        = "Jun (mm)",
  july        = "Jul (mm)",
  august      = "Aug (mm)",
  sept        = "Sep (mm)",
  #crm         = "CRM",
  vpd         = "VPD (kPa)")

coef_df$variable_label <- var_labels[coef_df$variable]


# Plot
p1 <- ggplot(coef_df, aes(x = reorder(variable_label, abs_coef), 
                           y = abs_coef, fill = abs_coef)) +
  geom_col() +
   coord_flip() +
   scale_fill_gradient(low = "lightblue", high = "steelblue") +
   labs(x = "Variable", y = "Absolute standardized effect size") +
   theme_classic(base_size = 25)
ggsave(plot = p1, "output/Lasso.tiff", width = 8, height = 7 , unit = "in", dpi = 600, bg= "white")


p2 <- ggplot(coef_df, aes(
  x = reorder(variable_label, coef),
  y = coef,
  fill = coef > 0)) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    values = c("TRUE" = "steelblue", "FALSE" = "tomato"),
    labels = c("Negative", "Positive"),
    name = "Effect direction") +
  coord_flip() +
  scale_y_continuous(breaks = seq(-0.4, 0.4, by = 0.1))+
  labs(x = "Variable", y = "Standardized effect size") +
  theme_classic(base_size = 22) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    panel.grid.major.x = element_line(color = "grey85"),
    panel.grid.minor = element_blank(),
    legend.position = "top")
ggsave(plot = p2, "output/Lasso2.tiff", width = 8, height = 7 , unit = "in", dpi = 600, bg= "white")
