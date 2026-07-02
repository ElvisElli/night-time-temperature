library(nlme)
library(ggplot2)
library(nlraa)

## reading the raw data
dat <- read.csv('../raw-data/Corn data organized-10-16-2024.csv')
meta <- read.csv('../raw-data/metadata.csv')
## dimensions of the data set
dim(dat)

## checking the variable names
names(dat)

## making sure that the numeric columns are, in fact, numeric
num.cols <- c('yield', 'CRM', 'SowingdateDOY')
dat[c(num.cols)] <- lapply(dat[num.cols], function(x) as.numeric(x))

## same for factor columns
fac.cols <- c('Location', 'Hybrid', 'Year', 'row_spacing')
dat[c(fac.cols)] <- lapply(dat[fac.cols], function(x) as.factor(x))

##  Looking at the data for NAs
summary(dat)

## it seems like the original data set contains a couple
## NAs in the yield and CRM columns. I am removing those here
dat <- na.omit(dat)


## calculating relative yield
dat$ry <- ave(dat$yield, dat$Location, dat$Year,
              FUN = function(x) 100 * x / max(x, na.rm = TRUE) )

## looking at correlations between variables
pairs(dat[sapply(dat, is.numeric)])



## the trend with CRM is visible even without controlling
## for other factors
with(dat, plot(ry ~ CRM, col = 4))

## the relationsip with planting data is not evident at first
with(dat, plot(ry ~ SowingdateDOY, col = 4))

## i am thinking that row spacing might have to be categorical
with(dat, plot(ry ~ row_spacing, col = 4))


## This model gives us a couple warnings that
## we might want to pay attention to.

## When is says that the fixed-effect model matrix is rank
## deficient, it is letting us know that this model is
## most likely overparameterized

fm.initial <- lme4::lmer(ry ~ CRM * SowingdateDOY +
                           Location * SowingdateDOY +
                           Year * SowingdateDOY +
                           CRM * Location +
                           row_spacing +
                           (1 | Hybrid),
                         data = dat)


## exploratory models

## we could fit these same models using lme4 or nlme
## here, I will follow with the nlme package simply because
## I believe that it offers some correlation structures for
## random effects that need to be handbuilt in lme4.

fmm1.1 <- lme4::lmer(yield ~ CRM * SowingdateDOY+
                       (1 | Location),
                     data = dat)

## Starting with a simple model in which we
## investigate the effect of CRM and Sowing date

fmm1 <- lme(yield ~ CRM * SowingdateDOY,
            random = ~ 1 | Location,
            data = dat)

## Adding year as a random effect nested within location
## seems to improve the AIC a lot.
fmm2 <- lme(yield ~ CRM * SowingdateDOY,
            random = ~ 1 | Year/Location,
            data = dat)

IC_tab(fmm1, fmm2)

## Adding hybrid as an effect within location and within ## year seems to improve the model even more

fmm3 <- lme(yield ~ CRM * SowingdateDOY,
            random = ~ 1 | Year/Location/Hybrid,
            data = dat)

IC_tab(fmm1, fmm2, fmm3)

## This is the same model but with a different syntax
## we can specify different structures at different levels
## of the random effects. This might be important later.

fmm3.1 <- lme(yield ~ CRM * SowingdateDOY,
            random = list(Year = ~ 1,
                          Location = ~ 1,
                          Hybrid = ~ 1 ),
            data = dat)

IC_tab(fmm3, fmm3.1)


## I think, for now, we can leave the random effects
## alone and start focusing on the fixed effects. We
## will address the random effects again later.
## I am under the impression that temperature, radiation,
## and preciptation might be important here as well. Let's
## look at them and see how they help us predict yield.


pairs(dat[c('cum_radn', 'cum_rain', 'mean_maxt', 'mean_mint')])

fmm4.1 <- lme(yield ~ CRM * SowingdateDOY * cum_rain,
            random = list(Year = ~ 1,
                          Location = ~ 1,
                          Hybrid = ~ 1 ),
            data = dat)


## Adding cumulative rain to the model didn't really
## improve our model. For now, we can focus on other variables.
## this is a little surprising though... were these trials irrigated?
IC_tab(fmm1, fmm2, fmm3, fmm4.1)


fmm4.2 <- lme(yield ~ CRM * SowingdateDOY * cum_radn,
            random = list(Location = ~ 1,
                          Year = ~ 1,
                          Hybrid = ~ 1 ),
            data = dat)

## Adding cumulative radiation to the model didn't really
## improve our model. For now, we can focus on other variables.
IC_tab(fmm1, fmm2, fmm3, fmm4.2)

fmm4.3 <- lme(yield ~ CRM * SowingdateDOY * mean_maxt,
            random = list(Location = ~ 1,
                          Year = ~ 1,
                          CRM = ~ 1 ),
            data = dat)

## Adding mean maximum temperature seems to have improved
## our model
IC_tab(fmm1, fmm2, fmm3, fmm4.3)

## checking that the predictor variables aren't correlated
pairs(dat[c('CRM', 'SowingdateDOY', 'mean_maxt')])


fmm4.4 <- lme(yield ~ CRM * SowingdateDOY * mean_maxt *
                mean_mint,
            random = list(Location = ~ 1,
                          Year = ~ 1,
                          CRM = ~ 1 ),
            data = dat)

## Adding cumulative radiation to the model didn't really
## improve our model. For now, we can focus on other variables.
IC_tab(fmm1, fmm2, fmm3, fmm4.3, fmm4.4)

## checking that the predictor variables aren't correlated.
## it seems like max and min temperatures are highly correlated
pairs(dat[c('CRM', 'SowingdateDOY', 'mean_maxt', 'mean_mint')])

## this can be a problem for the model
cor(dat$mean_maxt, dat$mean_mint)


fmm4.5 <- lme(yield ~ CRM * SowingdateDOY * mean_mint,
              random = list(Location = ~ 1,
                            Year = ~ 1,
                            CRM = ~ 1 ),
              data = dat)

## For now, I think it's best we continue with only one
## of them in the model. It seems that the model that
## contains only the max temp is more adequate than
## the one with only min temp
IC_tab(fmm1, fmm2, fmm3, fmm4.3, fmm4.4, fmm4.5)




fmm5 <- lme(yield ~ CRM * SowingdateDOY *
              mean_maxt * row_spacing,
            random = list(Location = ~ 1,
                          Year = ~ 1,
                          CRM = ~ 1 ),
            data = dat)
IC_tab(fmm1, fmm2, fmm3, fmm4.5, fmm5)


random.effects(fmm5, level = 1)



emmeans::emmeans(fmm5, specs =  ~ CRM, at = list(CRM = 100:120))


met <- apsimx::get_iem_apsim_met(c(-89, 35),
                                 dates = c('2025-01-01', '2025-02-01'))
attributes(met)
