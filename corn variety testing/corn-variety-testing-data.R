
## Script name: ====

## Script objective: ====

## Author: Elvis F. Elli ====

## Script was created on: 2024-06-05====

## Cleaning up environment ====
rm(list=ls())

## Libraries ====
library(rstudioapi)
library(tidyverse)
library(readxl)
library(smatr)
library(nlstools)  # for residual plots
library(ggplot2)
library(car)
library(lme4)
library(lattice)
library(lmerTest)


## Set working directory ====
setwd(dirname(getActiveDocumentContext()$path))

##data
data <- read_excel("../../data/corn variety testing/Corn data organized-10-16-2024.xlsx")
#####
data.summary <- data %>% 
  group_by(`SowingdateDOY`) %>% 
  summarise(p90=quantile(yield,probs = 0.9,na.rm = T))

data.summary4 <- data %>% 
  group_by(CRM) %>% 
  summarise(p90=quantile(yield,probs = 0.9,na.rm = T))

data2 <- data %>% 
  left_join(data.summary) %>% 
  mutate(abovep90 = ifelse(yield>p90,"yes","no"))

data3 <- data2 %>% filter (abovep90=="yes")

data4 <- data %>% 
  left_join(data.summary4) %>% 
  mutate(abovep90 = ifelse(yield>p90,"yes","no"))

data5 <- data4%>% filter(yield>p90)
#####
mod1 <- sma(yield~`SowingdateDOY`, data = data3)

coef(mod1)[1] #beta0
coef(mod1)[2] #beta1
mod1$r2[[1]]
label <- paste0("y = ",round(coef(mod1)[1],2), " + ", round(coef(mod1)[2],2), "x") 

mod2 <- sma(yield~ CRM, data =data4%>% filter (abovep90=="yes"))
coef(mod2) [1] #beta0
coef(mod2) [2] #beta1
mod2$r2[[1]]
label1 <- paste0 ("y=",round(coef(mod2)[1],2), "+",round(coef(mod2)[2],2), "x")

p <- 
  data %>% 
 ggplot(aes(x=Year,y=yield,colour=`SowingdateDOY`)) +
  geom_point()+
  facet_wrap(~Location) 

  ggsave(plot = p, "../output/Figure 0.png",width = 8,height = 7, units = "in", dpi = 600)


  dfMod[,c("data3")] %>% 
    unnest(cols = c(data)) %>% 
    ggplot(aes(x=`SowingdateDOY`,y=yield))+
    geom_point(shape = 21, size=4)+
    geom_line(data = dfMod[,c("pred")] %>% 
                unnest(cols = c(pred)),
              aes(x = CPdif, y = .), linewidth = 1)

             
#mod2 <- sma(yield~`SowingdateDOY`, data = filtered_data_yield )

#coef(mod2)[1] #beta0
#coef(mod2)[2] #beta1
#mod1$r2[[1]]
#label2 <- paste0("y = ",round(coef(mod2)[1],2), " + ", round(coef(mod2)[2],2), "x") 



##Quadratic model (ASA)====
  

str(data3)
  plot(data3$`SowingdateDOY`,data3$yield,pch=16)
  plot(data4$CRM, data4$yield)
  plot(data5$CRM,data5$yield)
 `SowingdateDOY2` <- data3$`SowingdateDOY`^2

 wak2 <- lm(yield~`SowingdateDOY`+`SowingdateDOY2`, data = data3)
 summary(wak2)
 
 `SowingdateDOYvalues` <- seq(from =80, to =140, by =1)
 yieldpredict <- predict(wak2,list(`SowingdateDOY`= `SowingdateDOYvalues`,
                                   `SowingdateDOY2` = `SowingdateDOYvalues`^2))
 summary(yieldpredict)
 
 p25 <- 
   ggplot()+
   geom_point(data= data,
     aes(x= `SowingdateDOY`, y= yield), size = 1, alpha = .2,colour="grey")+
 geom_point(data= data3,
            aes(x=`SowingdateDOY`, y= yield),
            fill= "blue",shape=21, size= 3, alpha=.4) +
 geom_line(aes(x=`SowingdateDOYvalues`, y =  yieldpredict),
           col = "red",linewidth = 1) +
   scale_x_continuous(breaks= seq(80,140, by= 10))+
   scale_y_continuous(breaks = seq(6000,16000, by = 2000))+
   # scale_y_continuous(breaks = seq(min(data$yield,na.rm= TRUE), max(data$yield,na.rm= TRUE), by = 2000))+
   geom_point(aes(x = 126, y = max(data3$yield) ), 
              fill = "blue", shape = 21, size = 5) +
   geom_text(aes(x = 128, y = max(data3$yield), 
                 label = "Yields > 90%"), 
             hjust = 0, size= 6) +
   geom_point(aes(x = 126, y = max(data3$yield) -450 ), 
              fill = "grey", shape = 21, size = 5) +
   geom_text(aes(x = 128, y = max(data3$yield) - 450, 
                 label = "Total Data"), 
             hjust = 0,size = 6)+
   #geom_text(aes(x=90, y= 310,label = "p=3.6^-9"),
           #  size = 8 ,alpha=.40)+
 labs (x="Sowingdate-DOY", y = expression(paste("Grain yield (",Kg~ha^-1,")"))) +
 theme_classic(base_size = 25)
 ggsave(plot = p25, "../output/Figure 25.png",
              width = 8, height = 7 , unit = "in", dpi = 600)
 
##### 
 ##1:Marianna
 
 marianna.df <- data %>% 
   filter(Location=="Marianna")
 
 ##inspection
 marianna.df %>% ggplot(aes(x=SowingdateDOY,y=yield))+
   geom_point()+
   facet_wrap(~Location)
 
 sd.marianna <- marianna.df %>% 
   group_by(`SowingdateDOY`) %>% 
   summarise(p90=quantile(yield,probs = 0.9,na.rm = T))
 
 rm.marianna <- marianna.df %>%
   group_by(CRM) %>% 
   summarise(p90=quantile(yield,probs = 0.9,na.rm = T))
 
 marianna.df2 <- marianna.df %>% 
   left_join(sd.marianna) %>% 
   mutate(abovep90 = ifelse(yield>p90,"yes","no")) %>% 
   filter(!is.na(CRM))
 
 marianna.model.df <- marianna.df2 %>% filter (abovep90=="yes")
 
 
 marianna.model <- lm(yield~SowingdateDOY+I(SowingdateDOY^2), data = marianna.model.df)
 summary(marianna.model)
 coefficients(marianna.model)
 
 SowingdateDOYvalues <- seq(from =80, to =140, by =1)
 
 marianna.yield.predict <- predict(marianna.model,list(SowingdateDOY= SowingdateDOYvalues))
 
sd.marianna.plot <- 
   ggplot()+
   geom_point(data= marianna.df,
              aes(x= SowingdateDOY, y= yield), size = 1.5, alpha = .2)+
   geom_point(data= marianna.model.df,
              aes(x=SowingdateDOY, y= yield),
              fill= "blue",shape=21, size= 2, alpha=.4) +
   geom_line(aes(x=`SowingdateDOYvalues`, y =  marianna.yield.predict),
             col = "red",linewidth = 1) +
   geom_text(aes(x=90, y= 310,label = "p=3.6^-9"),
             size = 8 ,alpha=.40)+
   labs (x="SowingdateDOY", y = expression(paste("Grain yield (",Bu~a^-1,")"))) +
   theme_classic(base_size = 25)+
  scale_y_continuous(limits = c(100,320))+
  ggtitle("Marianna")
 ggsave(plot = sd.marianna.plot, "../output/sd.marianna.plot.png",
        width = 8, height = 7 , unit = "in", dpi = 600)
 
 ##2:Keiser
 
 keiser.df <- data %>% 
   filter(Location=="Keiser")
 
 ##inspection
 keiser.df %>% ggplot(aes(x=SowingdateDOY,y=yield))+
   geom_point()+
   facet_wrap(~Location)
 
 sd.keiser <- keiser.df %>% 
   group_by(`SowingdateDOY`) %>% 
   summarise(p90=quantile(yield,probs = 0.9,na.rm = T))
 
 rm.keiser <- keiser.df %>%
   group_by(CRM) %>% 
   summarise(p90=quantile(yield,probs = 0.9,na.rm = T))
 
 keiser.df2 <- keiser.df %>% 
   left_join(sd.keiser) %>% 
   mutate(abovep90 = ifelse(yield>p90,"yes","no")) %>% 
   filter(!is.na(CRM))
 
 keiser.model.df <- keiser.df2 %>% filter (abovep90=="yes")
 
 
 keiser.model <- lm(yield~SowingdateDOY+I(SowingdateDOY^2), data = keiser.model.df)
 summary(keiser.model)
 coefficients(keiser.model)
 
 SowingdateDOYvalues <- seq(from =80, to =140, by =1)
 
 keiser.yield.predict <- predict(keiser.model,list(SowingdateDOY= SowingdateDOYvalues))
 
 sd.keiser.plot <- 
 ggplot()+
   geom_point(data= keiser.df,
              aes(x= SowingdateDOY, y= yield), size = 1.5, alpha = .2)+
   geom_point(data= keiser.model.df,
              aes(x=SowingdateDOY, y= yield),
              fill= "blue",shape=21, size= 2, alpha=.4) +
   geom_line(aes(x=`SowingdateDOYvalues`, y =  keiser.yield.predict),
             col = "red",linewidth = 1) +
   geom_text(aes(x=90, y= 310,label = "p=3.6^-9"),
             size = 8 ,alpha=.40)+
   labs (x="SowingdateDOY", y = expression(paste("Grain yield (",Bu~a^-1,")"))) +
   theme_classic(base_size = 25)+
   scale_y_continuous(limits = c(100,320))+
   ggtitle("Keiser")
 ggsave(plot = sd.keiser.plot, "../output/sd.keiser.plot.png",
        width = 8, height = 7 , unit = "in", dpi = 600)
 
 
 ##3:Rohwer
 
 Rohwer.df <- data %>% 
   filter(Location=="Rohwer")
 
 ##inspection
 Rohwer.df %>% ggplot(aes(x=SowingdateDOY,y=yield))+
   geom_point()+
   facet_wrap(~Location)
 
 sd.Rohwer <- Rohwer.df %>% 
   group_by(`SowingdateDOY`) %>% 
   summarise(p90=quantile(yield,probs = 0.9,na.rm = T))
 
 rm.Rohwer <- Rohwer.df %>%
   group_by(CRM) %>% 
   summarise(p90=quantile(yield,probs = 0.9,na.rm = T))
 
 Rohwer.df2 <- Rohwer.df %>% 
   left_join(sd.Rohwer) %>% 
   mutate(abovep90 = ifelse(yield>p90,"yes","no")) %>% 
   filter(!is.na(CRM))
 
 Rohwer.model.df <- Rohwer.df2 %>% filter (abovep90=="yes")
 
 Rohwer.model <- lm(yield ~ SowingdateDOY + I(SowingdateDOY^2), data = Rohwer.model.df)
 summary(Rohwer.model)
 coefficients(Rohwer.model)
 
 SowingdateDOYvalues <- seq(from =80, to =140, by =1)
 
 Rohwer.yield.predict <- predict(Rohwer.model,list(SowingdateDOY= SowingdateDOYvalues))
 
 sd.Rohwer.plot <- 
   ggplot()+
   geom_point(data= Rohwer.df,
              aes(x= SowingdateDOY, y= yield), size = 1.5, alpha = .2)+
   geom_point(data= Rohwer.model.df,
              aes(x=SowingdateDOY, y= yield),
              fill= "blue",shape=21, size= 2, alpha=.4) +
   geom_line(aes(x=`SowingdateDOYvalues`, y =  Rohwer.yield.predict),
             col = "red",linewidth = 1) +
   geom_text(aes(x=90, y= 310,label = "p=3.6^-9"),
             size = 8 ,alpha=.40)+
   labs (x="SowingdateDOY", y = expression(paste("Grain yield (",Bu~a^-1,")"))) +
   theme_classic(base_size = 25)+
   scale_y_continuous(limits = c(100,320))+
   ggtitle("Rohwer")
 ggsave(plot = sd.Rohwer.plot, "../output/sd.Rohwer.plot.png",
        width = 8, height = 7 , unit = "in", dpi = 600)
 
 
 
 
 ##4:stutgart
 
 Stuttgart.df <- data %>% 
   filter(Location=="Stuttgart")
 
 ##inspection
 Stuttgart.df %>% ggplot(aes(x=SowingdateDOY,y=yield))+
   geom_point()+
   facet_wrap(~Location)
 
 sd.Stuttgart <- Stuttgart.df %>% 
   group_by(`SowingdateDOY`) %>% 
   summarise(p90=quantile(yield,probs = 0.9,na.rm = T))
 
 rm.Stuttgart <- Stuttgart.df %>%
   group_by(CRM) %>% 
   summarise(p90=quantile(yield,probs = 0.9,na.rm = T))
 
 Stuttgart.df2 <- Stuttgart.df %>% 
   left_join(sd.Stuttgart) %>% 
   mutate(abovep90 = ifelse(yield>p90,"yes","no")) %>% 
   filter(!is.na(CRM))
 
 Stuttgart.model.df <- Stuttgart.df2 %>% filter (abovep90=="yes")
 
 
 Stuttgart.model <- lm(yield~SowingdateDOY+I(SowingdateDOY^2), data = Stuttgart.model.df)
 summary(Stuttgart.model)
 coefficients(Stuttgart.model)
 anova(Stuttgart.model)
 
 SowingdateDOYvalues <- seq(from =80, to =140, by =1)
 
 Stuttgart.yield.predict <- predict(Stuttgart.model,list(SowingdateDOY= SowingdateDOYvalues))
 
 sd.Stuttgart.plot <- 
   ggplot()+
   geom_point(data= Stuttgart.df,
              aes(x= SowingdateDOY, y= yield), size = 1.5, alpha = .2)+
   geom_point(data= Stuttgart.model.df,
              aes(x=SowingdateDOY, y= yield),
              fill= "blue",shape=21, size= 2, alpha=.4) +
   geom_line(aes(x=`SowingdateDOYvalues`, y =  Stuttgart.yield.predict),
             col = "red",linewidth = 1) +
   geom_text(aes(x=90, y= 310,label = "p=3.6^-9"),
             size = 8 ,alpha=.40)+
   labs (x="SowingdateDOY", y = expression(paste("Grain Yield (",Bu~a^-1,")"))) +
   theme_classic(base_size = 25)+
 scale_y_continuous(limits = c(100,320))+
   ggtitle("Stuttgart")
 ggsave(plot = sd.Stuttgart.plot, "../output/sd.Stuttgart.plot.png",
        width = 8, height = 7 , unit = "in", dpi = 600)
 
 
 
 perlocation.plot <- plot_grid(sd.marianna.plot,sd.Rohwer.plot,sd.keiser.plot,sd.Stuttgart.plot,labels = "AUTO")
 
 final.plot <- plot_grid(p15,perlocation.plot)
 
 ggsave(plot = final.plot, "../output/final.plot.tiff",
        width = 20, height = 10, unit = "in", dpi = 600)
##### 
 
 ## Impact of Sowing Date on Grain yield by CRM 
 
 CRM2 <- data5$CRM^2
 
 wak3 <-lm(yield~CRM + CRM2, data = data5)
 summary(wak3)
 crmvalues <- seq(from= 100, to= 120, by= 1) 
 yieldpredict2 <- predict(wak3,list(CRM= CRMvalues,
                                    CRM2=CRMvalues^2))
 model_summary <- summary(wak3)
 coefifficients_model <- model_summary$coefficients
 print(coefifficients_model)
 
 beta_0 <- coefifficients_model[1,1]
 beta_1 <- coefifficients_model[2,1]
 beta_2 <- coefifficients_model[3,1]
 
 linear_model <- lm(yield ~ CRM, data = data5)
 summary(linear_model)
 quadratic_model <- lm(yield ~ CRM2, data = data5)
 summary(quadratic_model)
 label_crm_yield <- paste("Yield =", round(beta_0, 2), "+", round(beta_1, 2), "* CRM +", round(beta_2, 2), "* CRM^2")
 print(label_crm_yield)
 r2 <- model_summary$r.squared
 
 label_r2 <-paste("R² =", round(r2, 3))
 

 
p26 <-
ggplot()+
  geom_point(aes(x= data$CRM, y= data$yield), fill= "grey", size = 1, alpha= .1)+
  geom_point(aes(x=data5$CRM,y=data5$yield), fill="blue",shape=21, size= 3, alpha= .4)+
  geom_line(aes(x=CRMvalues,y=yieldpredict2),col="red",linewidth=1)+
  labs(x="Corn Relative Maturity (days)",
       y= expression(paste("Grain yield (",kg~ha^-1,")"))) +
  geom_point(aes(x = 102, y = max(data3$yield) ), 
             fill = "blue", shape = 21, size = 5) +
  geom_text(aes(x = 103, y = max(data3$yield), 
                label = "Yields > 90%"), 
            hjust = 0, size= 6) +
  geom_point(aes(x = 102, y = max(data3$yield) -450 ), 
             fill = "grey", shape = 21, size = 5) +
 geom_text(aes(x = 103, y = max(data3$yield) - 450, 
                label = "Total Data"), 
            hjust = 0,size = 6)+
  theme_classic(base_size = 25)+
  scale_x_continuous(breaks= seq(100,125,by = 5))+
  scale_y_continuous(breaks = seq(6000,16000, by = 2000))
ggsave(plot = p26, "../output/Figure 26.png",
       width = 8, height = 7 , unit = "in", dpi = 600)

##1) marianna crm
marianna.df2.crm <- marianna.df %>% 
  left_join(rm.marianna) %>% 
  mutate(abovep90 = ifelse(yield>p90,"yes","no")) %>% 
  filter(!is.na(CRM))

marianna.model.df.crm <- marianna.df2.crm %>% filter(abovep90== "yes")

marianna.model.crm <- lm(yield~CRM + I(CRM^2),data = marianna.model.df.crm)

summary(marianna.model.crm)

coefficients(marianna.model.crm)

crmvalues <- seq(from= 100, to= 125, by= 1) 

marianna.yield.predict.crm <- predict(marianna.model.crm,list(CRM=crmvalues))

rm.marianna.plot <- 
  ggplot()+
  geom_point(data= marianna.df,
             aes(x= CRM, y= yield),fill= "grey", size= 1.2,alpha= .1)+
  geom_point(data=  marianna.model.df.crm,
             aes(x= CRM, y= yield), fill= 'blue',size= 2, shape= 21,alpha= .4)+
  geom_line(aes(x= crmvalues,y= marianna.yield.predict.crm),
            col= "red", linewidth= 1)+
  labs(x= "Corn Relative Maturity",y= expression(paste("Grain yield (",Bu~a^-1,")")))+
  theme_classic(base_size= 25)+
  ggtitle("Marianna")+
  scale_y_continuous(limits = c(125,300))
ggsave(plot = rm.marianna.plot, "../output/rm.marianna.plot.png",
       width = 8, height = 7 , unit = "in", dpi = 600)

##2) Keiser.crm

keiser.df2.crm <- keiser.df %>% 
  left_join(rm.keiser) %>% 
  mutate(abovep90 = ifelse(yield>p90,"yes","no")) %>% 
  filter(!is.na(CRM))

keiser.model.df.crm <- keiser.df2.crm %>% filter(abovep90== "yes")

keiser.model.crm <- lm(yield~CRM + I(CRM^2),data = keiser.model.df.crm)

summary(keiser.model.crm)

coefficients(keiser.model.crm)

crmvalues <- seq(from= 100, to= 125, by= 1) 

keiser.yield.predict.crm <- predict(keiser.model.crm,list(CRM=crmvalues))

rm.keiser.plot <- 
  ggplot()+
  geom_point(data= keiser.df,
             aes(x= CRM, y= yield),fill= "grey",size= 1.2,alpha= .1)+
  geom_point(data=  keiser.model.df.crm,
             aes(x= CRM, y= yield), fill= 'blue',size= 2, shape= 21,alpha= .4)+
  geom_line(aes(x= crmvalues,y= keiser.yield.predict.crm),
            col= "red", linewidth= 1)+
  labs(x= "Corn Relative Maturity",y= expression(paste("Grain yield (",Bu~a^-1,")")))+
  theme_classic(base_size= 25)+
  ggtitle("Keiser")+
  scale_y_continuous(limits = c(125,300))
ggsave(plot = rm.keiser.plot, "../output/rm.keiser.plot.png",
       width = 8, height = 7 , unit = "in", dpi = 600)

## Rohwer.crm 

Rohwer.df2.crm <- Rohwer.df %>% 
  left_join(rm.Rohwer) %>% 
  mutate(abovep90 = ifelse(yield>p90,"yes","no")) %>% 
  filter(!is.na(CRM))

Rohwer.model.df.crm <- Rohwer.df2.crm %>% filter(abovep90== "yes")

Rohwer.model.crm <- lm(yield~CRM + I(CRM^2),data = Rohwer.model.df.crm)

summary(Rohwer.model.crm)

coefficients(Rohwer.model.crm)

crmvalues <- seq(from= 100, to= 125, by= 1) 

Rohwer.yield.predict.crm <- predict(Rohwer.model.crm,list(CRM=crmvalues))

rm.Rohwer.plot <- 
  ggplot()+
  geom_point(data= Rohwer.df,
             aes(x= CRM, y= yield),fill= "grey",size= 1.2,alpha= .1)+
  geom_point(data=  Rohwer.model.df.crm,
             aes(x= CRM, y= yield), fill= 'blue',size= 2, shape= 21,alpha= .4)+
  geom_line(aes(x= crmvalues,y= Rohwer.yield.predict.crm),
            col= "red", linewidth= 1)+
  labs(x= "Corn Relative Maturity",y= expression(paste("Grain yield (",Bu~a^-1,")")))+
  theme_classic(base_size= 25)+
  ggtitle("Rohwer")+
  scale_y_continuous(limits = c(125,300))
ggsave(plot = rm.Rohwer.plot, "../output/rm.Rohwer.plot.png",
       width = 8, height = 7 , unit = "in", dpi = 600)

##Stuttgart.crm

Stuttgart.df2.crm <- Stuttgart.df %>% 
  left_join(rm.Stuttgart) %>% 
  mutate(abovep90 = ifelse(yield>p90,"yes","no")) %>% 
  filter(!is.na(CRM))

Stuttgart.model.df.crm <- Stuttgart.df2.crm %>% filter(abovep90== "yes")

Stuttgart.model.crm <- lm(yield~CRM + I(CRM^2),data = Stuttgart.model.df.crm)

summary(Stuttgart.model.crm)

coefficients(Stuttgart.model.crm)

crmvalues <- seq(from= 100, to= 125, by= 1) 

Stuttgart.yield.predict.crm <- predict(Stuttgart.model.crm,list(CRM=crmvalues))

rm.Stuttgart.plot <- 
  ggplot()+
  geom_point(data= Stuttgart.df,
             aes(x= CRM, y= yield),fill= "grey",size= 1.2,alpha= .1)+
  geom_point(data=  Stuttgart.model.df.crm,
             aes(x= CRM, y= yield), fill= 'blue',size=2, shape= 21,alpha= .4)+
  geom_line(aes(x= crmvalues,y= Stuttgart.yield.predict.crm),
            col= "red", linewidth= 1)+
  labs(x= "Corn Relative Maturity",y= expression(paste("Grain yield (",Bu~a^-1,")")))+
  theme_classic(base_size= 25)+
  ggtitle("Stuttgart")+
  scale_y_continuous(limits = c(125,300))
ggsave(plot = rm.Stuttgart.plot, "../output/rm.Stuttgart.plot.png",
       width = 8, height = 7 , unit = "in", dpi = 600)


perlocation.plot.crm <- plot_grid(rm.marianna.plot,rm.Rohwer.plot,rm.keiser.plot,rm.Stuttgart.plot)

final.plot.rm <- plot_grid(p16,perlocation.plot.crm)

ggsave(plot = final.plot.rm, "../output/final.plot.rm.tiff",
       width = 20, height = 10, unit = "in", dpi = 600)






### Grouping by CRM
#####

group_1 <- data %>% 
  filter(CRM %in% c(103,107,108) & !is.na(CRM)& !is.na(yield)) 


group_2 <- data %>% 
  filter(CRM %in% c(110,111,112) & !is.na(CRM)& !is.na(yield)) 

group_3 <- data %>% 
  filter(CRM %in% c(113,114,115) & !is.na(CRM)& !is.na(yield)) 

group_4 <- data %>% 
  filter(CRM %in% c(116,117,118) & !is.na(CRM)& !is.na(yield)) 

group_5 <- data %>% 
  filter(CRM %in% c(119,130,131) & !is.na(CRM)& !is.na(yield)) 

mean_group_1 <- group_1 %>%
  group_by(Year) %>%
  summarize(mean_yield = mean(yield, na.rm = TRUE),
            sd_yield= sd (yield,na.rm = TRUE),
                          .groups = 'drop') %>% 
  mutate(group = "103 to 108")

mean_group_2 <- group_2 %>%
  group_by(Year) %>%
  summarize(mean_yield = mean(yield, na.rm = TRUE), sd_yield = sd(yield, na.rm= TRUE),
            .groups = 'drop') %>% 
  mutate(group = "110 to 112")

mean_group_3 <- group_3 %>%
  group_by(Year) %>%
  summarize(mean_yield = mean(yield, na.rm = TRUE), sd_yield= sd(yield,na.rm= TRUE), 
            .groups = 'drop') %>% 
  mutate(group = "113 to 115")

mean_group_4 <- group_4 %>%
  group_by(Year) %>%
  summarize(mean_yield = mean(yield, na.rm = TRUE), sd_yield= sd(yield,na.rm= TRUE),
            .groups = 'drop') %>% 
  mutate(group = "116 to 118")

mean_group_5 <- group_5 %>%
  group_by(Year) %>%
  summarize(mean_yield = mean(yield, na.rm = TRUE), sd_yield= sd(yield,na.rm= TRUE),
            .groups = 'drop') %>% 
  mutate(group = "119 to 121")

combined_means <-  bind_rows(mean_group_1,mean_group_2,mean_group_3,mean_group_4,mean_group_5)
custom_colors <- c("red", "blue", "limegreen", "black", "yellow3")

crmgroups_yield <- 
ggplot(combined_means, aes(x = Year, y = mean_yield, color = group)) +
  geom_line(linewidth= 1) +               
  geom_point(size= 3) +              
  #geom_errorbar(aes(ymin = mean_yield - sd_yield, ymax = mean_yield + sd_yield), width = 0.2) + 
  scale_x_continuous(limits = c(min(data$Year), max(data$Year)), 
                     breaks = seq(min(data$Year), max(data$Year), by = 1)) + 
  scale_y_continuous(breaks = seq(6000,16000, by = 1000))+
  labs(title = "Mean Yield by Year and Relative Maturity Range",
       x = "Year",
       y = expression("Grain Yield (kg" ~ ha^{-1} ~ ")")) +
  scale_color_manual(values = custom_colors) +
  theme_classic(base_size= 25) +            
  theme(legend.title = element_blank())  
ggsave(plot = crmgroups_yield, "../output/crmgroups_yield.tiff",
       width = 20, height = 10, unit = "in", dpi = 600)








# View the updated data frame
print(crm.groups)

crm.groups %>%
  filter(Location == "Marianna") %>%
  ggplot(aes(x=Year , y= yield))  +
  geom_smooth(method = "lm",
              formula = y~x+I(x^2),se=F,
              aes(colour= crm_groups), na.rm = TRUE)+ 
  labs (x="Year", y = expression(paste("Grain yield (",Bu~a^-1,")"))) +
  scale_x_continuous(breaks = seq(min(data$Year), max(data$Year)))  



## Violin + boxplot
#####
  p17 <-  
  ggplot(data4,aes(x = Location, y = yield)) + 
  geom_violin(bw = 2.5, fill= 'steelblue') +
  # add a transparent boxplot and shrink its width to 0.3
  geom_boxplot(alpha=0, width=0.3) +
  # Reset point size to default and set point shape to 95
  geom_point(alpha = 0, shape = 95) +
  # Supply a subtitle detailing the kernel widthz
  labs(subtitle = 'yield distribution/ Location')+z
  theme_minimal()+
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  ggsave(plot = p17, "../output/Figure 17.png",
         width = 8, height = 7 , unit = "in", dpi = 600)
  
  
  p18 <-
    ggplot(data4, aes(x = Year, y = CRM)) + 
    geom_violin(bw = 2.5, fill= 'steelblue') +
    # add a transparent boxplot and shrink its width to 0.3
    geom_boxplot(alpha=0, width=0.3) +
    # Reset point size to default and set point shape to 95
    geom_point(alpha = 0, shape = 95) +
    # Supply a subtitle detailing the kernel width
    labs(subtitle = 'CRM distribution across Years') +
    theme_classic (base_size= 25)
  ggsave(plot = p18, "../output/Figure 18.png",
         width = 8, height = 7 , unit = "in", dpi = 600)
  
  p19 <-
    ggplot(data4,aes(x= CRM))+
    geom_histogram(fill= 'steelblue')+ 
    labs(title = "CRM distribution")+
    theme_classic(base_size = 25)
  ggsave(plot = p19,"../output/Figure 19.png", 
         width= 8, height= 7, unit= "in", dpi= 600)
  
  ggplot(data4,aes(x= yield))+
    geom_boxplot()
summary(data4)
#####

# CORRELATION MATRIX

data_clean <- data %>%
  filter(
    !is.na(yield) & !is.na(CRM) & !is.na(Year) & !is.na(SowingdateDOY) & 
      !is.na(Location) & !is.na(Hybrid) & !is.na(soil_textural_class)
  ) %>%
  filter(Location != "Greenfield" & Location != "Harrisburg")

# Standardizing Numeric Variables
data_clean$mean_maxt <- scale(data_clean$mean_maxt)
data_clean$row_spacing <- scale(data_clean$row_spacing)

data_clean$soil_textural_class_numeric <- as.numeric(as.factor(data_clean$soil_textural_class))


data_clean_all <- data_clean %>%
  select(yield, CRM, Year, SowingdateDOY, mean_maxt, mean_mint, cum_radn, 
         april, may, june, july, august, sept, Precipitation, 
         soil_textural_class_numeric, row_spacing, above_32) %>%
  select_if(is.numeric)  

# Remove Variables with Zero Variance
data_clean_all <- data_clean_all %>%
  select_if(function(x) sd(x, na.rm = TRUE) > 0)  

# Compute Correlation Matrix
correlation_matrix_all <- cor(data_clean_all, use = "complete.obs")
yield_correlations_all <- correlation_matrix_all["yield", ]
yield_correlations_all <- yield_correlations_all[-which(names(yield_correlations_all) == "yield")]


correlation_df <- data.frame(
  Variable = names(yield_correlations_all),
  Correlation = yield_correlations_all
)


ggplot(correlation_df, aes(x = Correlation, y = Variable, fill = Correlation > 0)) +
  geom_bar(stat = "identity") +  
  scale_fill_manual(values = c("lightcoral", "lightblue")) +  
  labs(title = "Correlation with Yield (All Locations Combined)",
       x = "Correlation", y = "Variable") +
  theme_classic() +
  theme(axis.text.x = element_text(size= 14),
        axis.text.y = element_text(size= 14),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size= 16),
        legend.position = "none")
ggsave("All_locations_correlation_plot.tiff",
       width = 25, height = 20, units = "cm",
       dpi = 600, compression = "lzw", bg = "white")
##Keiser
data_clean <- data %>%
  filter(
    !is.na(yield) & !is.na(CRM) & !is.na(Year) & !is.na(SowingdateDOY) & 
      !is.na(Location) & !is.na(Hybrid)
  ) %>%
  
  filter(Location != "Greenfield" & Location != "Harrisburg")
data_clean_keiser <- data_clean %>%
  filter(Location== "Keiser")%>%
  select(yield, CRM, Year, SowingdateDOY, mean_maxt, mean_mint, cum_radn, 
         april, may, june, july, august, sept, Precipitation, soil_textural_class, 
         row_spacing, above_32) %>%
  select_if(is.numeric)  

data_clean_keiser <- data_clean_keiser %>%
  select_if(function(x) sd(x, na.rm = TRUE) > 0)  

correlation_matrix_keiser <- cor(data_clean_keiser, use = "complete.obs")
yield_correlations_keiser <- correlation_matrix_keiser["yield", ]
yield_correlations_keiser <- yield_correlations_keiser[-which(names(yield_correlations_keiser) == "yield")]

# Convert to a data frame for ggplot
correlation_df.keiser <- data.frame(
  Variable = names(yield_correlations_keiser),
  Correlation = yield_correlations_keiser
)


ggplot(correlation_df.keiser, aes(x = Correlation, y = Variable, fill = Correlation > 0)) +
  geom_bar(stat = "identity") +  
  scale_fill_manual(values = c("lightcoral", "lightblue")) +  
  labs(title = "Keiser",
       x = "Correlation", y = "Variable") +
  theme_classic() +
  theme(axis.text.y = element_text(angle = 0),  
        legend.position = "none")  

ggsave("Keiser correlation.tiff",
       width=25,height=20,units ="cm",
       dpi=600,compression="lzw",bg="white")

##Marianna
data_clean <- data %>%
  filter(
    !is.na(yield) & !is.na(CRM) & !is.na(Year) & !is.na(SowingdateDOY) & 
      !is.na(Location) & !is.na(Hybrid)
  ) %>%
  
  filter(Location != "Greenfield" & Location != "Harrisburg")
data_clean_marianna <- data_clean %>%
  filter(Location== "Marianna")%>%
  select(yield, CRM, Year, SowingdateDOY, mean_maxt, mean_mint, cum_radn, 
         april, may, june, july, august, sept, Precipitation, soil_textural_class, 
         row_spacing, above_32) %>%
  select_if(is.numeric)  

data_clean_marianna <- data_clean_marianna %>%
  select_if(function(x) sd(x, na.rm = TRUE) > 0)  

correlation_matrix_marianna <- cor(data_clean_marianna, use = "complete.obs")
yield_correlations_marianna <- correlation_matrix_marianna["yield", ]
yield_correlations_marianna <- yield_correlations_marianna[-which(names(yield_correlations_marianna) == "yield")]

# Convert to a data frame for ggplot
correlation_df.marianna <- data.frame(
  Variable = names(yield_correlations_marianna),
  Correlation = yield_correlations_marianna
)


ggplot(correlation_df.marianna, aes(x = Correlation, y = Variable, fill = Correlation > 0)) +
  geom_bar(stat = "identity") +  
  scale_fill_manual(values = c("lightcoral", "lightblue")) +  
  labs(title = "Marianna",
       x = "Correlation", y = "Variable") +
  theme_classic() +
  theme(axis.text.y = element_text(angle = 0),  
        legend.position = "none")  

ggsave("marianna correlation.tiff",
       width=25,height=20,units ="cm",
       dpi=600,compression="lzw",bg="white")


##Stuttgart
data_clean <- data %>%
  filter(
    !is.na(yield) & !is.na(CRM) & !is.na(Year) & !is.na(SowingdateDOY) & 
      !is.na(Location) & !is.na(Hybrid)
  ) %>%
  
  filter(Location != "Greenfield" & Location != "Harrisburg")
data_clean_stuttgart <- data_clean %>%
  filter(Location== "Stuttgart")%>%
  select(yield, CRM, Year, SowingdateDOY, mean_maxt, mean_mint, cum_radn, 
         april, may, june, july, august, sept, Precipitation, soil_textural_class, 
         row_spacing, above_32) %>%
  select_if(is.numeric)  

data_clean_stuttgart <- data_clean_stuttgart %>%
  select_if(function(x) sd(x, na.rm = TRUE) > 0)  

correlation_matrix_stuttgart <- cor(data_clean_stuttgart, use = "complete.obs")
yield_correlations_stuttgart <- correlation_matrix_stuttgart["yield", ]
yield_correlations_stuttgart <- yield_correlations_stuttgart[-which(names(yield_correlations_stuttgart) == "yield")]

# Convert to a data frame for ggplot
correlation_df.stuttgart <- data.frame(
  Variable = names(yield_correlations_stuttgart),
  Correlation = yield_correlations_stuttgart
)


ggplot(correlation_df.stuttgart, aes(x = Correlation, y = Variable, fill = Correlation > 0)) +
  geom_bar(stat = "identity") +  
  scale_fill_manual(values = c("lightcoral", "lightblue")) +  
  labs(title = "Stuttgart",
       x = "Correlation", y = "Variable") +
  theme_classic() +
  theme(axis.text.y = element_text(angle = 0),  
        legend.position = "none")  

ggsave("stuttgart correlation.tiff",
       width=25,height=20,units ="cm",
       dpi=600,compression="lzw",bg="white")

##Rohwer
data_clean <- data %>%
  filter(
    !is.na(yield) & !is.na(CRM) & !is.na(Year) & !is.na(SowingdateDOY) & 
      !is.na(Location) & !is.na(Hybrid)
  ) %>%
  
  filter(Location != "Greenfield" & Location != "Harrisburg")
data_clean_rohwer <- data_clean %>%
  filter(Location== "Rohwer")%>%
  select(yield, CRM, Year, SowingdateDOY, mean_maxt, mean_mint, cum_radn, 
         april, may, june, july, august, sept, Precipitation, soil_textural_class, 
         row_spacing, above_32) %>%
  select_if(is.numeric)  

data_clean_rohwer <- data_clean_rohwer %>%
  select_if(function(x) sd(x, na.rm = TRUE) > 0)  

correlation_matrix_rohwer <- cor(data_clean_rohwer, use = "complete.obs")
yield_correlations_rohwer <- correlation_matrix_rohwer["yield", ]
yield_correlations_rohwer <- yield_correlations_rohwer[-which(names(yield_correlations_rohwer) == "yield")]

# Convert to a data frame for ggplot
correlation_df.rohwer <- data.frame(
  Variable = names(yield_correlations_rohwer),
  Correlation = yield_correlations_rohwer
)


ggplot(correlation_df.rohwer, aes(x = Correlation, y = Variable, fill = Correlation > 0)) +
  geom_bar(stat = "identity") +  
  scale_fill_manual(values = c("lightcoral", "lightblue")) +  
  labs(title = "Rohwer",
       x = "Correlation", y = "Variable") +
  theme_classic() +
  theme(axis.text.y = element_text(angle = 0),  
        legend.position = "none")  

ggsave("rohwer correlation.tiff",
       width=25,height=20,units ="cm",
       dpi=600,compression="lzw",bg="white")

##Bellfarm
data_clean <- data %>%
  filter(
    !is.na(yield) & !is.na(CRM) & !is.na(Year) & !is.na(SowingdateDOY) & 
      !is.na(Location) & !is.na(Hybrid)
  ) %>%
  
  filter(Location != "Greenfield" & Location != "Harrisburg")
data_clean_bellfarm <- data_clean %>%
  filter(Location== "Bellfarm")%>%
  select(yield, CRM, Year, SowingdateDOY, mean_maxt, mean_mint, cum_radn, 
         april, may, june, july, august, sept, Precipitation, soil_textural_class, 
         row_spacing, above_32) %>%
  select_if(is.numeric)  

data_clean_bellfarm <- data_clean_bellfarm %>%
  select_if(function(x) sd(x, na.rm = TRUE) > 0)  

correlation_matrix_bellfarm <- cor(data_clean_bellfarm, use = "complete.obs")
yield_correlations_bellfarm <- correlation_matrix_bellfarm["yield", ]
yield_correlations_bellfarm <- yield_correlations_bellfarm[-which(names(yield_correlations_bellfarm) == "yield")]

# Convert to a data frame for ggplot
correlation_df.bellfarm <- data.frame(
  Variable = names(yield_correlations_bellfarm),
  Correlation = yield_correlations_bellfarm
)


ggplot(correlation_df.bellfarm, aes(x = Correlation, y = Variable, fill = Correlation > 0)) +
  geom_bar(stat = "identity") +  
  scale_fill_manual(values = c("lightcoral", "lightblue")) +  
  labs(title = "Bellfarm",
       x = "Correlation", y = "Variable") +
  theme_classic() +
  theme(axis.text.y = element_text(angle = 0),  
        legend.position = "none")  

ggsave("bellfarm correlation.tiff",
       width=25,height=20,units ="cm",
       dpi=600,compression="lzw",bg="white")



##model 1
model_clean1 <- lmer(yield ~ CRM + Year + SowingdateDOY + Location + mean_maxt + mean_mint + cum_radn + cum_rain + above_32 +(1 | Hybrid), data = data_clean)
summary(model_clean1)
anova(model_clean1)
AIC(model_clean1)

##emmeans

library(emmeans)

data_clean$SowingdateDOY <- as.factor(data_clean$SowingdateDOY)

levels(data_clean$SowingdateDOY)

model.a <- lmer(yield ~ SowingdateDOY + (1 | Hybrid) + (1 | Year)+ (1 | Location), data = data_clean)

AIC(model.a)

emmeans(model.a,"SowingdateDOY")
pairs(emmeans(model.a,"SowingdateDOY"))


##model 2
model_clean2 <- lmer(yield ~ CRM * Location + Year + SowingdateDOY +
                       scale(mean_maxt) +
                       scale(mean_mint) + 
                       log(cum_radn +1) + 
                       log(cum_rain + 1) + 
                       above_32 +
                       (1 | Hybrid), data = data_clean)
summary(model_clean2)

AIC(model_clean2)
BIC(model_clean2)
residuals2 <- resid(model_clean2)
summary(residuals2)

qqnorm(residuals2)
qqline(residuals2, col = "red")


ggplot(data_clean,aes(x=cum_rain))+
  geom_histogram()

ggplot(data_clean,aes(x= log(cum_radn)))+
  geom_histogram()


model_clean2.6 <- lmer(yield ~ CRM + SowingdateDOY + row_spacing +
                         soil_textural_class +
                         mean_maxt +
                         mean_mint +
                         cum_radn +
                         Precipitation +
                         above_32 +
                         (1 | Hybrid) + (1|Year),
                       data = data_clean)
summary(model_clean2.6)
AIC(model_clean2.6)
BIC(model_clean2.6)
vif(model_clean2.6)

model_clean2.7 <- lmer(yield ~ CRM + Year + SowingdateDOY + row_spacing +
                         soil_textural_class +
                         Precipitation +
                         mean_maxt +
                         mean_mint +
                         cum_radn +
                         above_32 +
                         (1 | Hybrid) + (1|Location),
                       data = data_clean)
summary(model_clean2.7)
AIC(model_clean2.7)
BIC(model_clean2.7)
vif(model_clean2.7)


model_clean2.8 <- lmer(yield ~ CRM + Location * Year + SowingdateDOY + row_spacing +
                         soil_textural_class +
                         Precipitation +
                         mean_maxt +
                         mean_mint +
                         cum_radn +
                         above_32 +
                         (1 | Hybrid),
                       data = data_clean)
summary(model_clean2.8)
AIC(model_clean2.8)
BIC(model_clean2.8)
vif(model_clean2.8)


model_clean3 <- lmer(yield ~ CRM * SowingdateDOY + Location + Year + row_spacing +
                         soil_textural_class +
                         Precipitation +
                         mean_maxt +
                         mean_mint +
                         cum_radn +
                         above_32 +
                         (1 | Hybrid),
                       data = data_clean)
summary(model_clean3)
AIC(model_clean3)
BIC(model_clean3)
vif(model_clean3)

model_clean3.1 <- lmer(yield ~ CRM * SowingdateDOY + Location + Year + row_spacing +
                       soil_textural_class +
                       (1 | Precipitation) +
                       (1 | mean_maxt) +
                       (1 |mean_mint) +
                       (1 |cum_radn) +
                       (1 | above_32) +
                       (1 | Hybrid),
                     data = data_clean)
summary(model_clean3.1)
AIC(model_clean3.1)
BIC(model_clean3.1)
vif(model_clean3.1)

model_clean3.2 <- lmer(relative_yield ~ CRM + SowingdateDOY + Location + Year + row_spacing +
                         soil_textural_class +
                         (1 | Precipitation) +
                         (1 | mean_maxt) +
                         (1 |mean_mint) +
                         (1 |cum_radn) +
                         (1 | above_32) +
                         (1 | Hybrid),
                       data = data_clean)
summary(model_clean3.2)
AIC(model_clean3.2)
BIC(model_clean3.2)
vif(model_clean3.2)

str(data_clean)
summary(data_clean$SowingdateDOY)
data_clean$SowingdateDOY_std <- scale(data_clean$SowingdateDOY)

model_clean3.3 <- lmer(yield ~ CRM + SowingdateDOY_std + Location + Year + row_spacing +
                         soil_textural_class +
                         (1 | Precipitation) +
                         (1 | mean_maxt) +
                         (1 |mean_mint) +
                         (1 |cum_radn) +
                         (1 | above_32) +
                         (1 | Hybrid),
                       data = data_clean)
summary(model_clean3.3)
AIC(model_clean3.3)
BIC(model_clean3.3)
vif(model_clean3.3)

model_clean3.4 <- lmer(yield ~ CRM + SowingdateDOY_std + Location + Year + row_spacing +
                         soil_textural_class +
                         (1 | Precipitation) +
                         (1 | mean_maxt) +
                         (1 | Hybrid),
                       data = data_clean)
summary(model_clean3.4)
AIC(model_clean3.4)
BIC(model_clean3.4)
vif(model_clean3.4)

data_clean1 <- data_clean %>%
  group_by(Location, CRM, Year) %>%
  mutate(max_yield_per_crm = max(yield, na.rm = TRUE)) %>% 
  ungroup() 
data_clean1 <- data_clean1 %>% 
  mutate(relative_yield = yield / max_yield_per_crm)

data_clean1$SowingdateDOY <- scale(data_clean1$SowingdateDOY)
data_clean1$row_spacing <- scale(data_clean1$row_spacing)
data_clean1$CRM <- scale(data_clean1$CRM)
data_clean1$Year <- as.factor(data_clean1$Year)


model_final <- lmer(relative_yield ~ CRM * SowingdateDOY + Location * SowingdateDOY + 
                      Year * SowingdateDOY + CRM * Location + row_spacing + 
                      (1 | Hybrid), data = data_clean1)

summary(model_final)

AIC(model_final)
BIC(model_final)
vif(model_final)
residuals.model_final <- residuals(model_final)
qqnorm(residuals.model_final)
qqline(residuals.model_final)

cor(data_clean[,c("CRM","SowingdateDOY","row_spacing", "mean_maxt","cum_radn")])

model_clean3.2 <- lmer(relative_yield ~ CRM + SowingdateDOY + Location + Year + row_spacing +
                         soil_textural_class +
                         (1 | Precipitation) +
                         (1 | mean_maxt) +
                         (1 |mean_mint) +
                         (1 |cum_radn) +
                         (1 | above_32) +
                         (1 | Hybrid),
                       data = data_clean1)
summary(model_clean3.2)
AIC(model_clean3.2)
BIC(model_clean3.2)
vif(model_clean3.2)


##Elli, 1/15/25
write.csv(data_clean,"data_clean.csv")