## R code for "Avian colouration in polluted world: a meta-analysis" 

setwd("C:/Users/katar/Documents/Meta-analiza")

library(ggplot2)
library(ggpubr)
library(packcircles)
library(metafor)
library(devtools)
library(orchaRd)
library(patchwork)
library(tidyverse)
library(metafor)
library(ape)
library(metaAidR)

install.packages("pacman")
pacman::p_load(devtools, tidyverse, metafor, patchwork, R.rsp, emmeans)
devtools::install_github("daniel1noble/orchaRd", force = TRUE)
vignette("orchaRd")

# load data
mydata <- read.csv(here("C:/Users/katar/Documents/Meta-analiza/ES_biol_2023.csv"), sep = ";")

# load phylogenetic tree
mytree <- multi2di(read.tree("bird_tree.tre"))
mytree$edge.length[which(mytree$edge.length == 0)] = 1e-10

#subset tree to the species list
subtree <- drop.tip(mytree, setdiff(mytree$tip.label, mydata$tip_label))
cor <- vcv(subtree, cor = T)

#calculate weights
mydata$Z_var <- 1/(mydata$n-3)
mydata$weight <- 1/mydata$Z_var
summary(mydata)

# Model 0.0
M0.0 <- rma.mv(yi = ES_biol, V = Z_var,
                        method = "REML",
                        mods = ~ 1,
                        random = list(~1|tip_label,
                                      ~1|research_group,
                                      ~1|Study_no,
                                      ~1|es_ID),
                        R = list(tip_label = cor),
                        control=list(optimizer="optim",optmethod="Nelder-Mead"),
                        data = mydata)
summary(M0.0)

# Model 0.1
M0.1 <- rma.mv(yi = ES_biol, V = Z_var,
               method = "REML",
               mods = ~ 1,
               random = list(~1|tip_label,
                             ~1|research_group,
                             ~1|Study_no,
                             ~1|es_ID),
               control=list(optimizer="optim",optmethod="Nelder-Mead"),
               data = mydata)
summary(M0.1)
i2_ml(M0.1)

# Model 1.0
M1.0<- rma.mv(yi = ES_biol, V = Z_var,
              method = "REML",
              mods = ~ colour_type-1+sex+age+measured_part+study_type+raw_or_model+method ,
              random = list(~1|tip_label,
                            ~1|research_group,
                            ~1|Study_no,
                            ~1|es_ID),
              control=list(optimizer="optim",optmethod="Nelder-Mead"),
              data = mydata)
summary(M1.0)

# Figure 2. Orchard plot of the impact of pollution on avian colour traits as a function of the colour-producing mechanism conditioned for the measurement method (Nakagawa et al. 2023).

orchard_plot(M1.0, mod = "colour_type", group = "Study_no", by="method", xlab = "Correlation coefficient", transfm = "none", alpha = 0.4, angle = 0,  cb = FALSE, legend.pos = "top.left", condition.lab = "Method") + theme(legend.direction = "vertical")  
#extraction of marginal estimates for colour type from the M1.0:
res_colour_type <- orchaRd::mod_results(M1.0, mod = "colour_type",  group = "Study_no") 
res_colour_type

#ad. Figure S2
orchard_plot(M1.0, mod = "sex", group = "Study_no", xlab = "Correlation coefficient (r)", angle = 0,legend.pos = "top.left",  cb = FALSE) + theme(legend.direction = "vertical") 
res_sex <- orchaRd::mod_results(M1.0, mod = "sex",  group = "Study_no")
res_sex

orchard_plot(M1.0, mod = "age", group = "Study_no", xlab = "Correlation coefficient (r)", angle = 0,legend.pos = "top.left",  cb = FALSE) + theme(legend.direction = "vertical") 
res_age <- orchaRd::mod_results(M1.0, mod = "age",  group = "Study_no")
res_age

orchard_plot(M1.0, mod = "measured_part", group = "Study_no", xlab = "Correlation coefficient (r)", angle = 0,legend.pos = "top.left",  cb = FALSE) + theme(legend.direction = "vertical") 
res_measured_part <- orchaRd::mod_results(M1.0, mod = "measured_part",  group = "Study_no")
res_measured_part

#ad. Figure S3
orchard_plot(M1.0, mod = "raw_or_model", group = "Study_no", xlab = "Correlation coefficient (r)", angle = 0,legend.pos = "top.left",  cb = FALSE) + theme(legend.direction = "vertical") 
res_raw_or_model <- orchaRd::mod_results(M1.0, mod = "raw_or_model",  group = "Study_no")
res_raw_or_model

orchard_plot(M1.0, mod = "study_type", group = "Study_no", xlab = "Correlation coefficient (r)", angle = 0,legend.pos = "top.left",  cb = TRUE) + theme(legend.direction = "vertical") 
res_study_type <- orchaRd::mod_results(M1.0, mod = "study_type",  group = "Study_no")
res_study_type

orchard_plot(M1.0, mod = "method", group = "Study_no", xlab = "Correlation coefficient (r)", angle = 0,legend.pos = "top.left",  cb = TRUE) + theme(legend.direction = "vertical") 
res_method <- orchaRd::mod_results(M1.0, mod = "method",  group = "Study_no")
res_method

# Models 2.1 - 2.4

# Model 2.1 (on the subset of 'hue')
hue <- subset(mydata, HSB == "hue")
summary(hue)
mydata$Z_var <- 1/(mydata$n-3)
mydata$weight <- 1/mydata$Z_var
M2.1 <- rma.mv(yi = ES_biol, V = Z_var,
                    method = "REML",
                    mods = ~colour_type+pollution_type+colour_type:pollution_type+
                    sex+age+measured_part+study_type +raw_or_model+ method,
                    random = list(~1|tip_label,
                                  ~1|Study_no,
                                  ~1|es_ID),
                    control=list(optimizer="optim",optmethod="Nelder-Mead"),
                    data = hue)
summary(M2.1)# interaction not significant

M2.1 <- rma.mv(yi = ES_biol, V = Z_var,
                    method = "REML",
                    mods = ~colour_type+pollution_type+sex+age+measured_part
                    + study_type +raw_or_model+ method,
                    random = list(~1|tip_label,
                                  ~1|Study_no,
                                  ~1|es_ID),
                    control=list(optimizer="optim",optmethod="Nelder-Mead"),
                    data = hue)

summary(M2.1)

# Model 2.2 (on the subset of 'saturation')
saturation <- subset(mydata, HSB == "saturation")
summary(saturation)
M2.2 <- rma.mv(yi = ES_biol, V = Z_var,
                           method = "REML",
                           mods = ~colour_type+pollution_type+colour_type:pollution_type
                           +sex+age+measured_part+ study_type +raw_or_model+ method,
                           random = list(~1|tip_label,
                                         ~1|Study_no,
                                         ~1|es_ID),
                           control=list(optimizer="optim",optmethod="Nelder-Mead"),
                           data = saturation)

summary(M2.2)# interaction not significant
M2.2 <- rma.mv(yi = ES_biol, V = Z_var,
                           method = "REML",
                           mods = ~colour_type+pollution_type+sex+age+measured_part
                           +study_type +raw_or_model+ method,
                           random = list(~1|tip_label,
                                         ~1|Study_no,
                                         ~1|es_ID),
                           control=list(optimizer="optim",optmethod="Nelder-Mead"),
                           data = saturation)
summary(M2.2)

# Model 2.3 (on the subset of 'brightness')
brightness <- subset(mydata, HSB == "brightness")
summary(brightness)
M2.3<- rma.mv(yi = ES_biol, V = Z_var,
                          method = "REML",
                          mods = ~colour_type+pollution_type+colour_type:pollution_type+sex+age+measured_part+ study_type +raw_or_model+ method,
                          random = list(~1|tip_label,
                                        ~1|Study_no,
                                        ~1|es_ID),
                          control=list(optimizer="optim",optmethod="Nelder-Mead"),
                          data = brightness)
summary(M2.3)# interaction significant

# Model 2.4 (on the subset of 'other')
other <- subset(mydata, HSB == "other")
summary(other)
M2.4<- rma.mv(yi = ES_biol, V = Z_var,
                     method = "REML",
                     mods = ~colour_type+pollution_type+colour_type:pollution_type
                     +sex+age+measured_part+ study_type +raw_or_model+ method,
                     random = list(~1|tip_label,
                                   ~1|Study_no,
                                   ~1|es_ID),
                     control=list(optimizer="optim",optmethod="Nelder-Mead"),
                     data = other)
summary(M2.4)# interaction not significant

M2.4<- rma.mv(yi = ES_biol, V = Z_var,
                     method = "REML",
                     mods = ~colour_type+pollution_type+sex+age+measured_part+study_type
                     +raw_or_model+ method,
                     random = list(~1|tip_label,
                                   ~1|Study_no,
                                   ~1|es_ID),
                     control=list(optimizer="optim",optmethod="Nelder-Mead"),
                     data = other)
summary(M2.4)

# Model 3.0 (on the subset of carotenoid-based colouration)
carotenoids <- subset(mydata, colour_type == "carotenoid-based")
M3.0 <- rma.mv(yi = ES_biol, V = Z_var,
                           method = "REML",
                           mods = ~ pollution_type + measured_part:carotenoids_type+HSB+sex+age +measured_part
                           +study_type + raw_or_model+carotenoids_type,
                           random = list(~1|tip_label,
                                         ~1|Study_no,
                                         ~1|es_ID),
                           control=list(optimizer="optim",optmethod="Nelder-Mead"),
                           data = carotenoids)
summary(M3.0) # interaction not significant

M3.0 <- rma.mv(yi = ES_biol, V = Z_var,
                           method = "REML",
                           mods = ~ pollution_type +HSB+sex+age +measured_part
                           +study_type + raw_or_model+carotenoids_type,
                           random = list(~1|tip_label,
                                         ~1|Study_no,
                                         ~1|es_ID),
                           control=list(optimizer="optim",optmethod="Nelder-Mead"),
                           data = carotenoids)
summary(M3.0)

#Model 4.0 - Publication bias

M4.0 <- run.model(mydata, ~sqrt(Z_var))
summary(M4.0)


# Model 5.0 Time lag bias

M5.0 <- run.model(mydata, ~date)
summary(M5.0)


# Function to create Figure 3 and Figure S4

run.model<-function(mydata,formula){
  
  mydata<-as.data.frame(mydata)
  Z<-make_VCV_matrix(mydata, V="Z_var", cluster = "Study_no", obs="es_ID")
  
  rma.mv(yi=ES_biol, V=Z,
         method="REML",
         mods= formula,
         random = list(~1|tip_label,
                       ~1|Study_no,
                       ~1|es_ID,
                       ~1|research_group),
         control=list(optimizer="optim",optmethod="Nelder-Mead"),
         data = mydata)
}

plot_continuous<-function(data, model, moderator, xlab){
  pred<-predict.rma(model)
  data %>% mutate(fit=pred$pred, 
                  ci.lb=pred$ci.lb,
                  ci.ub=pred$ci.ub,
                  pr.lb=pred$cr.lb,
                  pr.ub=pred$cr.ub) %>% 
    ggplot(aes(x = moderator, y = ES_biol)) +
    geom_ribbon(aes(ymin = pr.lb, ymax = pr.ub, color = NULL), alpha = .1) +
    geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub, color = NULL), alpha = .3) +
    geom_point(aes(size=(1/sqrt(mydata$Z_var))), shape= 21, alpha=0.7, fill="#6F94B7", col="#305980") +
    geom_line(aes(y = fit), size = 1.5)+  
    labs(x = xlab, y = "Effect size", size = "Precision (1/SE)") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_size_continuous(range=c(1,12))+
    geom_hline(yintercept = 0,linetype = 2, colour = "black",alpha=0.5)+   
    theme(text = element_text(size = 18, colour = "black", hjust = 0.5),
          legend.text=element_text(size=14),
          legend.position=c(0,0), 
          legend.justification = c(0,0),
          legend.background = element_blank(), 
          legend.direction="horizontal",
          legend.title = element_text(size=15), 
          panel.border=element_rect(colour="black", fill=NA, size=1.2))
}

# Figure 3. Bubble plot showing a positive significant temporal trend in trend in published effect size estimates
pred_year_model <- predict.rma(M5.0)
plot_continuous(mydata, year_model, mydata$date, "Publication year")

# Figure S4. Funnel plot of effect sizes from the biology ES data set versus standard errors of effect sizes
pred_Egger_reg <- predict.rma(M4.0)  
plot_continuous(mydata, M4.0, sqrt(mydata$Z_var), "Standard error")