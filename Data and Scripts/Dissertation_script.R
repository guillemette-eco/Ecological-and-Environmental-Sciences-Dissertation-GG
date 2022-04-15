# Guillemette Gandon                 
# School of GeoSciences, University of Edinburgh  
# Ecological and Environmental Sciences Dissertation:                    
# A DECADE OF CHANGE IN PLANT COMMUNITIES OF THE GREAT SMOKY MOUNTAINS, USA

# Libraries

library(tidyverse)
library(png)
library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)
library(vegan)
library(stargazer)
library(lme4)
library(rstan)
library(gdata)
library(bayesplot)
library(rstanarm)
library(RColorBrewer)
library(stats)
library(blmeco) 
library(pscl)
library(boot)
library(MuMIn)
library(DALEX)
library(ggeffects) 

## Load dataframe
setwd("-")
herb <- read.csv("herb_raw_data.csv")

# Tidying and filtering dataset ----

## For cover, richness and diversity 

Herb_new <- dplyr::select(herb, -c(26:109))
Herb_new <- dplyr::select(Herb_new, -c(1,2,3,5,7,8,10,11,12,13,26,27))
Herb_long <- gather(Herb_new, Year, Cover, 5:15) # long format
Herb_long$Cover[Herb_long$Cover == 999] <- "NA" 
Herb_long$Year <- parse_number(Herb_long$Year)
Herb_long <- filter(Herb_long, Cover != "NA") # do not take into account missing values (see method)
Herb_long$Cover <- as.numeric(Herb_long$Cover)
Herb_long <- filter(Herb_long, BroaderGroup != "ZZbare")
Herb_long <- filter(Herb_long, BroaderGroup != "woody")
Herb_long <- filter(Herb_long, Taxa4Analysis != "Zunknown")
Herb_long <- filter(Herb_long, Year != "2007") 
str(Herb_long)
Herb_long$SapPlot <- as.factor(Herb_long$SapPlot)
herb_cover <- Herb_long
str(herb_cover)

## For NMDS matrix

herb_co <- dplyr::select(herb, -c(26:111))
herb_co <- dplyr::select(herb_co, -c(15))
herb_co <- gather(herb_co, Year, Cover, 15:24)
herb_co$Year <- parse_number(herb_co$Year)

herb_co <- filter(herb_co, BroaderGroup != "ZZbare")
herb_co <- filter(herb_co, BroaderGroup != "woody")
herb_co <- filter(herb_co, Taxa4Analysis != "Zunknown")

herb_co <- herb_co %>% dplyr::select(-c("AllYrs","CountInds","Sort","corner","Northing","Easting","BroaderGroup","Genus","Species"))
herb_co[is.na(herb_co)] <- 0
str(herb_co)

# Make matrix for NMDS

herb_matrix <- aggregate( Cover ~ Year + Taxa4Analysis + PlotID + SapPlot, herb_co, mean )

herb_matrix <- spread(herb_matrix, Taxa4Analysis, Cover)

str(herb_matrix)
herb_matrix$Year <- as.factor(herb_matrix$Year)
herb_matrix$SapPlot <- as.factor(herb_matrix$SapPlot)

# Community composition assemblages: NMDS ----

herb.NMDS <- metaMDS(herb_matrix [,4:41], distance = "bray", k = 8, trymax=300)

par(mfrow=c(1,1))
herb.NMDS$stress 
plot(herb.NMDS)

## Plotting community assemblages across years

group <- as.character(herb_matrix$Year)
colours <- as.character(herb_matrix$Year)

as.character(herb_matrix$Year) %>% 
  replace(colours=="2010", "#ddffa2") %>% 
  replace(colours=="2011", "#a9c3e1") %>% 
  replace(colours=="2012", "#ff4f4f") %>% 
  replace(colours=="2013", "#7CCD7C") %>% 
  replace(colours=="2014", "#CD69C9") %>% 
  replace(colours=="2015", "#104E8B") %>%
  replace(colours=="2016", "#EEB4B4") %>% 
  replace(colours=="2017", "#CD69C9") %>% 
  replace(colours=="2018", "#104E8B") -> colours

par(mfrow=c(1,1))
ordiplot(herb.NMDS, type = "n", cex.axis = 2, cex.lab=2)

for(i in unique(group)) {
  ordihull(herb.NMDS$point[grep(i, group),], draw="polygon",
           groups = group[group == i],col = colours[grep(i,group)],label=F) } 

#orditorp(herb.NMDS, display = "species", col = "black", air = 0.01)
orditorp(herb.NMDS, display = "sites", label=F, air = 0.01, cex = 1.25)

## Plotting community assemblages across plots

group <- as.character(herb_matrix$SapPlot)
colours <- as.character(herb_matrix$SapPlot)

as.character(herb_matrix$SapPlot) %>% 
  replace(colours=="1", "#ddffa2") %>% 
  replace(colours=="2", "#a9c3e1") %>% 
  replace(colours=="3", "#ff4f4f") %>% 
  replace(colours=="4", "#7CCD7C") %>% 
  replace(colours=="5", "#CD69C9") %>% 
  replace(colours=="6", "#104E8B") %>%
  replace(colours=="7", "#EEB4B4") -> colours

par(mfrow=c(1,1))
ordiplot(herb.NMDS, type = "n", cex.axis = 1.4, cex.lab=1.4)

for(i in unique(group)) {
  ordihull(herb.NMDS$point[grep(i, group),], draw="polygon",
           groups = group[group == i],col = colours[grep(i,group)],label=F) } 

orditorp(herb.NMDS, display = "species", col = "black", air = 0.01)
orditorp(herb.NMDS, display = "sites", label=F, air = 0.01, pch="+", cex = 0.4)

## Statistics: PERMANOVA ----

#multivariate analysis of variance used to compare groups of objects

# testing for the effect of years 
(herb.fit <- adonis(herb_matrix[,4:43] ~ Year, herb_matrix, 
                    permutations = 999, method = "bray"))

# testing for the effect of location
(herb.fit <- adonis(herb_matrix[,4:43] ~ SapPlot, herb_matrix, 
                    permutations = 999, method = "bray"))

# Cover Statistics ----

## Temporal changes ----

### Distribution

hist(herb_cover$Cover, col = "#F0E68C", breaks=100, xlim = c(0, 40))

#### Null model

null.cover <- glm.nb(Cover ~ 1, data = herb_cover)
AICc(null.cover)

### Fitting glmer nb model

cover.glmer <- glmer.nb(Cover ~ Year + (1|PlotID) + (1 + Year|Taxa4Analysis), data = herb_cover)
summary(cover.glmer)
plot(cover.glmer)

AICc(cover.glmer) # calculating AICc
dispersion_glmer(cover.glmer) # checking for overdispersion
r.squaredGLMM(cover.glmer) # extracting marginal and conditional R squares

stargazer(cover.glmer, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

#### Plotting results

pred.mm <- ggpredict(cover.glmer, terms = c("Year"))

(ggplot(pred.mm) + 
    geom_line(aes(x = x, y = predicted), col="black") +          # slope
    geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                fill = "#9BCD9B", alpha = 0.5) +  # error band
    geom_point(data = herb_cover,                      # adding the raw data (scaled values)
               aes(x = as.numeric(Year), y = Cover), col = "#9BCD9B") + 
    labs(x = "Year", y = "Species Cover", 
         title = "Change in Herbaceous Species Cover Over Time") + 
    theme_minimal())

## Effect of Plot location ----

cover.glmer <- glmer.nb(Cover ~ SapPlot + (1|Year) + (1 + Year|Taxa4Analysis), data = herb_cover)

## Effect of climate ----

### Temperature ----

####  Adding mean temperature values to create a mixed dataset

herb_climate <- filter(herb_cover, Year != "2018")
colnames(herb_climate)[colnames(herb_climate) == "Year"] <- "year"
herb_climate <- herb_climate %>% mutate(mean_temp = year)

herb_climate$mean_temp <- gsub("2008", "6.30", herb_climate$mean_temp)
herb_climate$mean_temp <- gsub("2009", "6.35", herb_climate$mean_temp)
herb_climate$mean_temp <- gsub("2010", "3.73", herb_climate$mean_temp)
herb_climate$mean_temp <- gsub("2011", "7.32", herb_climate$mean_temp)
herb_climate$mean_temp <- gsub("2012", "8.63", herb_climate$mean_temp)
herb_climate$mean_temp <- gsub("2013", "6.09", herb_climate$mean_temp)
herb_climate$mean_temp <- gsub("2014", "4.92", herb_climate$mean_temp)
herb_climate$mean_temp <- gsub("2015", "7.50", herb_climate$mean_temp)
herb_climate$mean_temp <- gsub("2016", "6.78", herb_climate$mean_temp)
str(herb_climate)
herb_climate$mean_temp <- as.numeric(herb_climate$mean_temp)

#### Cover vs Temperature Model

cover.temp.glmer <- glmer.nb(Cover ~ mean_temp + (1|year) + (1|PlotID) + (1+year|Taxa4Analysis), data = herb_climate)

summary(cover.temp.glmer)
plot(cover.temp.glmer) 
AICc(cover.temp.glmer)
r.squaredGLMM(cover.temp.glmer)
dispersion_glmer(cover.temp.glmer)
stargazer(cover.temp.glmer, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

#### Plotting cover against temperature

herb <- aggregate(Cover ~ year + mean_temp, herb_climate, mean )
cover.temp.glmer <- lm(Cover ~ mean_temp, data = herb)
pred.mm2 <- ggpredict(cover.temp.glmer, terms = c("mean_temp"))  # this gives overall predictions for the model

(ggplot(pred.mm2) + 
    geom_line(aes(x = x, y = predicted), col="black") +          # slope
    geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                fill = "#9BCD9B", alpha = 0.5) +  # error band
    geom_point(data = herb,                      # adding the raw data (scaled values)
               aes(x = mean_temp, y = Cover), col = "#9BCD9B")) + 
  labs(x = "Mean Temperatures", y = "Mean Percentage Cover", 
       title = "Effect of Temperatures on Herbaceous Species Cover") + 
  theme_minimal()

### Precipitation ----

#### Adding mean precipitation values to create a mix dataset

herb_climate <- herb_cover
colnames(herb_climate)[colnames(herb_climate) == "Year"] <- "year"
herb_climate <- herb_climate %>% mutate(rain = year)
herb_climate$rain <- gsub("2008", "107", herb_climate$rain)
herb_climate$rain <- gsub("2009", "105", herb_climate$rain)
herb_climate$rain <- gsub("2010", "97", herb_climate$rain)
herb_climate$rain <- gsub("2011", "121", herb_climate$rain)
herb_climate$rain <- gsub("2012", "108", herb_climate$rain)
herb_climate$rain <- gsub("2013", "156", herb_climate$rain)
herb_climate$rain <- gsub("2014", "82", herb_climate$rain)
herb_climate$rain <- gsub("2015", "98", herb_climate$rain)
herb_climate$rain <- gsub("2016", "75", herb_climate$rain)
herb_climate$rain <- gsub("2017", "99", herb_climate$rain)
herb_climate$rain <- gsub("2018", "90", herb_climate$rain)
str(herb_climate)
herb_climate$rain <- as.numeric(herb_climate$rain)

#### Precipitation Model

cover.rain.glmer <- glmer.nb(Cover ~ rain + (1|year) + (1|PlotID) + (1+year|Taxa4Analysis), data = herb_climate)

summary(cover.rain.glmer)
plot(cover.rain.glmer) 
AICc(cover.temp.glmer)
r.squaredGLMM(cover.temp.glmer)
dispersion_glmer(cover.temp.glmer) 
stargazer(cover.temp.glmer, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

# Individual Species Cover ----

# Select single species

stel_cover <- herb_cover[herb_cover$Taxa4Analysis == "Stellaria_corei",]

## Look at distribution

hist(stel_cover$Cover, col = "#F0E68C", breaks=40, xlim = c(0, 40))

## Fitting Model

str(stel_cover)
stel_cover$SapPlot <- as.factor(stel_cover$SapPlot)

stel.glmer <- glmer(Cover ~ as.numeric(Year) + (1|SapPlot/PlotID), data = stel_cover)
summary(stel.glmer)
plot(stel.glmer)

AICc(stel.glmer)
dispersion_glmer(stel.glmer)
r.squaredGLMM(stel.glmer)
stargazer(stel.glmer, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

# Richness Statistics ----

## Calculate richness ----

herb_rich <- herb_cover[!herb_cover$Cover == 0, ]
herb_rich <- filter(herb_rich, Cover != "NA") 

richness <- herb_rich %>%
  group_by(SapPlot, PlotID, Year) %>%
  summarise(richness = length(unique(Taxa4Analysis)))

richness <- as.data.frame(richness)
richness$Year <- as.factor(richness$Year)
richness$SapPlot <- as.factor(richness$SapPlot)
herb_richness <- aggregate( richness ~ SapPlot + PlotID + Year, richness, mean )
str(herb_richness)

## Fitting model -----

hist(herb_richness$richness, col = "#F0E68C", breaks=10)

#### Null model

null.richness <- lm(richness ~ 1, data = herb_richness)
AICc(null.richness)

### Mixed model 

mixed.year <- lmer(richness ~ as.numeric(Year) + (1|SapPlot/PlotID), data = herb_richness)
summary(mixed.year)
qqnorm(resid(mixed.year))
qqline(resid(mixed.year))  
AICc(mixed.year)
r.squaredGLMM(mixed.year)

car::Anova(mixed.year, type=3)

stargazer(mixed.year, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

## Plot richness against time

### Extract the prediction data frame
pred.mm <- ggpredict(mixed.year, terms = c("Year"))  # this gives overall predictions for the model

### Plot the predictions

(ggplot(pred.mm) + 
    geom_line(aes(x = x, y = predicted), col="black") +          # slope
    geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                fill = "#9BCD9B", alpha = 0.5) +  # error band
    geom_point(data = herb_richness,                      # adding the raw data (scaled values)
               aes(x = as.numeric(Year), y = richness), col = "#9BCD9B") + 
    labs(x = "Year", y = "Species richness", 
         title = "Change in Herbaceous Species Richness Over Time") + 
    theme_minimal()+
    scale_x_continuous(breaks = seq(1,10,1), labels = c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2018")))

### Plot mean richness against time

herby <- aggregate(richness ~ Year, herb_richness, mean)
tot_rich <- lm(richness ~ as.numeric(Year), data = herb_richness)
summary(tot_rich)
pred.mc <- ggpredict(tot_rich, terms = c("Year"))

(ggplot(pred.mc) + 
    geom_line(aes(x = x, y = predicted), col="#494A49") +          # slope
    geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                fill = "#9BCD9B", alpha = 0.5) +  # error band
    geom_point(data = herby,                      # adding the raw data (scaled values)
               aes(x = as.numeric(Year), y = richness), col = "black") + 
    labs(x = "Year", y = "Mean Species Richness") + 
    geom_line(data = herby,                      # adding the raw data (scaled values)
              aes(x = as.numeric(Year), y = richness), col = "#87AB87DF") + 
    theme_minimal()+
    scale_x_continuous(breaks = seq(1,10,1), labels = c("2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2018")))

## Effect of climate ----

### Adding climate data to diversity dataframe

# Temperatures

herb_climate <- herb_richness
herb_climate <- herb_climate %>% mutate(mean_temp = Year)
herb_climate <- filter(herb_climate, Year != "2018")
herb_climate$mean_temp <- gsub("2008", "6.30", herb_climate$mean_temp)
herb_climate$mean_temp <- gsub("2009", "6.35", herb_climate$mean_temp)
herb_climate$mean_temp <- gsub("2010", "3.73", herb_climate$mean_temp)
herb_climate$mean_temp <- gsub("2011", "7.32", herb_climate$mean_temp)
herb_climate$mean_temp <- gsub("2012", "8.63", herb_climate$mean_temp)
herb_climate$mean_temp <- gsub("2013", "6.09", herb_climate$mean_temp)
herb_climate$mean_temp <- gsub("2014", "4.92", herb_climate$mean_temp)
herb_climate$mean_temp <- gsub("2015", "7.50", herb_climate$mean_temp)
herb_climate$mean_temp <- gsub("2016", "6.78", herb_climate$mean_temp)
str(herb_climate)
herb_climate$mean_temp <- as.numeric(herb_climate$mean_temp)
herb_climate$Year <- as.factor(herb_climate$Year)
herb_climate$SapPlot <- as.factor(herb_climate$SapPlot)

# Precipitation

herb_climate <- herb_climate %>% mutate(rain = Year)
herb_climate$rain <- gsub("2008", "107", herb_climate$rain)
herb_climate$rain <- gsub("2009", "105", herb_climate$rain)
herb_climate$rain <- gsub("2010", "97", herb_climate$rain)
herb_climate$rain <- gsub("2011", "121", herb_climate$rain)
herb_climate$rain <- gsub("2012", "108", herb_climate$rain)
herb_climate$rain <- gsub("2013", "156", herb_climate$rain)
herb_climate$rain <- gsub("2014", "82", herb_climate$rain)
herb_climate$rain <- gsub("2015", "98", herb_climate$rain)
herb_climate$rain <- gsub("2016", "75", herb_climate$rain)
herb_climate$rain <- gsub("2017", "99", herb_climate$rain)
herb_climate$rain <- gsub("2018", "90", herb_climate$rain)
str(herb_climate)
herb_climate$rain <- as.numeric(herb_climate$rain)

### Temperature stats ----

mixed.temp <- lmer(richness ~ mean_temp + (1|Year) + (1|SapPlot/PlotID), data = herb_climate)
#mixed.temp <- lmerTest::lmer(richness ~ mean_temp + (1|Year) + (1|SapPlot/PlotID), data = herb_climate)
plot(mixed.temp)
qqnorm(resid(mixed.temp))
qqline(resid(mixed.temp)) 
summary(mixed.temp)
AICc (mixed.temp)
r.squaredGLMM(mixed.temp)
stargazer(mixed.temp, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

#car::Anova(mixed.temp, type=3)

### Precipitation stats ----

mixed.rain <- lmer(richness ~ rain + (1|Year) + (1|SapPlot/PlotID), data = herb_climate)
#mixed.rain <- lmerTest::lmer(richness ~ rain + (1|Year) + (1|SapPlot/PlotID), data = herb_climate)
plot(mixed.rain)
qqnorm(resid(mixed.rain))
qqline(resid(mixed.rain)) 
summary(mixed.rain)
AICc (mixed.rain)
stargazer(mixed.rain, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")
r.squaredGLMM(mixed.rain)

#car::Anova(mixed.rain, type=3)

### Plot location -----

richness_plot <- lmer(richness ~ SapPlot + (1|Year) + (1|SapPlot/PlotID), data = herb_richness)
summary(richness_plot)
qqnorm(resid(richness_plot))
qqline(resid(richness_plot))  
AICc(richness_plot)
r.squaredGLMM(richness_plot)

stargazer(richness_plot, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

#car::Anova(richness_plot, type=3)

# Diversity ----

## Measuring diversity ----

herb_div <- spread(herb_cover, Taxa4Analysis, Cover)
str(herb_div)
herb_matrix <- herb_div %>% dplyr::select(-c("Year", "PlotID", "SapPlot", "BroaderGroup"))
herb_matrix[is.na(herb_matrix)] <- 0
str(herb_matrix)

shannon <- diversity(herb_matrix) # Shannon's Index is the default one so no need to specify which method
shannon

herb_div$diversity <- shannon    
str(herb_div)
herb_div$SapPlot <- as_factor(herb_div$SapPlot)
herb_div[is.na(herb_div)] <- 0

## Statistics ----

hist(herb_div$diversity, col = "#F0E68C", breaks=10)

#### Null model

null.div <- glm(diversity ~ 1, data = herb_div)
AICc(null.div)

### Mixed model 

mixed.div <- glmer(diversity ~ Year + (1|SapPlot/PlotID), data = herb_div)

summary(mixed.div)
plot(mixed.div)
AICc(mixed.div)
dispersion_glmer(mixed.div)
r.squaredGLMM(mixed.div)

car::Anova(mixed.div, type=3)

stargazer(mixed.div, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")
## Plot

pred.div <- ggpredict(mixed.div, terms = c("Year"))  # this gives overall predictions for the model

(ggplot(pred.div) + 
    geom_line(aes(x = x, y = predicted), col= "black") +          # slope
    geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                fill = "#8DB6CD", alpha = 0.5) +  # error band
    geom_point(data = herb_div,                      # adding the raw data (scaled values)
               aes(x = Year, y = diversity), col = "#8DB6CD") + 
    labs(x = "Year", y = "Shannon's Diversity Index", 
         title = "Change in Herbaceous Species Diversity Over Time") + 
    theme_minimal()+
    scale_x_continuous(breaks=seq(2008, 2018, 1))) 

### Plot mean diversity against time

herbo <- aggregate(diversity ~ Year, herb_div, mean)
tot_div <- glm(diversity ~ as.numeric(Year), data = herbo)
summary(tot_div)
pred.mc <- ggpredict(tot_div, terms = c("Year"))

(ggplot(pred.mc) + 
    geom_line(aes(x = x, y = predicted), col="black") +          # slope
    geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                fill = "#8DB6CD", alpha = 0.5) +  # error band
    geom_point(data = herbo,                      # adding the raw data (scaled values)
               aes(x = as.numeric(Year), y = diversity), col = "black") + 
    labs(x = "Year", y = "Mean Species Richness") + 
    geom_line(data = herbo,                      # adding the raw data (scaled values)
              aes(x = as.numeric(Year), y = diversity), col = "#8DB6CD") + 
    theme_minimal()+
    scale_x_continuous(breaks=seq(2008, 2018, 1)))
  
## Effect of climate ----

### Adding climate data to diversity dataframe

# Temperatures

herb_div <- herb_div %>% dplyr::select(-c("Group", "Duration"))
herb_climate <- herb_div
herb_climate <- herb_climate %>% mutate(mean_temp = Year)
herb_climate <- filter(herb_climate, Year != "2018")
herb_climate$mean_temp <- gsub("2008", "6.30", herb_climate$mean_temp)
herb_climate$mean_temp <- gsub("2009", "6.35", herb_climate$mean_temp)
herb_climate$mean_temp <- gsub("2010", "3.73", herb_climate$mean_temp)
herb_climate$mean_temp <- gsub("2011", "7.32", herb_climate$mean_temp)
herb_climate$mean_temp <- gsub("2012", "8.63", herb_climate$mean_temp)
herb_climate$mean_temp <- gsub("2013", "6.09", herb_climate$mean_temp)
herb_climate$mean_temp <- gsub("2014", "4.92", herb_climate$mean_temp)
herb_climate$mean_temp <- gsub("2015", "7.50", herb_climate$mean_temp)
herb_climate$mean_temp <- gsub("2016", "6.78", herb_climate$mean_temp)
str(herb_climate)
herb_climate$mean_temp <- as.numeric(herb_climate$mean_temp)
herb_climate$Year <- as.factor(herb_climate$Year)
herb_climate$SapPlot <- as.factor(herb_climate$SapPlot)

# Precipitation

herb_climate <- herb_climate %>% mutate(rain = Year)
herb_climate$rain <- gsub("2008", "107", herb_climate$rain)
herb_climate$rain <- gsub("2009", "105", herb_climate$rain)
herb_climate$rain <- gsub("2010", "97", herb_climate$rain)
herb_climate$rain <- gsub("2011", "121", herb_climate$rain)
herb_climate$rain <- gsub("2012", "108", herb_climate$rain)
herb_climate$rain <- gsub("2013", "156", herb_climate$rain)
herb_climate$rain <- gsub("2014", "82", herb_climate$rain)
herb_climate$rain <- gsub("2015", "98", herb_climate$rain)
herb_climate$rain <- gsub("2016", "75", herb_climate$rain)
herb_climate$rain <- gsub("2017", "99", herb_climate$rain)
herb_climate$rain <- gsub("2018", "90", herb_climate$rain)
str(herb_climate)
herb_climate$rain <- as.numeric(herb_climate$rain)

### Temperature stats ----

mixed.temp <- lmer(diversity ~ mean_temp + (1|Year) + (1|SapPlot/PlotID), data = herb_climate)
plot(mixed.temp)
summary(mixed.temp)
AICc (mixed.temp)
r.squaredGLMM(mixed.temp)
stargazer(mixed.temp, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

car::Anova(mixed.temp, type=3)

### Precipitation stats ----

mixed.rain <- lmer(diversity ~ rain + (1|Year) + (1|SapPlot/PlotID), data = herb_climate)
plot(mixed.rain)
summary(mixed.rain)
AICc (mixed.rain)
stargazer(mixed.rain, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")
r.squaredGLMM(mixed.rain)

#car::Anova(mixed.rain, type=3)

### Plot location -----

div_plot <- lmer(diversity ~ SaPlot + (1|Year) + (1|SapPlot/PlotID), data = herb_div)
summary(div_plot)
qqnorm(resid(div_plot))
qqline(resid(div_plot))  
AICc(div_plot)
r.squaredGLMM(div_plot)

stargazer(div_plot, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

car::Anova(div_plot, type=3)
