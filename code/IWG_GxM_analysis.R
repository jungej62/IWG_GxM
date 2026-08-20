#IWG GxM analysis
#Objective: Identify IWG genotypes that yield well under variable N conditions
#Compare NUE among lines
#Look to see if N rates help specific lines yield well through time
library(Rmisc)
library(magrittr)
library(tidyverse)
library(nlme)
library(emmeans)
library(multcomp)
library(googlesheets4)
library(lme4)
library(lmerTest)

dat<-read_sheet("https://docs.google.com/spreadsheets/d/1MaD6Ob8TwMPMW4QM13nPGQBHT5CSNYcD0qj32dqdsSY/edit?gid=1198514056#gid=1198514056", sheet="Analysis", na="NA")
str(dat)
dat %<>% mutate(across(where(is.character), as_factor))
dat$fyear<-as.factor(dat$year)
dat$fnrate<-as.factor(dat$nrate)
dat$frep<-as.factor(dat$rep)
pps = position_dodge(width = .75)

dat %>%
  group_by(fyear, var, fnrate) %>% 
  summarize(sdyld=sd(grain_yld, na.rm=T),
            nyld=n(),
            myld=mean(grain_yld, na.rm=T)) %>% 
  mutate(seyld=sdyld/sqrt(nyld)) %>% 
ggplot(aes(x=var, y=myld, color=var))+
  facet_grid(fnrate~fyear)+
  geom_point(size=1.5, position=pps)+
  geom_errorbar(aes(ymin=myld-seyld, ymax=myld+seyld), width=.2, linewidth=.8, color="black", position=pps)+
  theme_minimal()

#Modeling N rate response as BLUPs. 
#Suggested model won't converge, so trying to scale N rates

dat$nscaled<-as.numeric(scale(dat$nrate))
blup_model<-lmer(
  grain_yld ~ nscaled * fyear + 
    (1 | rep) + 
    (1 | rep:fnrate) + 
    #(1 | rep:fnrate:var) +
    #(1 + nrate | var:fyear) +
    (1 + nscaled | var),
  data = dat,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

summary(blup_model)

#back-scale n rates
sd_N <- sd(dat$nrate, na.rm = TRUE)
mean_N <- mean(dat$nrate, na.rm = TRUE)
fixed_intercept <- fixef(blup_model)["(Intercept)"]
fixed_slope_scaled <- fixef(blup_model)["nscaled"]
genotype_ranefs <- ranef(blup_model)$var

genotype_blups <- ranef(blup_model)$var %>%
  tibble::rownames_to_column(var = "var") %>%
  rename(
    Intercept_Dev = `(Intercept)`, 
    Slope_Dev = nscaled
  ) %>%
  mutate(
    # Total yield at average N rate
    Total_Intercept = fixed_intercept + Intercept_Dev,
    
    # Scaled total slope (per 1 SD of N)
    Total_Slope_Scaled = fixed_slope_scaled + Slope_Dev,
    
    # UNSCALED slope: Yield increase per 1 physical unit of N
    Total_Slope_RealUnits = Total_Slope_Scaled / sd_N
  ) %>%
  arrange(desc(Total_Slope_RealUnits))

head(genotype_blups)

#Trying to plot these in 4 quadrats
avg_intercept <- mean(genotype_blups$Total_Intercept)
avg_slope <- mean(genotype_blups$Total_Slope_RealUnits)

# 2D Selection Matrix
ggplot(genotype_blups, aes(x = Total_Slope_RealUnits, y = Total_Intercept)) +
  geom_point(color = "navy", size = 2, alpha = 0.8) +
  geom_vline(xintercept = avg_slope, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = avg_intercept, linetype = "dashed", color = "gray50") +
  ggrepel::geom_text_repel(
    data = filter(genotype_blups, Total_Slope_RealUnits > avg_slope & Total_Intercept > avg_intercept),
    aes(label = var),
    size = 3,
    max.overlaps = 15
  ) +
  labs(
    title = "Genotype Selection Matrix: Nitrogen Responsiveness vs. Yield",
    subtitle = "Top-right quadrant: High yield potential + high N-responsiveness",
    x = "Nitrogen Slope (Yield Increase per Unit N)",
    y = "Predicted Yield at Mean Nitrogen Rate"
  ) +
  theme_bw()
