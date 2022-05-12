### Script to generate regularized estimates of BCR population trends 

### includes downloaded 2019 BBS trends from usgs trends website: https://www.mbr-pwrc.usgs.gov/
## two downloaded files
## BBS_1993-2019_expanded_trend_best.csv = best estimate of trends for expanded survey area 1993 - 2019
## BBS_1966-2019_core_best_trend.csv = best estimate of trends for core survey area 1993 - 2019 (excludes some northern areas where data collection was insufficient prior to 1993)

### NOTE: this script uses the 1966-2019 trends, but it would work for the expanded survey area as well



library(tidyverse)
library(jagsUI)
library(tidybayes)
library(ggrepel)
source("Functions/posterior_summary_functions.R")


tbcr <- read.csv("data/BBS_1966-2019_core_best_trend.csv")
#remove the "all forms" qualification from the lumped BBS species
tbcr <- tbcr %>% 
  mutate(Species.Name = gsub(x = trimws(Species.Name),pattern = " (all forms)",replacement = "",fixed = TRUE),
         Credibility.Code = factor(trimws(Credibility.Code),ordered = TRUE,
                                   levels = c("","Y","R"),
                                   labels = c("G","Y","R")))

### funciton to transform %-change trends to log-scale geometric mean rates of change
log_trans <- function(x){
  log((x/100)+1)
}
## transforming the 95% CIs into an estimate of variance fo the trend estimate
log_var_func <- function(uci,lci){
  ((log_trans(uci)-log_trans(lci))/(1.96*2))^2
}

sps <- unique(tbcr$Species.Name)

lumped <- sps[grepl(sps,pattern = " & ")]

sps_spl <- data.frame(Species.Name = rep(lumped,each = 2),
                      species = c("Western Grebe","Clark's Grebe",
                                  "Alder Flycatcher","Willow Flycatcher",
                                  "Cordilleran Flycatcher","Pacific-slope Flycatcher",
                                  "California Scrub-Jay","Woodhouse's Scrub-Jay",
                                  "Sagebrush Sparrow","Bell's Sparrow"))

# Lumped species are complicated because BBS doesn't provide single species estim --------
### this could be fixed if one were to identify which of the species occurs in each BCR
### or alternatively, if one were to decide to include only one in the groupings

regs <- tbcr %>% 
  select(Region,Region.Name) %>% 
  filter(grepl(x = Region,pattern = "^[AS][[:digit:]]")) %>% 
  distinct()

# model -------------------------------------------------------------------

modl <- "

### model for shrinking population trend estimates towards a global mean

## required data
## varhat = log-scale estaimtes of the variance of the trend estimates for each species
## betahat = log-scale estimates of the trends for each species
## nspecies = number of species


model{

for(s in 1:nspecies) {
  tau.betahat[s] <- 1/varhat[s] #transform variance to precision
	betahat[s] ~ dnorm(beta[s],tau.betahat[s] ) #betahat = data = trend estimates, tau.betahat = data = precision of trend estimate

 beta[s] ~ dnorm(mu,tau) #

} # end of s loop


mu ~ dnorm(0,0.5) #fixed effect mean for group
tau ~ dgamma(2,0.01)
sd <- 1/pow(tau,0.5) #precision of above










} # end of model

   "  

trend_comp = "models/trend_shrinkage_model_simple.txt"
cat(modl,file = trend_comp)   




pdf("BCR_trend_shrinkage.pdf")

out_trends <- NULL

for(i in 1:nrow(regs)){
  bcr = regs[i,"Region.Name"]
  b = regs[i,"Region"]
  
dat <- tbcr %>% 
  filter(Region.Name == bcr) %>% 
  mutate(log_trend = log_trans(Trend),
         var_trend = log_var_func(X97.5.CI,X2.5.CI))

# table(dat$BreedingBiome)
# table(dat$NonbreedingBiome)



   jags_data <- list(nspecies = nrow(dat),
                     betahat = dat$log_trend,
                     varhat = dat$var_trend)
   
   
  





params <- c("beta",
            "sd",
            "mu")


burnInSteps = 5000            # Number of steps to "burn-in" the samplers. this is sufficient for testing, but you'll want to increase this
nChains = 3                   # Number of chains to run.
numSavedSteps=1000         # Total number of steps in each chain to save. this is sufficient for testing, but you'll want to increase this
thinSteps=50                   # Number of steps to "thin" (1=keep every step).
nIter = ceiling( ( (numSavedSteps * thinSteps )+burnInSteps)) # Steps per chain.



out = jagsUI(data = jags_data,
             parameters.to.save = params,
             n.chains = 3,
             n.burnin = burnInSteps,
             n.thin = thinSteps,
             n.iter = nIter,
             parallel = T,
             model.file = trend_comp)


betas_summary <- as.data.frame(out$summary[paste0("beta[",1:jags_data$nspecies,"]"),])
# hist(betas_summary$Rhat)
# hist(betas_summary$n.eff)





beta_samples <- posterior_samples(out$samples,
                                  parm = "beta",
                                  dims = "sp",
                                  is_mcmc = TRUE)

sps <- data.frame(Species.Name = dat$Species.Name,
                  sp = 1:nrow(dat))
shrunk_trends <- beta_samples %>% group_by(sp) %>% 
  summarise(shrunk_trend = mean((exp(.value)-1)*100),
            shrunk_trend_lci = quantile((exp(.value)-1)*100,0.025),
            shrunk_trend_uci = quantile((exp(.value)-1)*100,0.975)) %>% 
  left_join(.,sps,by = "sp") %>% 
  left_join(.,dat,by = "Species.Name") %>% 
  mutate(t_shrink = abs(shrunk_trend - Trend),
         Half_CI_Width_original = (X97.5.CI-X2.5.CI)/2)

labl <- shrunk_trends %>% 
  filter(t_shrink > 2.5)

comp_plot <- ggplot(data = shrunk_trends,aes(x = Trend,y = shrunk_trend))+
  geom_point(aes(colour = Credibility.Code))+
  geom_errorbarh(data = labl,
                aes(y = shrunk_trend,xmin = X2.5.CI,xmax = X97.5.CI),
                alpha = 0.1,height = 0)+
  coord_cartesian(xlim = c(-15,20))+
  geom_abline(slope = 1,intercept = 0,alpha = 0.2)+
  scale_size_continuous(range = c(1,3),
                        trans = scales::trans_new("sd_prec",
                                                  transform = function(x){1/(x^2)},
                                                  inverse = function(x){1/(x^0.5)}))+
  geom_text_repel(data = labl,
                  aes(x = Trend,y = shrunk_trend,
                      label = Species.Name),
                  inherit.aes = FALSE,
                  min.segment.length = 0,
                  size = 4,
                  direction = "both",
                  nudge_x = 3)+
scale_colour_viridis_d(end = 0.9)+
  labs(title = bcr,
       caption = "Labeled species trends shrink > 2%/year and show error bars.\n Distance from 1:1 diagonal line reflects shrinkage")+
  xlab("Original Trend Estimate")+
  ylab("Posterior (shrunk) Trend Estimate")+
  theme_classic()
  
print(comp_plot)



out_trends <- bind_rows(out_trends,shrunk_trends)

print(round(i/nrow(regs),2))
}

dev.off()

tren <- out_trends %>% 
  select(Species.Name, AOU, Region, Region.Name,
         shrunk_trend,shrunk_trend_lci,shrunk_trend_uci,
         Trend,X2.5.CI,X97.5.CI) %>% 
  rename(Original_Trend = Trend,
         Original_Trend_lci = X2.5.CI,
         Original_Trend_uci = X97.5.CI) %>% 
  arrange(Region.Name,AOU)

write.csv(tren,
          file = "Shrunken_trends_BBS_by_BCR_global_mean.csv",
          row.names = FALSE)
