### Script to generate regularized estimates of population trends using the
### grouping structures from the Rosenberg et al 2019 paper
### Breeding and wintering biome

library(tidyverse)
library(jagsUI)
library(tidybayes)



data <- read.csv("data/NAbirdbase_select.csv")
rosen = read.csv("data/Rosenberg et al species list.csv")


### funciton to transform %-change trends to log-scale geometric mean rates of change
log_trans <- function(x){
  log((x/100)+1)
}
## transforming the 95% CIs into an estimate of variance fo the trend estimate
log_var_func <- function(uci,lci){
  ((log_trans(uci)-log_trans(lci))/(1.96*2))^2
}

dat <- data[,c("CommonName","trendyr","lcl_90CI","ucl_90CI")]
 
## I've dropped the Breeding Biome labels from the original data because they don't match the groupings in Rosenberg et al.

dat <- dat %>% drop_na()

rosen_sp <- rosen$species
dat_sp <- dat$CommonName

rosen_sp[-which(rosen_sp %in% dat_sp)]
#Not sure why these 53 species are missing trend estimates

dat_sp[-which(dat_sp %in% rosen_sp)]
#Not sure why these 53 species are missing trend estimates


rosen_groups <- rosen %>% 
  select(species,Breeding.Biome,Winter.Biome) %>% 
  rename(BreedingBiome = Breeding.Biome,
         NonbreedingBiome = Winter.Biome)

dat <- dat %>% 
  left_join(.,rosen_groups,by = c("CommonName" = "species")) %>% 
  drop_na() %>% 
  mutate(log_trend = log_trans(trendyr),
         var_trend = log_var_func(ucl_90CI,lcl_90CI),
         BB_f = as.integer(factor(BreedingBiome)),
         NBB_f = as.integer(factor(NonbreedingBiome)))

# table(dat$BreedingBiome)
# table(dat$NonbreedingBiome)



   jags_data <- list(nspecies = nrow(dat),
                     ngroups1 = max(dat$BB_f),
                     ngroups2 = max(dat$NBB_f),
                     group1 = dat$BB_f,
                     group2 = dat$NBB_f,
                     betahat = dat$log_trend,
                     varhat = dat$var_trend)
   
   
  
   # model -------------------------------------------------------------------
   
   modl <- "

### model for shrinking population trend estimates based on 
## breeding biome (group1) and wintering biome (group2)

## required data
## varhat = estaimtes of the variance of the trend estimates for each species
## betahat = estimates of the trends for each species
## nspecies = number of species
## ngroups1 = number of breeding biome groups
## ngroups2 = number of wintering biome groups
## group1 = vector of breeding biome indices for each species
## group2 = vector of wintering biome indices for each species
## 

model{

for(s in 1:nspecies) {
  tau.betahat[s] <- 1/varhat[s] #transform variance to precision
	betahat[s] ~ dnorm(beta[s],tau.betahat[s] ) #betahat = data = trend estimates, tau.betahat = data = precision of trend estimate

 grp1[s,group1[s]] ~ dnorm(mu_grp1[group1[s]],tau_g1[group1[s]]) #random effect
 grp2[s,group2[s]] ~ dnorm(mu_grp2[group2[s]],tau_g2[group2[s]]) #random effect

  beta[s] <- grp1[s,group1[s]]+grp2[s,group2[s]] #shrunken trends

} # end of s loop


for(g in 1:ngroups1){
mu_grp1[g] ~ dnorm(0,0.1) #fixed effect mean for group
sd_g1[g] ~ dnorm(0,0.1) T(0,) #sd of species trends within group
tau_g1[g] <- 1/pow(sd_g1[g],2) #precision of above
}

for(g in 1:ngroups2){
mu_grp2[g] ~ dnorm(0,0.1) #fixed effect mean for group
sd_g2[g] ~ dnorm(0,0.1) T(0,) #sd of species trends within group
tau_g2[g] <- 1/pow(sd_g2[g],2) #precision of above
}







} # end of model

   "  

trend_comp = "models/trend_shrinkage_model.txt"
cat(modl,file = trend_comp)   






params <- c("beta",
            "sd_g2",
            "mu_grp2",
            "mu_grp1",
            "sd_g1",
            "grp1",
            "grp2")


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
hist(betas_summary$Rhat)
hist(betas_summary$n.eff)



beta_samples <- gather_draws(out$samples,beta[sp])

sps <- data.frame(CommonName = dat$CommonName,
                  sp = 1:nrow(dat))
new_trends <- beta_samples %>% group_by(sp) %>% 
  summarise(new_trend = mean((exp(.value)-1)*100),
            new_trend_lci = quantile((exp(.value)-1)*100,0.025),
            new_trend_uci = quantile((exp(.value)-1)*100,0.975)) %>% 
  left_join(.,sps,by = "sp") %>% 
  left_join(.,dat,by = "CommonName")

comp_plot <- ggplot(data = new_trends,aes(x = trendyr,y = new_trend))+
  geom_point(aes(colour = BreedingBiome,size = var_trend))+
  coord_cartesian(xlim = c(-7,10),ylim = c(-7,10))+
  geom_abline(slope = 1,intercept = 0,alpha = 0.2)+
  theme_classic()
  
print(comp_plot)
