## load packages
library(ape)
library(phytools)
library(brms)
library(rethinking)
library(rstan)

# rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


## set working directory
setwd("//files.iem.uzh.ch/Data/Institute/Human_Ecology/ajaegg/Private/My Documents/GitHub/Coalitions-by-sex-mammals")

## read files
d<- read.csv("sexbiascoalitions.csv", header=TRUE)
trees<- read.nexus("TreeBlockSexBiasCooperation.nex")

# convert trees to covariance matrix (see https://cran.r-project.org/web/packages/brms/vignettes/brms_phylogenetics.html)
A<- list() 
for(i in 1:100){
  A[[i]]<- vcv.phylo(trees[[i]])
}




### Model 1: Intercept only ###

## set prior
prior.m1<- c(
  prior(normal(0,1), "Intercept", dpar="muFemales"), # favors moderate probabilities, penalizes extremes (near 0 or 1)
  prior(normal(0,1), "Intercept", dpar="muMales"),
  prior(exponential(1), "sd", dpar="muFemales"), 
  prior(exponential(1), "sd", dpar="muMales")
)

m1_loop1<- brm(coalitions ~ 1 + (1|Genus_species), data = d, family = categorical(), cov_ranef = list(Genus_species = A[[1]]), 
         prior = prior.m1, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.95, max_treedepth = 15))
# no divergent transitions
plot(m1_loop1) # good convergence, chains overlapping and stationary
summary(m1_loop1) # all Rhat 1, ESS >1000


# extract samples
post_m1_looped<- posterior_samples(m1_loop1)

# loop and add samples
for(i in 2:100){
  m1_loopi<- brm(coalitions ~ 1 + (1|Genus_species), data = d, family = categorical(), cov_ranef = list(Genus_species = A[[i]]), 
                 prior = prior.m1, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.95, max_treedepth = 15))

post_m1_loopi<- posterior_samples(m1_loopi)
post_m1_looped<- rbind(post_m1_looped, post_m1_loopi)
save(post_m1_looped, file="post_m1_looped.robj")
}

nrow(post_m1_looped)/8000 # check all 100 models added 



### Model 2: Sex-segregation ###

get_prior(coalitions ~ SexComposition + (1|Genus_species), data = d, family = categorical(), cov_ranef = list(Genus_species = A.cons))
## adjust prior
prior.m2<- c(
  prior(normal(0,1), "Intercept", dpar="muFemales"), 
  prior(normal(0,1), "Intercept", dpar="muMales"),
  prior(normal(0,0.5), "b", coef = "SexCompositionsegregated", dpar="muFemales"),
  prior(normal(0,0.5), "b", coef = "SexCompositionsegregated", dpar="muMales"),
  prior(exponential(1), "sd", dpar="muFemales"), 
  prior(exponential(1), "sd", dpar="muMales")
)

m2_loop1<- brm(coalitions ~ SexComposition + (1|Genus_species), data = d, family = categorical(), cov_ranef = list(Genus_species = A[[1]]), 
         prior = prior.m2, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.95, max_treedepth = 15))
# no divergent transitions
plot(m2_loop1) # good convergence, chains overlapping and stationary
summary(m2_loop1) # all Rhat 1, ESS >1000

# extract samples
post_m2_looped<- posterior_samples(m2_loop1)

# loop and add samples
for(i in 2:100){
  m2_loopi<- brm(coalitions ~ SexComposition + (1|Genus_species), data = d, family = categorical(), cov_ranef = list(Genus_species = A[[i]]), 
                 prior = prior.m2, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.95, max_treedepth = 15))
  
  post_m2_loopi<- posterior_samples(m2_loopi)
  post_m2_looped<- rbind(post_m2_looped, post_m2_loopi)
  save(post_m2_looped, file="post_m2_looped.robj")
}

nrow(post_m2_looped)/8000 # check all 100 models added 




### Model 3: Primates vs others ###

# create dummy
d$Primate<- 0
d$Primate[d$Order=="Primata"]<- 1

# adjust prior
prior.m3<- c(
  prior(normal(0,1), "Intercept", dpar="muFemales"), 
  prior(normal(0,1), "Intercept", dpar="muMales"),
  prior(normal(0,0.5), "b", coef = "Primate", dpar="muFemales"),
  prior(normal(0,0.5), "b", coef = "Primate", dpar="muMales"),
  prior(exponential(1), "sd", dpar="muFemales"), 
  prior(exponential(1), "sd", dpar="muMales")
)

m3_loop1<- brm(coalitions ~ Primate + (1|Genus_species), data = d, family = categorical(), cov_ranef = list(Genus_species = A[[1]]), 
         prior = prior.m3, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.95, max_treedepth = 15))
# no divergent transitions
plot(m3_loop1) # good convergence, chains overlapping and stationary
summary(m3_loop1) # all Rhat 1, ESS >1000

# extract samples
post_m3_looped<- posterior_samples(m3_loop1)

# loop and add samples
for(i in 2:100){
  m3_loopi<-  brm(coalitions ~ Primate + (1|Genus_species), data = d, family = categorical(), cov_ranef = list(Genus_species = A[[i]]), 
                  prior = prior.m3, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.95, max_treedepth = 15))
  post_m3_loopi<- posterior_samples(m3_loopi)
  post_m3_looped<- rbind(post_m3_looped, post_m3_loopi)
  save(post_m3_looped, file="post_m3_looped.robj")
}

nrow(post_m3_looped)/8000 # check all 100 models added 



### Model 4: Main predictors ###

# transform predictors
d$f_philopatry<- 1
d$f_philopatry[d$Philopatry=="N" | d$Philopatry=="M"]<- 0

# re-code DV
d$f_coalitions<- 1
d$f_coalitions[d$coalitions=="Males"]<- 0

# how much variance is there actually?
table(d$f_coalitions, d$Food_resource_defendable)
table(d$f_coalitions, d$f_philopatry)

# adjust prior
get_prior(f_coalitions ~ Food_resource_defendable + f_philopatry + (1|Genus_species), data = d, family = bernoulli(), cov_ranef = list(Genus_species = A[[1]]))
prior.m4<- c(
  prior(normal(0,1), "Intercept"), 
  prior(normal(0,0.5), "b"), 
  prior(exponential(1), "sd"))

## run model
m4_loop1<- brm(f_coalitions ~ Food_resource_defendable + f_philopatry + (1|Genus_species), data = d, family = bernoulli(), cov_ranef = list(Genus_species = A[[1]]), 
         prior = prior.m4, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.95, max_treedepth = 15))
# no divergent transitions
plot(m4_loop1) # good convergence, chains overlapping and stationary
summary(m4_loop1) # all Rhat 1, ESS >1000

# extract samples
post_m4_looped<- posterior_samples(m4_loop1)

# loop and add samples
for(i in 2:100){
  m4_loopi<-  brm(f_coalitions ~ Food_resource_defendable + f_philopatry + (1|Genus_species), data = d, family = bernoulli(), cov_ranef = list(Genus_species = A[[i]]), 
                  prior = prior.m4, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.95, max_treedepth = 15))
  post_m4_loopi<- posterior_samples(m4_loopi)
  post_m4_looped<- rbind(post_m4_looped, post_m4_loopi)
  save(post_m4_looped, file="post_m4_looped.robj")
}

nrow(post_m4_looped)/8000 # check all 100 models added 



## robustness checks: each predictor on its own
m4a_loop1<- brm(f_coalitions ~ Food_resource_defendable + (1|Genus_species), data = d, family = bernoulli(), cov_ranef = list(Genus_species = A[[1]]), 
          prior = prior.m4, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.95, max_treedepth = 15))
# no divergent transitions
plot(m4a_loop1)
summary(m4a_loop1)

# extract samples
post_m4a_looped<- posterior_samples(m4a_loop1)

# loop and add samples
for(i in 2:100){
  m4a_loopi<-  brm(f_coalitions ~ Food_resource_defendable + (1|Genus_species), data = d, family = bernoulli(), cov_ranef = list(Genus_species = A[[i]]), 
                  prior = prior.m4, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.95, max_treedepth = 15))
  post_m4a_loopi<- posterior_samples(m4a_loopi)
  post_m4a_looped<- rbind(post_m4a_looped, post_m4a_loopi)
  save(post_m4a_looped, file="post_m4a_looped.robj")
}

nrow(post_m4a_looped)/8000 # check all 100 models added 



m4b_loop1<- brm(f_coalitions ~ f_philopatry + (1|Genus_species), data = d, family = bernoulli(), cov_ranef = list(Genus_species = A[[1]]), 
          prior = prior.m4, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.95, max_treedepth = 15))
# no divergent transitions
plot(m4b_loop1)
summary(m4b_loop1)

# extract samples
post_m4b_looped<- posterior_samples(m4b_loop1)

# loop and add samples
for(i in 2:100){
  m4b_loopi<-  brm(f_coalitions ~ f_philopatry + (1|Genus_species), data = d, family = bernoulli(), cov_ranef = list(Genus_species = A[[i]]), 
                   prior = prior.m4, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.95, max_treedepth = 15))
  post_m4b_loopi<- posterior_samples(m4b_loopi)
  post_m4b_looped<- rbind(post_m4b_looped, post_m4b_loopi)
  save(post_m4b_looped, file="post_m4b_looped.robj")
}

nrow(post_m4b_looped)/8000 # check all 100 models added 



### Model 5: Male coalitions ~ Dimorphism, Philopatry ###

# transform predictors
d$SexDim.c<- d$Sex_Dim-1 # mean = monomorphic
d$m_philopatry<- 1
d$m_philopatry[d$Philopatry=="N" | d$Philopatry=="F"]<- 0

# re-code DV
d$m_coalitions<- 1
d$m_coalitions[d$coalitions=="Females"]<- 0

# examine variation
table(d$m_coalitions, d$m_philopatry)

d$m_philopatry<- as.factor(d$m_philopatry)

# same prior as for m4
m5_loop1<- brm(m_coalitions ~ SexDim.c + m_philopatry + (1|Genus_species), data = d, family = bernoulli(), cov_ranef = list(Genus_species = A[[1]]), 
         prior = prior.m4, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.99, max_treedepth = 15))
# no divergent transitions
plot(m5_loop1) # good convergence, chains overlapping and stationary
summary(m5_loop1) # all Rhat 1, ESS >1000

# extract samples
post_m5_looped<- posterior_samples(m5_loop1)

# loop and add samples
for(i in 2:100){
  m5_loopi<-  brm(m_coalitions ~ SexDim.c + m_philopatry + (1|Genus_species), data = d, family = bernoulli(), cov_ranef = list(Genus_species = A[[i]]), 
                  prior = prior.m4, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.95, max_treedepth = 15))
  post_m5_loopi<- posterior_samples(m5_loopi)
  post_m5_looped<- rbind(post_m5_looped, post_m5_loopi)
  save(post_m5_looped, file="post_m5_looped.robj")
}

nrow(post_m5_looped)/8000 # check all 100 models added 


## robustness checks: each predictor on its own
m5a_loop1<- brm(m_coalitions ~ SexDim.c + (1|Genus_species), data = d, family = bernoulli(), cov_ranef = list(Genus_species = A[[1]]), 
          prior = prior.m4, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.99, max_treedepth = 15))
# no divergent transitions
plot(m5a_loop1)
summary(m5a_loop1)

# extract samples
post_m5a_looped<- posterior_samples(m5a_loop1)

# loop and add samples
for(i in 2:100){
  m5a_loopi<-  brm(m_coalitions ~ SexDim.c + (1|Genus_species), data = d, family = bernoulli(), cov_ranef = list(Genus_species = A[[i]]), 
                  prior = prior.m4, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.95, max_treedepth = 15))
  post_m5a_loopi<- posterior_samples(m5a_loopi)
  post_m5a_looped<- rbind(post_m5a_looped, post_m5a_loopi)
  save(post_m5a_looped, file="post_m5a_looped.robj")
}

nrow(post_m5a_looped)/8000 # check all 100 models added 


m5b_loop1<- brm(m_coalitions ~ m_philopatry + (1|Genus_species), data = d, family = bernoulli(), cov_ranef = list(Genus_species = A[[1]]), 
          prior = prior.m4, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.99, max_treedepth = 15))
# no divergent transitions
plot(m5b_loop1)
summary(m5b_loop1)

# extract samples
post_m5b_looped<- posterior_samples(m5b_loop1)

# loop and add samples
for(i in 2:100){
  m5b_loopi<-  brm(m_coalitions ~ m_philopatry + (1|Genus_species), data = d, family = bernoulli(), cov_ranef = list(Genus_species = A[[i]]), 
                  prior = prior.m4, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.95, max_treedepth = 15))
  post_m5b_loopi<- posterior_samples(m5b_loopi)
  post_m5b_looped<- rbind(post_m5b_looped, post_m5b_loopi)
  save(post_m5b_looped, file="post_m5b_looped.robj")
}

nrow(post_m5b_looped)/8000 # check all 100 models added 




##### Models 4 & 5 as ordinal #####

## m4
d$freq.females<- d$freq.females+1 # cumuative distro requires positive integers
get_prior(freq.females ~ Food_resource_defendable + f_philopatry + (1|Genus_species), data = d, family = cumulative(), cov_ranef = list(Genus_species = A[[1]]))
prior.m4<- c(
  prior(normal(0,1), "Intercept"), 
  prior(normal(0,0.5), "b"), 
  prior(exponential(1), "sd"))

m4.ord_loop1<- brm(freq.females ~ Food_resource_defendable + f_philopatry + (1|Genus_species), data = d, family = cumulative("logit"), cov_ranef = list(Genus_species = A[[1]]), 
         prior = prior.m4, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.95, max_treedepth = 15))
plot(m4.ord_loop1)
summary(m4.ord_loop1)
# all good

# extract samples
post_m4.ord_looped<- posterior_samples(m4.ord_loop1)

# loop and add samples
for(i in 2:100){
  m4.ord_loopi<-  brm(freq.females ~ Food_resource_defendable + f_philopatry + (1|Genus_species), data = d, family = cumulative("logit"), cov_ranef = list(Genus_species = A[[i]]), 
                  prior = prior.m4, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.95, max_treedepth = 15))
  post_m4.ord_loopi<- posterior_samples(m4.ord_loopi)
  post_m4.ord_looped<- rbind(post_m4.ord_looped, post_m4.ord_loopi)
  save(post_m4.ord_looped, file="post_m4.ord_looped.robj")
}

nrow(post_m4.ord_looped)/8000 # check all 100 models added 


## m5
d$freq.males<- d$freq.males+1 # cumuative distro requires positive integers
get_prior(freq.males ~ SexDim.c + m_philopatry + (1|Genus_species), data = d, family = cumulative(), cov_ranef = list(Genus_species = A[[1]]))
prior.m5<- c(
  prior(normal(0,1), "Intercept"), 
  prior(normal(0,0.5), "b"), 
  prior(exponential(1), "sd"))

m5.ord_loop1<- brm(freq.males ~ SexDim.c + m_philopatry + (1|Genus_species), data = d, family = cumulative("logit"), cov_ranef = list(Genus_species = A[[1]]), 
             prior = prior.m5, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.95, max_treedepth = 15))
plot(m5.ord_loop1)
summary(m5.ord_loop1)
# all good

# extract samples
post_m5.ord_looped<- posterior_samples(m5.ord_loop1)

# loop and add samples
for(i in 2:100){
  m5.ord_loopi<-  brm(freq.males ~ SexDim.c + m_philopatry + (1|Genus_species), data = d, family = cumulative("logit"), cov_ranef = list(Genus_species = A[[i]]), 
                      prior = prior.m4, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.95, max_treedepth = 15))
  post_m5.ord_loopi<- posterior_samples(m5.ord_loopi)
  post_m5.ord_looped<- rbind(post_m5.ord_looped, post_m5.ord_loopi)
  save(post_m5.ord_looped, file="post_m5.ord_looped.robj")
}

nrow(post_m5.ord_looped)/8000 # check all 100 models added 

