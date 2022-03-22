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


# create consensus tree
cons.tree<- consensus.edges(trees)
A.cons<- vcv.phylo(cons.tree)

d$coalitions<- relevel(as.factor(d$coalitions), ref="Both")



### Model 1: Intercept only ###

## set prior
prior.m1<- c(
  prior(normal(0,1), "Intercept", dpar="muFemales"), # favors moderate probabilities, penalizes extremes (near 0 or 1)
  prior(normal(0,1), "Intercept", dpar="muMales"),
  prior(exponential(1), "sd", dpar="muFemales"), 
  prior(exponential(1), "sd", dpar="muMales")
)

m1<- brm(coalitions ~ 1 + (1|Genus_species), data = d, family = categorical(), cov_ranef = list(Genus_species = A.cons), 
         prior = prior.m1, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.999, max_treedepth = 15))
# no divergent transitions
plot(m1) # good convergence, chains overlapping and stationary
summary(m1) # all Rhat 1, ESS >1000

# extract samples and save for post-processing
post_m1<- posterior_samples(m1)
save(post_m1, file="post_m1.robj")


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

m2<- brm(coalitions ~ SexComposition + (1|Genus_species), data = d, family = categorical(), cov_ranef = list(Genus_species = A.cons), 
         prior = prior.m2, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.999, max_treedepth = 15))
# no divergent transitions
plot(m2) # good convergence, chains overlapping and stationary
summary(m2) # all Rhat 1, ESS >1000

# extract samples and save for post-processing
post_m2<- posterior_samples(m2)
save(post_m2, file="post_m2.robj")


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

m3<- brm(coalitions ~ Primate + (1|Genus_species), data = d, family = categorical(), cov_ranef = list(Genus_species = A.cons), 
         prior = prior.m3, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.999, max_treedepth = 15))
# no divergent transitions
plot(m3) # good convergence, chains overlapping and stationary
summary(m3) # all Rhat 1, ESS >1000

# extract samples and save for post-processing
post_m3<- posterior_samples(m3)
save(post_m3, file="post_m3.robj")


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
get_prior(f_coalitions ~ Food_resource_defendable + f_philopatry + (1|Genus_species), data = d, family = bernoulli(), cov_ranef = list(Genus_species = A.cons))
prior.m4<- c(
  prior(normal(0,1), "Intercept"), 
  prior(normal(0,0.5), "b"), 
  prior(exponential(1), "sd"))

## run model
m4<- brm(f_coalitions ~ Food_resource_defendable + f_philopatry + (1|Genus_species), data = d, family = bernoulli(), cov_ranef = list(Genus_species = A.cons), 
         prior = prior.m4, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.999, max_treedepth = 15))
# no divergent transitions
plot(m4) # good convergence, chains overlapping and stationary
summary(m4) # all Rhat 1, ESS >1000

# extract samples and save for post-processing
post_m4<- posterior_samples(m4)
save(post_m4, file="post_m4.robj")


## robustness checks: each predictor on its own
m4a<- brm(f_coalitions ~ Food_resource_defendable + (1|Genus_species), data = d, family = bernoulli(), cov_ranef = list(Genus_species = A.cons), 
          prior = prior.m4, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.999, max_treedepth = 15))
# no divergent transitions
plot(m4a)
summary(m4a)
post_m4a<- posterior_samples(m4a)
save(post_m4a, file="post_m4a.robj")

m4b<- brm(f_coalitions ~ f_philopatry + (1|Genus_species), data = d, family = bernoulli(), cov_ranef = list(Genus_species = A.cons), 
          prior = prior.m4, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.999, max_treedepth = 15))
# no divergent transitions
plot(m4b)
summary(m4b)
post_m4b<- posterior_samples(m4b)
save(post_m4b, file="post_m4b.robj")


### Model 5: Male coalitions ~ Dimorphism, Philopatry ###

# transform predictors
d$SexDim.z<- (d$Sex_Dim-1)/sd(d$Sex_Dim) # mean = monomorphic, units = SD
d$m_philopatry<- 1
d$m_philopatry[d$Philopatry=="N" | d$Philopatry=="F"]<- 0

# re-code DV
d$m_coalitions<- 1
d$m_coalitions[d$coalitions=="Females"]<- 0

# examine variation
table(d$m_coalitions, d$m_philopatry)

d$m_philopatry<- as.factor(d$m_philopatry)

# same prior as for m4
m5<- brm(m_coalitions ~ SexDim.z + m_philopatry + (1|Genus_species), data = d, family = bernoulli(), cov_ranef = list(Genus_species = A.cons), 
         prior = prior.m4, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.999, max_treedepth = 15))
# no divergent transitions
plot(m5) # good convergence, chains overlapping and stationary
summary(m5) # all Rhat 1, ESS >1000

# extract samples and save for post-processing
post_m5<- posterior_samples(m5)
save(post_m5, file="post_m5.robj")

## robustness checks: each predictor on its own
m5a<- brm(m_coalitions ~ SexDim.z + (1|Genus_species), data = d, family = bernoulli(), cov_ranef = list(Genus_species = A.cons), 
          prior = prior.m4, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.999, max_treedepth = 15))
# no divergent transitions
plot(m5a)
summary(m5a)
post_m5a<- posterior_samples(m5a)
save(post_m5a, file="post_m5a.robj")

m5b<- brm(m_coalitions ~ m_philopatry + (1|Genus_species), data = d, family = bernoulli(), cov_ranef = list(Genus_species = A.cons), 
          prior = prior.m4, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.999, max_treedepth = 15))
# no divergent transitions
plot(m5b)
summary(m5b)
post_m5b<- posterior_samples(m5b)
save(post_m5b, file="post_m5b.robj")
