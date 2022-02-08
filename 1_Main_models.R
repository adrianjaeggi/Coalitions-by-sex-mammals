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
setwd("//files.iem.uzh.ch/Data/Institute/Human_Ecology/ajaegg/Private/My Documents/GitHub/sex-bias-cooperation")

## read files
d<- read.csv("sexbiascoalitions.csv", header=TRUE)
trees<- read.tree("TreeBlockSexBiasCooperation.tre")

# convert trees to covariance matrix (see https://cran.r-project.org/web/packages/brms/vignettes/brms_phylogenetics.html)
A<- list() 
for(i in 1:100){
  A[[i]]<- vcv.phylo(trees[[i]])
}


# create consensus tree
cons.tree<- consensus.edges(trees)
A.cons<- vcv.phylo(cons.tree)

d$coalitions<- relevel(d$coalitions, ref="Unbiased")
d$Philopatry<- relevel(d$Philopatry, ref="N")


  ### Model 1: Intercept only ###

## set prior (same as in previous paper on sex bias)
prior.m1<- c(
  prior(normal(0,5), "Intercept", dpar="muFemalebiased"), # reasonably narrow prior (favors values from ~ -10 to 10 in log-odds space, i.e. basically 0-1 in probability space) --> no need to sample infinite space beyond that
  prior(normal(0,5), "Intercept", dpar="muMalebiased"),
  prior(cauchy(0,1), "sd", dpar="muFemalebiased"), # pretty strong prior on variance components because phylo variance with huge CIs
  prior(cauchy(0,1), "sd", dpar="muMalebiased")
)

# simulate from prior
m1.prior<- brm(coalitions ~ 1 + (1|Genus_species), data = d, family = categorical(), cov_ranef = list(Genus_species = A.cons), 
            prior = prior.m1, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.90, max_treedepth = 15, stepsize = 0.01),
            sample_prior = "only")
post_m1.prior<- posterior_samples(m1.prior)
dens(logistic(post_m1.prior$b_muFemalebiased_Intercept))
# this prior actually puts most weight on 0 and 1 --> choose more reasonable range
dens(post_m1.prior$sd_Genus_species__muFemalebiased_Intercept)
summary(post_m1.prior$sd_Genus_species__muFemalebiased_Intercept)
# vast majority of samples near 0, but some insanely large

## adjusted prior
prior.m1<- c(
  prior(normal(0,1), "Intercept", dpar="muFemalebiased"), # reasonably narrow prior (favors values from ~ -10 to 10 in log-odds space, i.e. basically 0-1 in probability space) --> no need to sample infinite space beyond that
  prior(normal(0,1), "Intercept", dpar="muMalebiased"),
  prior(cauchy(0,0.1), "sd", dpar="muFemalebiased"), # pretty strong prior on variance components because phylo variance with huge CIs
  prior(cauchy(0,0.1), "sd", dpar="muMalebiased")
)

# simulate from prior
m1.prior<- brm(coalitions ~ 1 + (1|Genus_species), data = d, family = categorical(), cov_ranef = list(Genus_species = A.cons), 
               prior = prior.m1, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.90, max_treedepth = 15, stepsize = 0.01),
               sample_prior = "only")
post_m1.prior<- posterior_samples(m1.prior)
dens(logistic(post_m1.prior$b_muFemalebiased_Intercept))
# much more reasonable: most weight around 0.5 and very little at 0 or 1
dens(post_m1.prior$sd_Genus_species__muFemalebiased_Intercept)
summary(post_m1.prior$sd_Genus_species__muFemalebiased_Intercept)
# now vast majority of samples <1, some still quite large
# this prior is fine

m1<- brm(coalitions ~ 1 + (1|Genus_species), data = d, family = categorical(), cov_ranef = list(Genus_species = A.cons), 
               prior = prior.m1, chains = 4, cores = 4, iter = 4000, warmup = 3000, control = list(adapt_delta = 0.999, max_treedepth = 15, stepsize = 0.0001))
# no divergent transitions
plot(m1) # good convergence, chains overlapping and stationary; some extreme values for sd's
summary(m1) # all Rhat 1, ESS >1000
# virtually same probability for male-bias and female-bias

# extract samples and save for post-processing
post_m1<- posterior_samples(m1)
save(post_m1, file="post_m1.robj")


    ### Model 2: Sex-segregation ###

get_prior(coalitions ~ SexComposition + (1|Genus_species), data = d, family = categorical(), cov_ranef = list(Genus_species = A.cons))
## adjust prior
prior.m2<- c(
  prior(normal(0,1), "Intercept", dpar="muFemalebiased"), # reasonably narrow prior (favors values from ~ -10 to 10 in log-odds space, i.e. basically 0-1 in probability space) --> no need to sample infinite space beyond that
  prior(normal(0,1), "Intercept", dpar="muMalebiased"),
  prior(normal(0,1), "b", coef = "SexCompositionsegregated", dpar="muFemalebiased"),
  prior(normal(0,1), "b", coef = "SexCompositionsegregated", dpar="muMalebiased"),
  prior(cauchy(0,0.1), "sd", dpar="muFemalebiased"), # pretty strong prior on variance components because phylo variance with huge CIs
  prior(cauchy(0,0.1), "sd", dpar="muMalebiased")
)

# simulate from prior
m2.prior<- brm(coalitions ~ SexComposition + (1|Genus_species), data = d, family = categorical(), cov_ranef = list(Genus_species = A.cons), 
               prior = prior.m2, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.90, max_treedepth = 15, stepsize = 0.01),
               sample_prior = "only")
post_m2.prior<- posterior_samples(m2.prior)
dens(exp(post_m2.prior$b_muFemalebiased_SexCompositionsegregated)) # check prior for categorical coefficients on odd's ratio scale 
summary(exp(post_m2.prior$b_muFemalebiased_SexCompositionsegregated)) # median = 1, but mean is 1.7, and very fat tail
# try narrower prior on coefficient

prior.m2<- c(
  prior(normal(0,1), "Intercept", dpar="muFemalebiased"), # reasonably narrow prior (favors values from ~ -10 to 10 in log-odds space, i.e. basically 0-1 in probability space) --> no need to sample infinite space beyond that
  prior(normal(0,1), "Intercept", dpar="muMalebiased"),
  prior(normal(0,0.5), "b", coef = "SexCompositionsegregated", dpar="muFemalebiased"),
  prior(normal(0,0.5), "b", coef = "SexCompositionsegregated", dpar="muMalebiased"),
  prior(cauchy(0,0.1), "sd", dpar="muFemalebiased"), # pretty strong prior on variance components because phylo variance with huge CIs
  prior(cauchy(0,0.1), "sd", dpar="muMalebiased")
)

# simulate from prior
m2.prior<- brm(coalitions ~ SexComposition + (1|Genus_species), data = d, family = categorical(), cov_ranef = list(Genus_species = A.cons), 
               prior = prior.m2, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.90, max_treedepth = 15, stepsize = 0.01),
               sample_prior = "only")
post_m2.prior<- posterior_samples(m2.prior)
dens(exp(post_m2.prior$b_muFemalebiased_SexCompositionsegregated)) # check prior for categorical coefficients on odd's ratio scale 
summary(exp(post_m2.prior$b_muFemalebiased_SexCompositionsegregated)) # median = 1, mean = 1.1, more reasonable tail (max = 5) --> go with this

m2<- brm(coalitions ~ SexComposition + (1|Genus_species), data = d, family = categorical(), cov_ranef = list(Genus_species = A.cons), 
               prior = prior.m2, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.999, max_treedepth = 15, stepsize = 0.01))
# adjust sampler settings until no divergent transitions
plot(m2) # good convergence, chains overlapping and stationary; some extreme values for sd's
summary(m2) # all Rhat 1, ESS >1000
# no effect of sex segregation -> safe to ignore in future models

# extract samples and save for post-processing
post_m2<- posterior_samples(m2)
save(post_m2, file="post_m2.robj")


    ### Model 3: Primates vs others ###

# create dummy
d$Primate<- 0
d$Primate[d$Order=="Primata"]<- 1

# adjust prior
prior.m3<- c(
  prior(normal(0,1), "Intercept", dpar="muFemalebiased"), # reasonably narrow prior (favors values from ~ -10 to 10 in log-odds space, i.e. basically 0-1 in probability space) --> no need to sample infinite space beyond that
  prior(normal(0,1), "Intercept", dpar="muMalebiased"),
  prior(normal(0,0.5), "b", coef = "Primate", dpar="muFemalebiased"),
  prior(normal(0,0.5), "b", coef = "Primate", dpar="muMalebiased"),
  prior(cauchy(0,0.1), "sd", dpar="muFemalebiased"), # pretty strong prior on variance components because phylo variance with huge CIs
  prior(cauchy(0,0.1), "sd", dpar="muMalebiased")
)

m3<- brm(coalitions ~ Primate + (1|Genus_species), data = d, family = categorical(), cov_ranef = list(Genus_species = A.cons), 
         prior = prior.m3, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.999, max_treedepth = 15, stepsize = 0.01))
# adjust sampler settings until no divergent transitions
plot(m3) # good convergence, chains overlapping and stationary; some extreme values for sd's
summary(m3) # all Rhat 1, ESS >1000
# slightly lower probabilities in primates, but very uncertain

# extract samples and save for post-processing
post_m3<- posterior_samples(m3)
save(post_m3, file="post_m3.robj")


    ### Model 4: Main predictors ###

# transform predictors
d$SexDim.z<- (d$Sex_Dim-1)/sd(d$Sex_Dim) # mean = monomorphic, units = SD

# adjust prior
get_prior(coalitions ~ Food_resource_defendable + Philopatry + SexDim.z + (1|Genus_species), data = d, family = categorical(), cov_ranef = list(Genus_species = A.cons))
prior.m4<- c(
  prior(normal(0,1), "Intercept", dpar="muFemalebiased"), # reasonably narrow prior (favors values from ~ -10 to 10 in log-odds space, i.e. basically 0-1 in probability space) --> no need to sample infinite space beyond that
  prior(normal(0,1), "Intercept", dpar="muMalebiased"),
  prior(normal(0,0.5), "b", dpar="muFemalebiased"), # same prior for all slopes
  prior(normal(0,0.5), "b", dpar="muMalebiased"),
  prior(cauchy(0,0.1), "sd", dpar="muFemalebiased"), # pretty strong prior on variance components because phylo variance with huge CIs
  prior(cauchy(0,0.1), "sd", dpar="muMalebiased")
)

# simulate from prior to check slopes on non-categorical predictors
m4.prior<- brm(coalitions ~ Food_resource_defendable + Philopatry + SexDim.z + (1|Genus_species), data = d, family = categorical(), cov_ranef = list(Genus_species = A.cons), 
               prior = prior.m4, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.999, max_treedepth = 15, stepsize = 0.01),
               sample_prior = "only")

post_m4.prior<- posterior_samples(m4.prior)
plot(NULL, xlim=c(-3,3), ylim=c(-5,5)) # -5/5 is roughly 0/1 on logistic scale
for(i in sample(1:nrow(post_m4.prior), 100)){
  abline(post_m4.prior[i,1],post_m4.prior[i,7], col=col.alpha("lightblue", 0.2), lwd=0.5)
}
# reasonable prior, most associations not too strong

m4<- brm(coalitions ~ Food_resource_defendable + Philopatry + SexDim.z + (1|Genus_species), data = d, family = categorical(), cov_ranef = list(Genus_species = A.cons), 
               prior = prior.m4, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.999, max_treedepth = 15, stepsize = 0.01))
# adjust sampler settings until no divergent transitions --> 1 divergent transition, fine for now
plot(m4) # good convergence, chains overlapping and stationary; some extreme values for sd's
summary(m4) # all Rhat 1, ESS >1000
# strongest effect is of sexual dimorphism

# extract samples and save for post-processing
post_m4<- posterior_samples(m4)
save(post_m4, file="post_m4.robj")


  ## robustness checks: each predictor on its own
m4a<- brm(coalitions ~ Food_resource_defendable + (1|Genus_species), data = d, family = categorical(), cov_ranef = list(Genus_species = A.cons), 
         prior = prior.m4, chains = 4, cores = 4, iter = 5000, warmup = 3000, control = list(adapt_delta = 0.999, max_treedepth = 15, stepsize = 0.01))
# 40 divergent transitions -> increase warmup, iter -> no more divergent transitions
plot(m4a)
summary(m4a)
# all good, no strong associations
post_m4a<- posterior_samples(m4a)
save(post_m4a, file="post_m4a.robj")

m4b<- brm(coalitions ~ Philopatry + (1|Genus_species), data = d, family = categorical(), cov_ranef = list(Genus_species = A.cons), 
         prior = prior.m4, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.999, max_treedepth = 15, stepsize = 0.01))
# no divergent transitions
plot(m4b)
summary(m4b)
# all good, no strong effects, though trend to more female coalitions with female philopatry
post_m4b<- posterior_samples(m4b)
save(post_m4b, file="post_m4b.robj")

m4c<- brm(coalitions ~ SexDim.z + (1|Genus_species), data = d, family = categorical(), cov_ranef = list(Genus_species = A.cons), 
         prior = prior.m4, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.999, max_treedepth = 15, stepsize = 0.01))
# no divergent transitions
plot(m4c)
summary(m4c)
# all good, still equally strong effect of dimorphism on male bias
post_m4c<- posterior_samples(m4c)
save(post_m4c, file="post_m4c.robj")



