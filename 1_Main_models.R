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

d$coalitions<- relevel(d$coalitions, ref="Both")
d$Philopatry<- relevel(d$Philopatry, ref="N")


  ### Model 1: Intercept only ###

## set prior (same as in previous paper on sex bias)
prior.m1<- c(
  prior(normal(0,5), "Intercept", dpar="muFemales"), # reasonably narrow prior (favors values from ~ -10 to 10 in log-odds space, i.e. basically 0-1 in probability space) --> no need to sample infinite space beyond that
  prior(normal(0,5), "Intercept", dpar="muMales"),
  prior(cauchy(0,1), "sd", dpar="muFemales"), # pretty strong prior on variance components because phylo variance with huge CIs
  prior(cauchy(0,1), "sd", dpar="muMales")
)

# simulate from prior
m1.prior<- brm(coalitions ~ 1 + (1|Genus_species), data = d, family = categorical(), cov_ranef = list(Genus_species = A.cons), 
            prior = prior.m1, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.90, max_treedepth = 15, stepsize = 0.01),
            sample_prior = "only")
post_m1.prior<- posterior_samples(m1.prior)
dens(logistic(post_m1.prior$b_muFemales_Intercept))
# this prior actually puts most weight on 0 and 1 --> choose more reasonable range
dens(post_m1.prior$sd_Genus_species__muFemales_Intercept)
summary(post_m1.prior$sd_Genus_species__muFemales_Intercept)
# vast majority of samples near 0, but some insanely large

## adjusted prior
prior.m1<- c(
  prior(normal(0,1), "Intercept", dpar="muFemales"), # reasonably narrow prior (favors values from ~ -10 to 10 in log-odds space, i.e. basically 0-1 in probability space) --> no need to sample infinite space beyond that
  prior(normal(0,1), "Intercept", dpar="muMales"),
  prior(cauchy(0,0.1), "sd", dpar="muFemales"), # pretty strong prior on variance components because phylo variance with huge CIs
  prior(cauchy(0,0.1), "sd", dpar="muMales")
)

# simulate from prior
m1.prior<- brm(coalitions ~ 1 + (1|Genus_species), data = d, family = categorical(), cov_ranef = list(Genus_species = A.cons), 
               prior = prior.m1, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.90, max_treedepth = 15, stepsize = 0.01),
               sample_prior = "only")
post_m1.prior<- posterior_samples(m1.prior)
dens(logistic(post_m1.prior$b_muFemales_Intercept))
# much more reasonable: most weight around 0.5 and very little at 0 or 1
dens(post_m1.prior$sd_Genus_species__muFemales_Intercept)
summary(post_m1.prior$sd_Genus_species__muFemales_Intercept)
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
  prior(normal(0,1), "Intercept", dpar="muFemales"), # reasonably narrow prior (favors values from ~ -10 to 10 in log-odds space, i.e. basically 0-1 in probability space) --> no need to sample infinite space beyond that
  prior(normal(0,1), "Intercept", dpar="muMales"),
  prior(normal(0,1), "b", coef = "SexCompositionsegregated", dpar="muFemales"),
  prior(normal(0,1), "b", coef = "SexCompositionsegregated", dpar="muMales"),
  prior(cauchy(0,0.1), "sd", dpar="muFemales"), # pretty strong prior on variance components because phylo variance with huge CIs
  prior(cauchy(0,0.1), "sd", dpar="muMales")
)

# simulate from prior
m2.prior<- brm(coalitions ~ SexComposition + (1|Genus_species), data = d, family = categorical(), cov_ranef = list(Genus_species = A.cons), 
               prior = prior.m2, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.90, max_treedepth = 15, stepsize = 0.01),
               sample_prior = "only")
post_m2.prior<- posterior_samples(m2.prior)
dens(exp(post_m2.prior$b_muFemales_SexCompositionsegregated)) # check prior for categorical coefficients on odd's ratio scale 
summary(exp(post_m2.prior$b_muFemales_SexCompositionsegregated)) # median = 1, but mean is 1.7, and very fat tail
# try narrower prior on coefficient

prior.m2<- c(
  prior(normal(0,1), "Intercept", dpar="muFemales"), # reasonably narrow prior (favors values from ~ -10 to 10 in log-odds space, i.e. basically 0-1 in probability space) --> no need to sample infinite space beyond that
  prior(normal(0,1), "Intercept", dpar="muMales"),
  prior(normal(0,0.5), "b", coef = "SexCompositionsegregated", dpar="muFemales"),
  prior(normal(0,0.5), "b", coef = "SexCompositionsegregated", dpar="muMales"),
  prior(cauchy(0,0.1), "sd", dpar="muFemales"), # pretty strong prior on variance components because phylo variance with huge CIs
  prior(cauchy(0,0.1), "sd", dpar="muMales")
)

# simulate from prior
m2.prior<- brm(coalitions ~ SexComposition + (1|Genus_species), data = d, family = categorical(), cov_ranef = list(Genus_species = A.cons), 
               prior = prior.m2, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.90, max_treedepth = 15, stepsize = 0.01),
               sample_prior = "only")
post_m2.prior<- posterior_samples(m2.prior)
dens(exp(post_m2.prior$b_muFemales_SexCompositionsegregated)) # check prior for categorical coefficients on odd's ratio scale 
summary(exp(post_m2.prior$b_muFemales_SexCompositionsegregated)) # median = 1, mean = 1.1, more reasonable tail (max = 5) --> go with this

m2<- brm(coalitions ~ SexComposition + (1|Genus_species), data = d, family = categorical(), cov_ranef = list(Genus_species = A.cons), 
               prior = prior.m2, chains = 4, cores = 4, iter = 4000, warmup = 3000, control = list(adapt_delta = 0.999, max_treedepth = 15, stepsize = 0.01))
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
  prior(normal(0,1), "Intercept", dpar="muFemales"), # reasonably narrow prior (favors values from ~ -10 to 10 in log-odds space, i.e. basically 0-1 in probability space) --> no need to sample infinite space beyond that
  prior(normal(0,1), "Intercept", dpar="muMales"),
  prior(normal(0,0.5), "b", coef = "Primate", dpar="muFemales"),
  prior(normal(0,0.5), "b", coef = "Primate", dpar="muMales"),
  prior(cauchy(0,0.1), "sd", dpar="muFemales"), # pretty strong prior on variance components because phylo variance with huge CIs
  prior(cauchy(0,0.1), "sd", dpar="muMales")
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
  prior(cauchy(0,0.1), "sd"))

# simulate from prior to check slopes on non-categorical predictors
m4.prior<- brm(f_coalitions ~ Food_resource_defendable + f_philopatry + (1|Genus_species), data = d, family = bernoulli(), cov_ranef = list(Genus_species = A.cons), 
               prior = prior.m4, chains = 4, cores = 4, iter = 4000, warmup = 2000, control = list(adapt_delta = 0.999, max_treedepth = 15, stepsize = 0.0001), inits=0,
               sample_prior = "only")
post_m4.prior<- posterior_samples(m4.prior)
plot(NULL, xlim=c(-3,3), ylim=c(-5,5)) # -5/5 is roughly 0/1 on logistic scale
for(i in sample(1:nrow(post_m4.prior), 100)){
  abline(post_m4.prior[i,1],post_m4.prior[i,7], col=col.alpha("lightblue", 0.2), lwd=0.5)
}
# reasonable prior, most associations not too strong

## run model
m4<- brm(f_coalitions ~ Food_resource_defendable + f_philopatry + (1|Genus_species), data = d, family = bernoulli(), cov_ranef = list(Genus_species = A.cons), 
         prior = prior.m4, chains = 4, cores = 4, iter = 5000, warmup = 3000, control = list(adapt_delta = 0.999, max_treedepth = 15, stepsize = 0.0001), inits=0)
# 9 divergent transitions, adjusting settings doesn't help much --> re-parameterize
d$f_philopatry<- as.factor(d$f_philopatry)
d$Food_resource_defendable<- as.factor(d$Food_resource_defendable)
get_prior(f_coalitions ~ -1 + Food_resource_defendable * f_philopatry + (1|Genus_species), data = d, family = bernoulli(), cov_ranef = list(Genus_species = A.cons))
prior.m4<- c(
  prior(normal(0,0.5), "b"), 
  prior(cauchy(0,0.1), "sd"))
m4<- brm(f_coalitions ~ -1 + Food_resource_defendable * f_philopatry + (1|Genus_species), data = d, family = bernoulli(), cov_ranef = list(Genus_species = A.cons), 
         prior = prior.m4, chains = 4, cores = 4, iter = 5000, warmup = 3000, control = list(adapt_delta = 0.999, max_treedepth = 15, stepsize = 0.0001), inits=0)
# 3 divergent transitions, that's better
plot(m4) # good convergence, chains overlapping and stationary; some extreme values for sd's
summary(m4) # all Rhat 1, ESS >1000
# no effect of defensibility, philopatry

# extract samples and save for post-processing
post_m4<- posterior_samples(m4)
save(post_m4, file="post_m4.robj")


## robustness checks: each predictor on its own
m4a<- brm(f_coalitions ~ -1 + Food_resource_defendable + (1|Genus_species), data = d, family = bernoulli(), cov_ranef = list(Genus_species = A.cons), 
          prior = prior.m4, chains = 4, cores = 4, iter = 5000, warmup = 3000, control = list(adapt_delta = 0.999, max_treedepth = 15, stepsize = 0.0001), inits=0)
# 12 divergent transitions
plot(m4a)
summary(m4a)
# all good, no strong association
post_m4a<- posterior_samples(m4a)
save(post_m4a, file="post_m4a.robj")

m4b<- brm(f_coalitions ~ -1 + f_philopatry + (1|Genus_species), data = d, family = bernoulli(), cov_ranef = list(Genus_species = A.cons), 
          prior = prior.m4, chains = 4, cores = 4, iter = 5000, warmup = 3000, control = list(adapt_delta = 0.999, max_treedepth = 15, stepsize = 0.0001), inits=0)
# 5 divergent transitions
plot(m4b)
summary(m4b)
# all good, no strong effects
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
m5<- brm(m_coalitions ~ -1 + SexDim.z + m_philopatry + (1|Genus_species), data = d, family = bernoulli(), cov_ranef = list(Genus_species = A.cons), 
         prior = prior.m4, chains = 4, cores = 4, iter = 5000, warmup = 3000, control = list(adapt_delta = 0.999, max_treedepth = 15, stepsize = 0.0001), inits=0)
# 2 divergent transitions
plot(m5) # good convergence, chains overlapping and stationary; some extreme values for sd's
summary(m5) # all Rhat 1, ESS >1000
# no effect of dimorphism, philopatry

# extract samples and save for post-processing
post_m5<- posterior_samples(m5)
save(post_m5, file="post_m5.robj")

prior.m5<- c(
  prior(normal(0,1), "Intercept"), 
  prior(normal(0,0.5), "b"), 
  prior(cauchy(0,0.1), "sd"))

## robustness checks: each predictor on its own
m5a<- brm(m_coalitions ~ SexDim.z + (1|Genus_species), data = d, family = bernoulli(), cov_ranef = list(Genus_species = A.cons), 
          prior = prior.m5, chains = 4, cores = 4, iter = 5000, warmup = 3000, control = list(adapt_delta = 0.999, max_treedepth = 15, stepsize = 0.0001), inits=0)
# 1 divergent transitions
plot(m5a)
summary(m5a)
# all good, weak association
post_m5a<- posterior_samples(m5a)
save(post_m5a, file="post_m5a.robj")

m5b<- brm(m_coalitions ~ -1 + m_philopatry + (1|Genus_species), data = d, family = bernoulli(), cov_ranef = list(Genus_species = A.cons), 
          prior = prior.m4, chains = 4, cores = 4, iter = 5000, warmup = 3000, control = list(adapt_delta = 0.999, max_treedepth = 15, stepsize = 0.0001), inits=0)
# no divergent transitions
plot(m5b)
summary(m5b)
# all good, no strong effects
post_m5b<- posterior_samples(m5b)
save(post_m5b, file="post_m5b.robj")

