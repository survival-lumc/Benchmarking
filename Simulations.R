library(survival)
library(progress)
library(mvtnorm)
library(dplyr)
options(dplyr.summarise.inform = FALSE) # suppress summarise info
library(dynpred) # only if progress=2
library(tictoc)
library(prodlim)
library(geepack)

source("Function_definitions.R")

eps <- 1e-10

simsettings <- list(centers = list(M = 300, # number of centers, was 302
                                   mean = 200, # average center size
                                   var = 22500), # variance of center sizes, should be >= mean
                    event = list(shape = list(mean = log(0.94), # log Weibull shape
                                              var = eps),
                                 rate = list(mean = log(0.032), # log Weibull rate
                                             var = eps), # frailty variance
                                 cor = 0,
                                 tau = 12), # horizon time point
                    fup = list(shape = list(mean = 0.4, # log Weibull shape
                                            var = 0.24^2),
                               rate = list(mean = -4.8, # log Weibull rate
                                           var = 1.72^2),
                               cor = -0.87),
                    x = list(varw = 0.224, # within centers (residual) variance
                             varb = 0.056, # between centers variance
                             beta = 1), # effect (log hazard ratio)
                    seed = 20210203,
                    nsim = 50
)
simsettings_base <- simsettings

#
# Commented out, because this will take a VERY long time
#
# tic("Simulations (50 replications) from base scenario")
# res <- dosim(simsettings, progress=2)
# toc()
# hist(res$Z)
# mean(res$Z)
# sd(res$Z)
# table(res$Performance)

#
# Pseudo-observation method of Logan
#
# Also commented out because it takes a long time (though considerably
# less then above)
# res <- dosim.Logan(simsettings, progress=1)
# table(res$under)
# table(res$over)

#
# Original simulation was done by running five of these simulations
# with 10 replications each in parallel, using different random seeds
#

#
# Here are the other scenarios studied in the paper
#

# Base, same fup
simsettings <- simsettings_base
simsettings$fup$shape$var <- eps
simsettings$fup$rate$var <- eps

# Fewer centers
simsettings <- simsettings_base
simsettings$centers$M <- 30
simsettings$nsim <- 500

# Fewer patients
simsettings <- simsettings_base
simsettings$centers$mean <- 20
simsettings$centers$var <- 225
simsettings$nsim <- 500

# Non-PH
simsettings <- simsettings_base
simsettings$event$shape$var <- 0.15

# Small frailty
simsettings <- simsettings_base
simsettings$event$rate$var <- 0.15

# Large frailty
simsettings <- simsettings_base
simsettings$event$rate$var <- 0.3

