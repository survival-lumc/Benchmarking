my.rweibull <- function(n, a, b) # shape a and rate b
{
  return((-log(runif(n))/b)^(1/a))
}

rep.row <- function(x, n) {
  matrix(rep(x, each=n), nrow=n)
}

gendata <- function(simsettings)
{
  # This function generates data according to procedures described in
  # Benchmarking survival outcomes: A funnel plot for survival data
  #
  # Unpack simulation settings
  M <- simsettings$centers$M
  mu <- simsettings$centers$mean
  size <- ifelse(simsettings$centers$var == mu,
                 1/eps,
                 mu^2 / (simsettings$centers$var - mu))
  ni <- rnbinom(M, mu=mu, size=size)
  ni[ni<5] <- 5 # set minimum size of 5
  N <- sum(ni)
  me1 <- simsettings$event$shape$mean
  vare1 <- simsettings$event$shape$var
  me2 <- simsettings$event$rate$mean
  vare2 <- simsettings$event$rate$var
  core <- simsettings$event$cor
  cove <- core * sqrt(vare1 * vare2)
  mc1 <- simsettings$fup$shape$mean
  varc1 <- simsettings$fup$shape$var
  mc2 <- simsettings$fup$rate$mean
  varc2 <- simsettings$fup$rate$var
  corc <- simsettings$fup$cor
  covc <- corc * sqrt(varc1 * varc2)
  varw <- simsettings$x$varw
  varb <- simsettings$x$varb
  beta <- simsettings$x$beta

  ## Generate Weibull parameters for baseline event distributions in the different centers
  wepars <- matrix(NA, M, 2)
  sige <- matrix(c(vare1, cove, cove, vare2), 2, 2)
  wepars <- rmvnorm(M, mean = c(me1, me2), sigma = sige)
  wepars <- exp(wepars)
  ## Generate Weibull parameters for censoring distributions in the different centers
  wcpars <- matrix(NA, M, 2)
  sigc <- matrix(c(varc1, covc, covc, varc2), 2, 2)
  wcpars <- rmvnorm(M, mean = c(mc1, mc2), sigma = sigc)
  wcpars <- exp(wcpars)
  # First mean per center
  xc <- rnorm(M, mean=0, sd=sqrt(varb))
  # Reserve room for generated event and censoring time points and x
  tt <- cc <- x <- aa <- bb <- rep(NA, N)
  start <- 1
  for (m in 1:M) {
    n <- ni[m]
    end <- start-1 + n
    whm <- (start:end)
    # Generate x
    x[whm] <- rnorm(n, mean=xc, sd=sqrt(varw))
    # Generate time-to-event
    am <- wepars[m, 1]
    bm <- wepars[m, 2]
    tt[whm] <- (round(my.rweibull(n, a = am, b = bm * exp(beta * x[whm])) * 500) + 0.5)/ 500
    # Generate censoring
    am <- wcpars[m, 1]
    bm <- wcpars[m, 2]
    cc[whm] <- (round(my.rweibull(n, a=am, b=bm) * 500) + 1)/ 500
    # Save censoring parameters
    aa[whm] <- am
    bb[whm] <- bm
    start <- end + 1
  }
  dfr <- data.frame(ID=rep(1:M, ni), time=pmin(tt, cc), status=as.numeric(tt<=cc), x=x,
                    a=aa, b=bb)
  return(dfr)
}

gendata_nonPH <- function(simsettings)
{
  # This function generates data according to procedures described in
  # Benchmarking survival outcomes: A funnel plot for survival data
  # This is the non-proportional hazards setting
  #
  # Unpack simulation settings
  M <- simsettings$centers$M
  mu <- simsettings$centers$mean
  size <- ifelse(simsettings$centers$var == mu,
                 1/eps,
                 mu^2 / (simsettings$centers$var - mu))
  ni <- rnbinom(M, mu=mu, size=size)
  ni[ni<5] <- 5 # set minimum size of 5
  N <- sum(ni)
  me1 <- simsettings$event$shape$mean
  vare1 <- simsettings$event$shape$var
  me2 <- simsettings$event$rate$mean
  vare2 <- simsettings$event$rate$var
  core <- simsettings$event$cor
  cove <- core * sqrt(vare1 * vare2)
  mc1 <- simsettings$fup$shape$mean
  varc1 <- simsettings$fup$shape$var
  mc2 <- simsettings$fup$rate$mean
  varc2 <- simsettings$fup$rate$var
  corc <- simsettings$fup$cor
  covc <- corc * sqrt(varc1 * varc2)
  varw <- simsettings$x$varw
  varb <- simsettings$x$varb
  beta <- simsettings$x$beta
  
  ## Generate Weibull parameters for baseline event distributions in the different centers
  wepars <- matrix(NA, M, 2)
  sige <- matrix(c(vare1, cove, cove, vare2), 2, 2)
  wepars <- rmvnorm(M, mean = c(me1, me2), sigma = sige)
  wepars <- exp(wepars)
  # Reset rate parameter to ensure that the 12 months probability equals that for means
  S12 <- exp(-0.032 * 12^0.94)
  wepars[, 1] <- log(0.032 * 12^0.94 / wepars[, 2]) / log(12)
  
  ## Generate Weibull parameters for censoring distributions in the different centers
  wcpars <- matrix(NA, M, 2)
  sigc <- matrix(c(varc1, covc, covc, varc2), 2, 2)
  wcpars <- rmvnorm(M, mean = c(mc1, mc2), sigma = sigc)
  wcpars <- exp(wcpars)
  # First mean per center
  xc <- rnorm(M, mean=0, sd=sqrt(varb))
  # Reserve room for generated event and censoring time points and x
  tt <- cc <- x <- rep(NA, N)
  start <- 1
  for (m in 1:M) {
    n <- ni[m]
    end <- start-1 + n
    whm <- (start:end)
    # Generate x
    x[whm] <- rnorm(n, mean=xc, sd=sqrt(varw))
    # Generate time-to-event
    am <- wepars[m, 1]
    bm <- wepars[m, 2]
    tt[whm] <- (round(my.rweibull(n, a = am, b = bm * exp(beta * x[whm])) * 500) + 0.5)/ 500
    # Generate censoring
    am <- wcpars[m, 1]
    bm <- wcpars[m, 2]
    cc[whm] <- (round(my.rweibull(n, a=am, b=bm) * 500) + 1)/ 500
    start <- end + 1
  }
  dfr <- data.frame(ID=rep(1:M, ni), time=pmin(tt, cc), status=as.numeric(tt<=cc), x=x)
  return(dfr)
}

dosim <- function(simsettings, progress=1) {
  # This function performs simulations according to procedures described in
  # Benchmarking survival outcomes: A funnel plot for survival data
  #
  # Unpack simulation settings
  seed <- simsettings$seed
  nsim <- simsettings$nsim
  tau <- simsettings$event$tau
  
  set.seed(seed)
  all.agg <- NULL # this will contain all results
  if (progress==1) pb <- progress_bar$new(total = nsim)
  for (sim in 1:nsim) {
    if (progress==1) pb$tick()
    if (progress==2) deb(sim, method="cat")
    gd <- gendata(simsettings)
    gd$statustau <- gd$status
    gd$statustau[gd$time > tau] <- 0
    gd$timetau <- pmin(gd$time, tau)
    gd$pij <- NA
    
    cev <- coxph(Surv(timetau, statustau) ~ x, data=gd)
    sfev <- survfit(cev, newdata=gd[1, ]) # just for event time points
    tev <- sfev$time

    # Extract censoring parameters
    first.center <- which(!duplicated(gd$ID))
    am <- gd$a[first.center]
    bm <- gd$b[first.center]
    
    M <- max(gd$ID)
    ni <- table(gd$ID)
    N <- sum(ni)
    # Calculation of pij's
    if (progress==2) pb <- progress_bar$new(total = M)
    for (i in 1:M) {
      if (progress==2) pb$tick()
      whi <- which(gd$ID==i)
      gdi <- gd[whi, ]
      # Censoring distribution
      sfc <- survfit(Surv(time, status==0) ~ 1, data=gdi)
      tc <- sfc$time
      tt <- sort(unique(c(tc[tc<=tau], tev)))
      for (j in 1:ni[i]) {
        sfev <- survfit(cev, newdata=gdi[j, ])
        Sijt <- summary(sfev, times=tt, extend=TRUE)$surv
        fijt <- -diff(c(1, Sijt))
        Git <- summary(sfc, times=tt, extend=TRUE)$surv
        gd$pij[whi][j] <- sum(fijt * Git)
      }
    }
    
    gd <- as_tibble(gd)
    p0 <- mean(gd$statustau)
    agg <- gd %>% 
      group_by(ID) %>% 
      summarise(n=n(), Observed=sum(statustau), Expected=sum(pij), Variance=sum(pij*(1-pij))) %>%
      mutate(OE=Observed/Expected) %>% 
      mutate(Z=(Observed - Expected)/sqrt(Variance)) %>%
      mutate(best_half=Z<quantile(Z,0.5)) %>%
      mutate(NNN=(Expected^2/Variance)*(1-p0)/p0)
    
    crit1 <- qnorm(0.975)
    crit2 <- qnorm(1-0.025/nrow(agg))
    agg$Performance <- 1
    agg$Performance[agg$Observed < agg$Expected + crit2*sqrt(agg$Variance)] <- 2
    agg$Performance[agg$Observed < agg$Expected + crit1*sqrt(agg$Variance)] <- 3
    agg$Performance[agg$Observed < agg$Expected - crit1*sqrt(agg$Variance)] <- 4
    agg$Performance[agg$Observed < agg$Expected - crit2*sqrt(agg$Variance)] <- 5
    agg$Performance <- factor(agg$Performance, levels=1:5,
                              labels=c("Clearly worse than average", "Worse than average", "Within range",
                                       "Better than average", "Clearly better than average"))
    
    agg$OEplot <- agg$OE
    agg$OEplot[agg$OEplot>3] <- 3
    agg$toohigh <- 0
    agg$toohigh[agg$OE>3] <- 1
    agg$toohigh <- factor(agg$toohigh)
    agg$a <- am
    agg$b <- bm
    agg$sim <- sim
    all.agg <- rbind(all.agg, agg)
  }
  attr(all.agg, "simsettings") <- simsettings
  return(all.agg)
}

dosim_nonPH <- function(simsettings, progress=1) {
  # This function performs simulations according to procedures described in
  # Benchmarking survival outcomes: A funnel plot for survival data;
  # This is the non-proportional hazards setting, only difference actually
  # is that gendata_nonPH is used instead of gendata
  #
  # Unpack simsettings
  seed <- simsettings$seed
  nsim <- simsettings$nsim
  tau <- simsettings$event$tau
  
  set.seed(seed)
  all.agg <- NULL # this will contain all results
  if (progress==1) pb <- progress_bar$new(total = nsim)
  for (sim in 1:nsim) {
    if (progress==1) pb$tick()
    if (progress==2) deb(sim, method="cat")
    gd <- gendata_nonPH(simsettings)
    gd$statustau <- gd$status
    gd$statustau[gd$time > tau] <- 0
    gd$timetau <- pmin(gd$time, tau)
    gd$pij <- NA
    
    cev <- coxph(Surv(timetau, statustau) ~ x, data=gd)
    sfev <- survfit(cev, newdata=gd[1, ]) # just for event time points
    tev <- sfev$time

    M <- max(gd$ID)
    ni <- table(gd$ID)
    N <- sum(ni)
    # Calculation of pij's
    if (progress==2) pb <- progress_bar$new(total = M)
    for (i in 1:M) {
      if (progress==2) pb$tick()
      whi <- which(gd$ID==i)
      gdi <- gd[whi, ]
      # Censoring distribution
      sfc <- survfit(Surv(time, status==0) ~ 1, data=gdi)
      tc <- sfc$time
      tt <- sort(unique(c(tc[tc<=tau], tev)))
      for (j in 1:ni[i]) {
        sfev <- survfit(cev, newdata=gdi[j, ])
        Sijt <- summary(sfev, times=tt, extend=TRUE)$surv
        fijt <- -diff(c(1, Sijt))
        Git <- summary(sfc, times=tt, extend=TRUE)$surv
        gd$pij[whi][j] <- sum(fijt * Git)
      }
    }
    
    gd <- as_tibble(gd)
    p0 <- mean(gd$statustau)
    agg <- gd %>% 
      group_by(ID) %>% 
      summarise(n=n(), Observed=sum(statustau), Expected=sum(pij), Variance=sum(pij*(1-pij))) %>%
      mutate(OE=Observed/Expected) %>% 
      mutate(Z=(Observed - Expected)/sqrt(Variance)) %>%
      mutate(best_half=Z<quantile(Z,0.5)) %>%
      mutate(NNN=(Expected^2/Variance)*(1-p0)/p0)
    
    crit1 <- qnorm(0.975)
    crit2 <- qnorm(1-0.025/nrow(agg))
    agg$Performance <- 1
    agg$Performance[agg$Observed < agg$Expected + crit2*sqrt(agg$Variance)] <- 2
    agg$Performance[agg$Observed < agg$Expected + crit1*sqrt(agg$Variance)] <- 3
    agg$Performance[agg$Observed < agg$Expected - crit1*sqrt(agg$Variance)] <- 4
    agg$Performance[agg$Observed < agg$Expected - crit2*sqrt(agg$Variance)] <- 5
    agg$Performance <- factor(agg$Performance, levels=1:5,
                              labels=c("Clearly worse than average", "Worse than average", "Within range",
                                       "Better than average", "Clearly better than average"))
    
    agg$OEplot <- agg$OE
    agg$OEplot[agg$OEplot>3] <- 3
    agg$toohigh <- 0
    agg$toohigh[agg$OE>3] <- 1
    agg$toohigh <- factor(agg$toohigh)
    agg$sim <- sim
    all.agg <- rbind(all.agg, agg)
  }
  attr(all.agg, "simsettings") <- simsettings
  return(all.agg)
}

expit <- function(x) exp(x)/(1+exp(x))

dosim.Logan <- function(simsettings, B=1000, progress=1) {
  # This function performs simulations according to procedures described in
  # Benchmarking survival outcomes: A funnel plot for survival data
  # It implements the pseudo-observations approach of Logan et al.
  #
  # Unpack simulation settings
  seed <- simsettings$seed
  nsim <- simsettings$nsim
  tau <- simsettings$event$tau
  
  set.seed(seed)
  all.agg <- NULL # this will contain all results
  if (progress>0) pb <- progress_bar$new(total = nsim)
  for (sim in 1:nsim) {
    if (progress>0) pb$tick()
    gd <- gendata(simsettings)
    N <- nrow(gd)
    gd$statustau <- gd$status
    gd$statustau[gd$time > tau] <- 0
    gd$timetau <- pmin(gd$time, tau)
    
    M <- max(gd$ID)
    
    # Extract censoring parameters
    first.center <- which(!duplicated(gd$ID))
    am <- gd$a[first.center]
    bm <- gd$b[first.center]

    # Observed
    sf <- survfit(Surv(time, status) ~ ID, data=gd)
    p.O <- 1 - summary(sf, times=tau, extend=TRUE)$surv
    
    # Pseudo-observations (using prodlim)
    f <- prodlim(Hist(time, status) ~ 1, data=gd)
    gd$ps <- 1 - c(as.vector(jackknife(f, times = tau))) # pseudo-values for event, not survival
    gg <- geese(ps ~ x, data=gd, id=ID, family="gaussian", corstr="ind", mean.link="logit", variance="binomial", scale.fix=TRUE)
    lp <- gg$beta[1] + gg$beta[2] * gd$x
    gd$pij <- expit(lp)
    
    gd <- as_tibble(gd) %>%
      mutate(sdpij = sqrt(pij * (1-pij)),
             rij = (ps - pij) / sdpij) %>% # Pearson residuals
      group_by(ID)
    
    pihat.boot <- matrix(NA, M, B)
    
    for (b in 1:B) {
      idx <- sample(1:N, size=N, replace=TRUE)
      gd$Yij <- gd$pij + gd$rij[idx] * gd$sdpij
      tmp <- gd %>% summarise(pihat=mean(Yij))
      pihat.boot[, b] <- tmp$pihat
    }
    
    quantiles <- apply(pihat.boot, 1, quantile, c(0.025, 0.975))
    agg <- data.frame(ID = 1:M, t(quantiles), Observed=p.O)
    agg$under <- (quantiles[1, ] > p.O)
    agg$over <- (quantiles[2, ] < p.O)
    agg$a <- am
    agg$b <- bm
    agg$sim <- sim
    all.agg <- rbind(all.agg, agg)
  }
  attr(all.agg, "simsettings") <- simsettings
  return(all.agg)
}

dosim_nonPH.Logan <- function(simsettings, B=1000, progress=1) {
  # This function performs simulations according to procedures described in
  # Benchmarking survival outcomes: A funnel plot for survival data
  # It implements the pseudo-observations approach of Logan et al.
  # This is the non-proportional hazards setting
  #
  # Unpack simulation settings
  seed <- simsettings$seed
  nsim <- simsettings$nsim
  tau <- simsettings$event$tau
  
  set.seed(seed)
  all.agg <- NULL # this will contain all results
  if (progress>0) pb <- progress_bar$new(total = nsim)
  for (sim in 1:nsim) {
    if (progress>0) pb$tick()
    gd <- gendata_nonPH(simsettings)
    N <- nrow(gd)
    gd$statustau <- gd$status
    gd$statustau[gd$time > tau] <- 0
    gd$timetau <- pmin(gd$time, tau)
    
    M <- max(gd$ID)
    
    # Observed
    sf <- survfit(Surv(time, status) ~ ID, data=gd)
    p.O <- 1 - summary(sf, times=tau, extend=TRUE)$surv
    
    # Pseudo-observations
    f <- prodlim(Hist(time, status) ~ 1, data=gd)
    gd$ps <- 1 - c(as.vector(jackknife(f, times = tau))) # pseudo-values for event, not survival
    gg <- geese(ps ~ x, data=gd, id=ID, family="gaussian", corstr="ind", mean.link="logit", variance="binomial", scale.fix=TRUE)
    lp <- gg$beta[1] + gg$beta[2] * gd$x
    gd$pij <- expit(lp)
    
    gd <- as_tibble(gd) %>%
      mutate(sdpij = sqrt(pij * (1-pij)),
             rij = (ps - pij) / sdpij) %>% # Pearson residuals
      group_by(ID)
    
    pihat.boot <- matrix(NA, M, B)
    
    for (b in 1:B) {
      idx <- sample(1:N, size=N, replace=TRUE)
      gd$Yij <- gd$pij + gd$rij[idx] * gd$sdpij
      tmp <- gd %>% summarise(pihat=mean(Yij))
      pihat.boot[, b] <- tmp$pihat
    }
    
    quantiles <- apply(pihat.boot, 1, quantile, c(0.025, 0.975))
    agg <- data.frame(ID = 1:M, t(quantiles), Observed=p.O)
    agg$under <- (quantiles[1, ] > p.O)
    agg$over <- (quantiles[2, ] < p.O)
    agg$sim <- sim
    all.agg <- rbind(all.agg, agg)
  }
  attr(all.agg, "simsettings") <- simsettings
  return(all.agg)
}

dosim_censoring <- function(simsettings, progress=1) {
  # This function only keeps track of censoring proportions,
  # as per request of a reviewer
  #
  # Unpack simsettings
  seed <- simsettings$seed
  nsim <- simsettings$nsim
  tau <- simsettings$event$tau
  
  set.seed(seed)
  all_cens <- data.frame(status0=rep(NA, nsim), status1=rep(NA, nsim))
  if (progress==1) pb <- progress_bar$new(total = nsim)
  for (sim in 1:nsim) {
    if (progress==1) pb$tick()
    gd <- gendata(simsettings)
    tbl <- table(gd$status)
    all_cens$status0[sim] <- tbl[1]
    all_cens$status1[sim] <- tbl[2]
  }
  return(all_cens)
}

dosim_both <- function(simsettings, B=1000, progress=1) {
  # This function performs simulations according to procedures described in
  # Benchmarking survival outcomes: A funnel plot for survival data
  # It implements the pseudo-observations approach of Logan et al.
  # as well as the method proposed in the paper
  #
  # Unpack simulation settings
  seed <- simsettings$seed
  nsim <- simsettings$nsim
  tau <- simsettings$event$tau
  
  set.seed(seed)
  all.agg <- NULL # this will contain all results
  if (progress==1) pb <- progress_bar$new(total = nsim)
  for (sim in 1:nsim) {
    if (progress==1) pb$tick()
    if (progress==2) deb(sim, method="cat")
    gd <- gendata(simsettings)
    N <- nrow(gd)
    gd$statustau <- gd$status
    gd$statustau[gd$time > tau] <- 0
    gd$timetau <- pmin(gd$time, tau)
    gd$pij <- NA
    
    M <- max(gd$ID)
    
    if (progress==2) cat("survfit ... ")
    cev <- coxph(Surv(timetau, statustau) ~ x, data=gd)
    sfev <- survfit(cev, newdata=gd[1, ]) # just for event time points
    tev <- sfev$time

    # Extract censoring parameters
    first.center <- which(!duplicated(gd$ID))
    am <- gd$a[first.center]
    bm <- gd$b[first.center]
    
    M <- max(gd$ID)
    ni <- table(gd$ID)
    N <- sum(ni)
    # Calculation of pij's
    if (progress==2) pb <- progress_bar$new(total = M)
    for (i in 1:M) {
      if (progress==2) pb$tick()
      whi <- which(gd$ID==i)
      gdi <- gd[whi, ]
      # Censoring distribution
      sfc <- survfit(Surv(time, status==0) ~ 1, data=gdi)
      tc <- sfc$time
      tt <- sort(unique(c(tc[tc<=tau], tev)))
      for (j in 1:ni[i]) {
        sfev <- survfit(cev, newdata=gdi[j, ])
        Sijt <- summary(sfev, times=tt, extend=TRUE)$surv
        fijt <- -diff(c(1, Sijt))
        Git <- summary(sfc, times=tt, extend=TRUE)$surv
        gd$pij[whi][j] <- sum(fijt * Git)
      }
    }
    
    #
    # Logan
    #
    
    if (progress==2) cat("calculating pseudo-observations ... ")
    
    # Observed
    sf <- survfit(Surv(time, status) ~ ID, data=gd)
    p.O <- 1 - summary(sf, times=tau, extend=TRUE)$surv
    
    # Pseudo-observations
    f <- prodlim(Hist(time, status) ~ 1, data=gd)
    gd$ps <- 1 - c(as.vector(jackknife(f, times = tau))) # pseudo-values for event, not survival
    gg <- geese(ps ~ x, data=gd, id=ID, family="gaussian", corstr="ind", mean.link="logit", variance="binomial", scale.fix=TRUE)
    lp <- gg$beta[1] + gg$beta[2] * gd$x
    gd$psij <- expit(lp)
    
    gd <- as_tibble(gd) %>%
      mutate(sdpsij = sqrt(psij * (1-psij)),
             rij = (ps - psij) / sdpsij) %>% # Pearson residuals
      group_by(ID)
    
    pihat.boot <- matrix(NA, M, B)
    
    if (progress==2) cat("bootstrapping ... ")
    
    for (b in 1:B) {
      idx <- sample(1:N, size=N, replace=TRUE)
      gd$Yij <- gd$psij + gd$rij[idx] * gd$sdpsij
      tmp <- gd %>% summarise(pihat=mean(Yij))
      pihat.boot[, b] <- tmp$pihat
    }
    
    quantiles <- apply(pihat.boot, 1, quantile, c(0.025, 0.975))

    #
    # Funnel plot
    #
    
    if (progress==2) cat("aggregating ...\n")

    p0 <- mean(gd$statustau)
    agg <- gd %>% 
      group_by(ID) %>% 
      summarise(n=n(), Observed=sum(statustau), Expected=sum(pij), Variance=sum(pij*(1-pij))) %>%
      mutate(OE=Observed/Expected) %>% 
      mutate(Z=(Observed - Expected)/sqrt(Variance)) %>%
      mutate(best_half=Z<quantile(Z,0.5)) %>%
      mutate(NNN=(Expected^2/Variance)*(1-p0)/p0)
    
    crit1 <- qnorm(0.975)
    crit2 <- qnorm(1-0.025/nrow(agg))
    agg$Performance <- 1
    agg$Performance[agg$Observed < agg$Expected + crit2*sqrt(agg$Variance)] <- 2
    agg$Performance[agg$Observed < agg$Expected + crit1*sqrt(agg$Variance)] <- 3
    agg$Performance[agg$Observed < agg$Expected - crit1*sqrt(agg$Variance)] <- 4
    agg$Performance[agg$Observed < agg$Expected - crit2*sqrt(agg$Variance)] <- 5
    agg$Performance <- factor(agg$Performance, levels=1:5,
                              labels=c("Clearly worse than average", "Worse than average", "Within range",
                                       "Better than average", "Clearly better than average"))
    
    # Aggregating Logan results and adding to agg
    agg <- cbind(agg, data.frame(t(quantiles), KM=p.O))
    agg$under <- (quantiles[1, ] > p.O)
    agg$over <- (quantiles[2, ] < p.O)
    
    agg$a <- am
    agg$b <- bm
    agg$sim <- sim
    all.agg <- rbind(all.agg, agg)
  }
  attr(all.agg, "simsettings") <- simsettings
  return(all.agg)
}
