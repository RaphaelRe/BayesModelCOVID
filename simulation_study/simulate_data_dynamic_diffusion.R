###############################################################################
# The same skript as 'simulate_data_dynamic.R, but with 4 strata and possible 
# between the strata
###############################################################################

PATH_DATA = "../data/"
N_CORES <- 5 # the number of used cores to generate data. Simulation time is not 
             # THAT long but can save a view minutes if N_DATASETS is high

### libraries
library(data.table)
library(magrittr)
library(Rcpp)
library(lubridate)
library(pbapply)
library(ggplot2)

### overall seed
set.seed(123)

# number of generated datasets
N_DATASETS <- 100

### other parameters
N_COUNTRY <- 10
N_DAYS <- 600
n_strata <- 4 # number of different age groups in population
POPULATION <- 100000000 # arbitrary countries, can be the same over all countries - in each stratum
POPULATION <- rep(POPULATION/n_strata, n_strata)


### load time distributions
gamma <- fread(paste0(PATH_DATA, "gamma_generation_time.csv"))[[3]]
Xi_C <- fread(paste0(PATH_DATA, "Xi_C_incubation_period.csv"))[[3]]
Xi_repC <- fread(paste0(PATH_DATA, "Xi_R_reporting_times_weekdays_estimated_lgl.csv"))[,-c(1:2)] 
Xi_D <- fread(paste0(PATH_DATA, "Xi_D_symptoms_to_death_weekdays.csv"))[[3]]
Xi_H <- fread(paste0(PATH_DATA, "XiH_all.csv"))[[1]] 
Xi_Hicu <- fread(paste0(PATH_DATA, "XiH_all.csv"))[[3]]


## adapt reporting patterns for a part of the countries (assume that these countries use a reporting patter where they do not report anything on this specific day)
# reported cases
Xi_repC_2 <- lapply(1:7, function(day){
  xi_tmp <- Xi_repC[[day]]
  xi_tmp[seq(8-day, length(xi_tmp), by = 7)] <- 0
  xi_tmp <- xi_tmp/sum(xi_tmp)
}) %>% as.data.table %>% 
  set_names(names(Xi_repC))


### define parameters for tau, each country has its own tau deviating from the overall
mean_tau <- 10
sd_tau <- 2


### overdispersion negative binomial distributions, 1000 is rather small and gives stability for the simulation
sapply(c("inf", "dea", "rep", "hos", "icu"), function(x){
  name <- paste0("phi_", x)
  assign(name, 1000, envir = .GlobalEnv)
})


### define R0 and sd for sampling a country specific R0 which deviates by sd_R0
mean_R0 <- 3.25
sd_R0 <- 0.1


####################
###### Effects on Rt

### propability reinfection
p_reinf <- 0.16
# p_reinf <- 1 # would imply a correction factor1 of 0


### Variants of concern
# over-contagiousness of VOCs
beta_alpha <- 1.4
beta_delta <- 2.4
# beta_alpha <- 1 # would imply same contagiousness as wild type
# beta_delta <- 1

# load fractions from approximated data, use UK, could anything, but Uk has 
# a little bit more prevalence of delta
vocs <- fread(paste0(PATH_DATA, "variants_of_concern.csv"))[country == "UnitedKingdom" ,
                                          .(date, proportion_original, 
                                            proportion_alpha, proportion_delta)][1:N_DAYS]
# melt(vocs, id.vars="date") %>%
#   ggplot(aes(date, value, color=variable))+geom_line()

# define fractions
proportion_original <- vocs$proportion_original
proportion_alpha <- vocs$proportion_alpha
proportion_delta <- vocs$proportion_delta


### vaccinations
# vaccinations can also be defined arbitrary; defined in % of population
# start at day 280, second vaccination is 40 days later
vaccinations1 <- rep(0, N_DAYS)
vaccinations1[281:N_DAYS] <- seq(0,0.8, len=N_DAYS-280-40) %>% append(rep(0.8,40))

vaccinations2 <- c(rep(0,40),vaccinations1[1:(N_DAYS-40)])

# alternative where vaccinations2 converges to 80 at time 500, version above probably better
# vaccinations2 <- rep(0, N_DAYS)
# vaccinations2[(281+15):N_DAYS] <- seq(0,0.8, len=N_DAYS-280-15)
# vaccinations2[(281+40):N_DAYS] <- seq(0,0.8, len=N_DAYS-280-40)

# plot(1:N_DAYS, vaccinations1)
# points(1:N_DAYS,vaccinations2)

# effects of vaccinations on R0
vacc_on_R1 <- 0.6
vacc_on_R2 <- 0.3

# effects of vaccinations on ifr
vacc_on_ifr <- 0.8 # 80% reduced ifr


### define piD
################################## old versio n##################################

# note: every country can use the same piD since it is deterministic
piD <- rep(0.007, N_DAYS)
# adapt it after 280+14, i.e. reduce it, to simulate vaccination
piD[(281+14):N_DAYS] <- piD[(281+14):N_DAYS] * c(seq(1 ,0.1, len=N_DAYS-280-40-14), rep(0.1, 40))

# add effect of variats of concern:
# assume effect of alpha to be 50% and delta 100% more severe
piD <- piD*proportion_original + piD*proportion_alpha*1.5 + piD*proportion_delta*2
# plot(piD)

# individual pi_D for each strata
piD <- replicate(n_strata , piD)

################################## ##################################

################################## with different IFRs in each strata


piD <- rep(0.003, N_DAYS) %>% replicate(n_strata, .)
piD <- sapply(1:4, function(i) piD[, i]* i)
rowMeans(piD)

piD[(281+14):N_DAYS, ] <- piD[(281+14):N_DAYS, ] * c(seq(1 ,0.1, len=N_DAYS-280-40-14), rep(0.1, 40))
piD <- piD*proportion_original + piD*proportion_alpha*1.5 + piD*proportion_delta*2
# matplot(piD)

# mean over all strata. as all strata have the same proportion in the popultation 
# it is sufficient to take the simple average here
piD_overall_mean <- rowMeans(piD)
# lines(1:length(piD_overall_mean), piD_overall_mean)



### define piH & piHicu
# basic piH is 0.9
# basic piHicu is 0.2
# take ratio from piD to get actual piHs
piH <- 0.9
piHicu <- 0.2

piH <- apply(piD, 2, function(x){
  sapply(x, function(y) y/x[1])*piH
})
piHicu <- apply(piD, 2, function(x){
  sapply(x, function(y) y/x[1])*piHicu
})
# piH <- sapply(piD, function(x) x/piD[1])*piH
# piHicu <- sapply(piD, function(x) x/piD[1])*piHicu


################################## ##################################
################################## with different IFRs in each strata

piH <- rep(0.9, n_strata)
piHicu <- rep(0.2, n_strata)

piH <- rep(0.35, n_strata) * 1:4
piHicu <- rep(0.08, n_strata) * 1:4

scale_factor <- apply(piD, 2, function(x){
  sapply(x, function(y) y/x[1])
})

piH <- sweep(scale_factor, MARGIN=2, piH, `*`)
piHicu <- sweep(scale_factor, MARGIN=2, piHicu, `*`)

piH_overall_mean <- rowMeans(piH)
piHicu_overall_mean <- rowMeans(piHicu)
# matplot(piH)
# lines(1:length(piH_overall_mean), piH_overall_mean)



### rhos aka reporting ratio
### make only 5 periods. This helps with the sampling time in the algorithm, but does not affect results
rhos <- 1/5:1



cppFunction(
  ## calculates the partial convolution over all time points (for all strata)
  'NumericMatrix partial_conv_cpp(NumericMatrix x, 
                                NumericVector dist,
                                NumericMatrix p,
                                int phi) {
    int n = x.nrow();
    int K = x.ncol();
    NumericMatrix res(n, K);
    for (int k = 0; k < K; ++k) {
      for (int t = 0; t < n; ++t) {
        double sumut = 0.0;
        for(int u = 0; u < t; ++u) {
          sumut += (x(u, k) * dist[t-u]);
        }
        //Rcout << sumut << " ";
        int val = rnbinom_mu(1, phi, sumut*p(t, k))[0];
        res(t, k) = val;
      }
    }
    //Rcout << res[109, 3] << " ";
    return res;
  }')



cppFunction(
  # calculates the partial convolution for one time point over all strata
  'NumericVector calc_one_t_cpp(int t, NumericMatrix x, 
                                NumericVector dist,
                                NumericVector p,
                                int phi) {
    int K = x.ncol();
    
    NumericVector res(K);
    t -= 1;
    for (int k = 0; k < K; ++k) {
      double sumut = 0.0;
      for(int u = 0; u < t; ++u) {
        sumut += (x[u] * dist[t-u]);
      }
      res(k) = rnbinom_mu(1, phi, sumut*p[k])[0];
    }
    return res;
  }')



cppFunction(
  # calculates the partial convolution for one time point over all strata, with diffusion between them
  'NumericVector calc_renewal_diffusion_cpp(int t, NumericMatrix x, 
                                         NumericVector dist,
                                         NumericVector R,
                                         int phi, 
                                         NumericVector diffusion) {
  int K = x.ncol();
  
  t -= 1;
  NumericVector mu_k(K);
  NumericVector sumut(K);
  // Rcout << sumut << " ____";
  
  for(int u = 0; u < t; ++u) {
    for (int k = 0; k < K; ++k) {
      for (int ktmp = 0; ktmp < K; ++ktmp) {
        if (ktmp == k) {
          sumut[k] += (x(u, k) * dist[t-u]) * diffusion[0];
        } else {
          sumut[k] += (x(u, ktmp) * dist[t-u]) * diffusion[1];
        }
      }
    }
    mu_k = sumut*R;
  }
  // sample a negative binomial random variable for each mu_k
  NumericVector res(K);
  for (int k = 0; k < K; ++k) {
    res(k) = rnbinom_mu(1, phi, mu_k(k))[0];
  }
  return res;
}')

# short code to test
# xxx <- cbind(c(1,1,1), c(2,2,2), c(100,100,100), c(10,10,10))
# ppp <- rep(1, 4)
# calc_renewal_diffusion_cpp(3, xxx, rep(1,5), ppp, 10000, c(1,0))  # standard without diffusion to other groups
# calc_renewal_diffusion_cpp(3, xxx, rep(1,5), ppp, 10000, c(4,0.2))  # with diffusion



## define sigmoid function which defines the probability if NPI k gets active
sigmoid <- function(x, b0,b1){
  1/(1+exp(-(b0+b1*x)))
}

### define NPI effects
# define 5 mean effects NPIs which will be used for all countries
NPI_means <- c(0.22, 0.25, 0.3, 0.4, 0.45) 

# assume the same sd for all NPIs here
NPI_sd <- 0.01
# function to generate 5 effects from means
generate_NPIs <- function(){
  rnorm(length(NPI_means),mean=NPI_means, sd=NPI_sd) %>% 
    set_names(paste0("NPI_", 1:length(NPI_means)))
}

# check mean influence on Rt
3.25*exp(-sum(NPI_means[-5])) #  value is slightly > 1 implying a Rt > 1 but 
                              # correction factors will push this below 1


# this bbs define the probability which is used to set a NPI active
bb0 <- -5 # intercept
# for 5 npis
bb1 <- 0.001 # NPI1
bb2 <- 0.0005 # NPI2
bb3 <- 0.0004 # NPI3
bb4 <- 0.0003 # NPI4
bb5 <- 0.0001 # NPI5

# viz sigmoids to get an idea of the probs
# grid <- 0:20000
# plot(grid, sigmoid(grid, bb0, bb1), type = "l", ylim=c(0,1))
# lines(grid, sigmoid(grid, bb0, bb2), type = "l", col=2)
# lines(grid, sigmoid(grid, bb0, bb3), type = "l", col=3)
# lines(grid, sigmoid(grid, bb0, bb4), type = "l", col=4)
# lines(grid, sigmoid(grid, bb0, bb5), type = "l", col=5)



### Define function to simulate data
#  This function simulates the data for ONE country using the objects from above
simulate_data_country <- function(seed=321, diffusion=F){
  
  if (diffusion) {
    diffusion <- 0.2
    diffusion <- c(1-diffusion * (n_strata-1), diffusion)  # used to calculate the diffusion to other strata, first is dynamics in own strata, second is diffusion from other strata
  } else {
    diffusion <- c(1, 0)
  }
  
  set.seed(seed)
  
  # sample the used reporting scheme for the country
  reporting_scheme <- sample(c(T,F), 1)
  
  if (reporting_scheme) {
    Xi_repC_m <- Xi_repC
  } else {
    Xi_repC_m <- Xi_repC_2
  }
  
  ### sample NPI effects
  NPIs <- generate_NPIs()
  
  data <- data.table(
    date = seq.Date(as.Date("2020-01-01"), by="day", len=N_DAYS)
    )
  
  # add true NPI effects to know the truth
  data[, true_val_NPI1 := NPIs[1]]
  data[, true_val_NPI2 := NPIs[2]]
  data[, true_val_NPI3 := NPIs[3]]
  data[, true_val_NPI4 := NPIs[4]]
  data[, true_val_NPI5 := NPIs[5]]
  
  # sample tau for this country given tau_mean
  tau <- rnorm(1, mean_tau, sd_tau)
  
  # sample an R0 for this country given mean_R0 
  R0 <- rnorm(1, mean_R0, sd_R0)
  
  # add effect of VOCs
  R0_t <- R0*proportion_original + 
    proportion_alpha*R0*beta_alpha + 
    proportion_delta*R0*beta_delta

  ### sample 1st value, i.e. index cases from poisson with rate tau
  I_1 <- rnbinom(n_strata, mu=tau, size=phi_inf)

  I_t <- matrix(-1, N_DAYS, n_strata)
  R_t_vec <- matrix(-1, N_DAYS, n_strata)
  c1_t_vec <- matrix(-1, N_DAYS, n_strata)
  c2_t_vec <- matrix(-1, N_DAYS, n_strata)
  cfull_t_vec <- matrix(-1, N_DAYS, n_strata)
  C_t <- matrix(-1, N_DAYS, n_strata)
  Hicu_t <- matrix(-1, N_DAYS, n_strata)
  
  I_t[1, ] <- I_1
  R_t_vec[1,] <- R0
  Hicu_t[1,] <- 0
  c1_t_vec[1,] <- 0
  c2_t_vec[1,] <- 0
  cfull_t_vec[1,] <- 1
  
  #### generate infections dynamically
  # at each time point t, calculate a probability whether NPI_k gets online
  
  # assign states of NPIs
  state_NPI1 <- rep(F, N_DAYS)
  state_NPI2 <- rep(F, N_DAYS)
  state_NPI3 <- rep(F, N_DAYS)
  state_NPI4 <- rep(F, N_DAYS)
  state_NPI5 <- rep(F, N_DAYS)
  
  # NPI1 gets active after 30 days for the rest of the pandemic
  # state_NPI1[31:N_DAYS] <- T
  
  # calculate the first Ct, loop starts at t=2, need decimals and Ct[1]
  Sumut <-  I_t[1, ] * Xi_C[1]
  decimals <- Sumut-floor(Sumut) # 4x1
  
  # if(decimals<0.5){
  #   decimals_last_day <- decimals
  # } else{
  #   decimals_last_day <- decimals-1
  # }
  
  decimals_last_day <- ifelse(decimals<0.5, decimals, decimals-1)
  C_t[1, ] <- round(Sumut)
  
  for (t in 2:N_DAYS){
    # cat(t)
    # cat("-")
    
    ### check whether NPI_k gets active at t
    # if NPI gets active it stays active for 60 to 120 days
    for (npi in 1:5) {
      # paste("npi:", npi) %>% print
      NPI <- paste0("state_NPI",npi)
      if (!get(NPI)[t]) {
        # p_set_active <- sigmoid(Hicu_t[t-1], b0=bb0, b1=bb1)
        p_set_active <- sigmoid(sum(Hicu_t[t-1,]), b0=bb0, b1=get(paste0("bb",npi)))
        # print(Hicu_t[t-1,])
        # print(p_set_active)
        if (as.logical(rbinom(1, 1, p_set_active))){
          active_duration <- sample(60:120,1)
          assign(NPI,
                 replace(get(NPI),
                         t:(t+active_duration), T)[1:N_DAYS]
                 )
          # state_NPI2[t:(t+active_duration)] <- T # original version, only here as alterative to assign each NPI seperately without the outer loop over npis
        }
      }
    }
    
    # define NPIs as matrix, I know that is not efficient, since it is done at every t
    # currently too lazy to solve this
    # possible solution: define it before the loop and assign it in each column and row
    NPI_path <- cbind(state_NPI1, state_NPI2, state_NPI3, state_NPI4, state_NPI5)
    # calculate R without correction factors
    R_t <- R0_t[1:t] * exp(-(NPI_path %*% NPIs)[1:t]) # probably the last value would be enough in the calculate...then the following line is not necessary
    Rt <- R_t[t]
    
    # correction factors c1 and c2
    # c1_t <- sum(I_t[1:(t-1)] / POPULATION) * (1-p_reinf)
    # special case for t = 2 as R evaluates the matrix to a vector
    if (t == 2){
      c1_t <- (I_t[1:(t-1), ] / POPULATION) * (1-p_reinf)
    } else {
      c1_t <- colSums(I_t[1:(t-1), ] / POPULATION) * (1-p_reinf)  
    }
    c2_t <- (vaccinations1[t]) * vacc_on_R1 +
            (vaccinations2[t]) * vacc_on_R2
    
    cf_full <- (1 - c1_t - c2_t * (1-c1_t))
    
    Rt <- Rt * cf_full  # Rt is overloaded to a vector now
    
    c1_t_vec[t, ] <- c1_t
    c2_t_vec[t, ] <- c2_t
    
    cfull_t_vec[t, ] <- cf_full
    R_t_vec[t, ] <- Rt # just here to keep track and monitor it
    
    # given final Rt, finally calculate I_t and Hicu_t
    
    ### calc I_t
    # slow version
    # sumut <- 0
    # for (u in 1:(t-1)) {
    #   sumut <- sumut + I_t[u] * gamma[t-u+1]
    # }
    # I_t[t] <- rnbinom(1, mu = sumut*Rt, size = phi_inf)
    # I_t[t,] <- calc_one_t_cpp(t, I_t[1:t, ], gamma, Rt, phi_inf)
    if (any(Rt < 0)) {
      browser()
    }
    I_t[t,] <- calc_renewal_diffusion_cpp(t, I_t[1:t, ], gamma, Rt, phi_inf, diffusion)
    
    if (any(is.na(I_t[t,]))) {
      browser()
    }
    # print("test whether the comment version gives exctly the same results!!!, Maybe add a t-1 in the cpp call")
    
    #### generate Ct
  
    Sumut <- decimals_last_day
    for(u in 1:t){
      Sumut <- Sumut + I_t[u,] * Xi_C[t-u+1]
    }
    decimals <- Sumut-floor(Sumut)
    decimals_last_day <- ifelse(decimals<0.5, decimals, decimals-1)
    # if(decimals<0.5){
    #   decimals_last_day <- decimals
    # } else{
    #   decimals_last_day <- decimals-1
    # }
    C_t[t,] <- round(Sumut)
    
    ### calc Hicu_t
    #### generate Hicut
    # slow version
    # Sumut <- 0
    # for(u in 1:t){
    #   Sumut <- Sumut + C_t[u] * Xi_Hicu[t-u+1]
    # }
    # Hicu_t[t] <- rnbinom(1, mu = Sumut*piHicu[t], size = phi_icu)
    Hicu_t[t, ] <- calc_one_t_cpp(t, C_t, Xi_Hicu, piHicu[t, ], phi_icu) # Ct is longer, but only the first t values should be used, nice sanity check
  }
  
  data[, NPI1 := as.numeric(state_NPI1)]
  data[, NPI2 := as.numeric(state_NPI2)]
  data[, NPI3 := as.numeric(state_NPI3)]
  data[, NPI4 := as.numeric(state_NPI4)]
  data[, NPI5 := as.numeric(state_NPI5)]
  
  
  # assign the values to the data frame
  # also keep track of the sum over all strata for each t
  data[, paste0("I_t_", 1:n_strata) := lapply(seq_len(ncol(I_t)), function(i) I_t[,i])]
  data[, I_t := rowSums(I_t)]
  
  data[, paste0("C_t_", 1:n_strata) := lapply(seq_len(ncol(C_t)), function(i) C_t[,i])]
  data[, C_t := rowSums(C_t)]
  
  data[, paste0("Hicu_t_", 1:n_strata) := lapply(seq_len(ncol(Hicu_t)), function(i) Hicu_t[,i])]
  data[, Hicu_t := rowSums(Hicu_t)]
  
  
  # do the same for other params
  data[, paste0("R_t_", 1:n_strata) := lapply(seq_len(ncol(R_t_vec)), function(i) R_t_vec[,i])]
  data[, R_t := rowSums(R_t_vec)]
  
  data[, paste0("cf1_", 1:n_strata) := lapply(seq_len(ncol(c1_t_vec)), function(i) c1_t_vec[,i])]
  data[, cf1 := rowSums(c1_t_vec)]
  
  data[, cf2 := c2_t_vec[,1]]
  
  data[, paste0("cf_full_", 1:n_strata) := lapply(seq_len(ncol(cfull_t_vec)), function(i) cfull_t_vec[,i])]
  data[, cf_full := rowSums(cfull_t_vec)]
  
  # data[, cf_full := cfull_t_vec]
  # data[, I_t := I_t]
  # data[, C_t := C_t]
  # data[, Hicu_t := Hicu_t]
  # 
  # data[, R_t := R_t_vec]
  # data[, cf1 := c1_t_vec]
  # data[, cf2 := c2_t_vec]
  # data[, cf_full := cfull_t_vec]
  # 
  # just to check whether the simulation is doing some reasonable stuff
  # layout(1:2)
  # I_t %>% plot
  # R_t_vec %>% plot
  
  
  
  #### generate repCt
  
  # we need the weekdays since the model must know which column to pick from Xi_repC
  data[, weekday := lubridate::wday(date, week_start = T)]
  # also need identifier for rho period
  if (N_DAYS %% 10 != 0) {
    stop("N_DAYS should be a multiple of 10")
  }
  rho_periods <- rep(1:10, each=N_DAYS/10)
  rho_periods <- rep(1:5, each=N_DAYS/5)
  data[, rho_period := rho_periods]
  
  # calc reporting delay
  ErepC_t <- rep(-1, N_DAYS)
  repC_t <- rep(-1, N_DAYS)
  for(t in 1:N_DAYS){
    # need a weekday for each day
    Sumut <- 0
    for(u in 1:t){
      wd <- data$weekday[u]
      # reporting can be don eon all strata
      Sumut <- Sumut + data$C_t[u] * Xi_repC_m[[wd]][t-u+1]
    }
    
    rho_tmp <- rhos[rho_periods[t]]
    ErepC_t[t] <- Sumut*rho_tmp
    repC_t[t] <- rnbinom(1, mu = ErepC_t[t], size = phi_rep)
  }
  data[, repC_t := repC_t]
  data[, ErepC_t := ErepC_t]
 
  #### generate Dt
  # D_t <- rep(-1, N_DAYS)
  # 
  # for(t in 1:N_DAYS){
  #   Sumut <- 0
  #   for(u in 1:t){
  #     Sumut <- Sumut + C_t[u] * Xi_D[t-u+1]
  #   }
  #   D_t[t] <- rnbinom(1, mu = Sumut*piD[t], size = phi_dea)
  # }
  
  D_t <- partial_conv_cpp(C_t, Xi_D, piD, phi_dea)
  
  data[, paste0("D_t_", 1:n_strata) := lapply(seq_len(ncol(D_t)), function(i) D_t[,i])]
  data[, D_t := rowSums(D_t)]
  # data[, D_t := D_t]
  
    #### generate Ht
  # H_t <- rep(-1, N_DAYS)
  # 
  # for(t in 1:N_DAYS){
  #   Sumut <- 0
  #   for(u in 1:t){
  #     Sumut <- Sumut + C_t[u] * Xi_H[t-u+1]
  #   }
  #   H_t[t] <- rnbinom(1, mu = Sumut*piH[t], size = phi_hos)
  # }
  
  
  H_t <- partial_conv_cpp(C_t, Xi_H, piH, phi_hos)
  
  data[, paste0("H_t_", 1:n_strata) := lapply(seq_len(ncol(H_t)), function(i) H_t[,i])]
  data[, H_t := rowSums(H_t)]
  # data[, H_t := H_t]
  
  
  
  ### add seasonality
  # NOT REQUIRED IN SIMULATION
  # define_dummy <- function(x, months){
  #   return(as.numeric(x %in% months))
  # }
  # 
  # define_variable <- function(d, name, months){
  #   d[,(name) := define_dummy(month(date), months)]
  #   return(d)
  # }
  # 
  # define_variable(data, "winter", c(12,1,2))
  # define_variable(data, "spring", c(3,4,5))
  # # define_variable(data, "summer", c(6,7,8)) # actually not required
  # define_variable(data, "autumn", c(9,10,11))
  # 
  
  # add further required colmuns
  # data[, ifr_t_m := piD]
  data[, paste0("ifr_t_m_", 1:n_strata) := lapply(seq_len(ncol(piD)), function(i) piD[,i])]
  data[, ifr_t_m := piD_overall_mean]
  data[, first_vaccination := vaccinations1]
  data[, second_vaccination := vaccinations2]
  data[, alpha := proportion_alpha]
  data[, delta := proportion_delta]
  
  return(data[])
}



# test it for one country
data <- simulate_data_country(2463, diffusion = T)
data[NPI1 == 0, NPI1 := NA]
data[NPI2 == 0, NPI2 := NA]
data[NPI3 == 0, NPI3 := NA]
data[NPI4 == 0, NPI4 := NA]
data[NPI5 == 0, NPI5 := NA]
data[NPI1 == 1, NPI1 := NPI1]
data[NPI2 == 1, NPI2 := NPI2*1.01]
data[NPI3 == 1, NPI3 := NPI3*1.02]
data[NPI4 == 1, NPI4 := NPI4*1.03]
data[NPI5 == 1, NPI5 := NPI5*1.04]
g1 <- ggplot(data)+
  geom_line(aes(date, I_t_1), col = "grey20")+
  geom_line(aes(date, I_t_2), col = "grey40")+
  geom_line(aes(date, I_t_3), col = "grey60")+
  geom_line(aes(date, I_t_4), col = "grey80")+ 
  geom_line(aes(date, I_t), col = "grey0")
  # geom_line(aes(date, C_t), col="darkgrey")+
  # geom_line(aes(date, repC_t), col="darkorange")+
  # geom_line(aes(date, D_t), col="darkred")+
  # geom_line(aes(date, Hicu_t), col="darkblue")+
  # geom_line(aes(date, H_t), col="blue")
g2 <- ggplot(data)+
  geom_line(aes(date, NPI1),col=2)+
  geom_line(aes(date, NPI2),col=3)+
  geom_line(aes(date, NPI3),col=4)+
  geom_line(aes(date, NPI4),col=5)+
  geom_line(aes(date, NPI5),col=6)
egg::ggarrange(g1,g2)



### use simulate_data_country function to generate N_COUNTRY countries
generate_data <- function(seed=321, write=F, path=paste0(PATH_DATA, "data_simulated.csv"), diffusion=F){
  set.seed(seed)
  seeds <- sample(1:9999, N_COUNTRY)
  dds <- lapply(seeds, simulate_data_country, diffusion = diffusion) %>% 
    set_names(paste0("country", LETTERS[1:N_COUNTRY])) %>% 
    rbindlist(idcol = "country")
  if (write){
    fwrite(dds, path)
  }
  return(dds)
}


# test for many countries
d_tmp <- generate_data(diffusion = T)
ggplot(d_tmp[country == "countryA"])+
  geom_line(aes(date, I_t))+
  geom_line(aes(date, I_t_1), col = "green1")+
  geom_line(aes(date, I_t_2), col = "green2")+
  geom_line(aes(date, I_t_3), col = "green3")+
  geom_line(aes(date, I_t_4), col = "green4")+ 
  geom_line(aes(date, C_t), col="darkgrey")+
  geom_line(aes(date, repC_t), col="darkorange")+
  geom_line(aes(date, D_t_1*20), col="darkred")+
  geom_line(aes(date, D_t_2*20), col="darkred")+
  geom_line(aes(date, D_t_3*20), col="darkred")+
  geom_line(aes(date, D_t_4*20), col="darkred")+
  geom_line(aes(date, D_t*20), col="darkred", linetype = "dashed")+
  geom_line(aes(date, Hicu_t), col="darkblue")+
  geom_line(aes(date, H_t), col="blue")+
  facet_wrap(~country, scales = "free")



### generate N_DATASETS datasets with diffusion and without diffusion
tt <- Sys.time()
N_DATASETS <- 20
set.seed(123)
seeds <- sample(1:9999, N_DATASETS)

datasets <- parallel::mcmapply(function(ds, seed){
                              generate_data(seed, write=T, path=paste0(PATH_DATA, "/simulated_data/diffusion/",ds), diffusion = T) 
                            }, 
                            ds=paste0("data_sim_5NPIs_", 1:N_DATASETS, "_diffusion_varIFR.csv"),
                            seed=seeds,
                            mc.cores = N_CORES, SIMPLIFY = F)

Sys.time()-tt



plot_dsets <- function(ds){
  ggplot(ds)+
    geom_line(aes(date, I_t))+
    geom_line(aes(date, C_t), col="darkgrey")+
    geom_line(aes(date, repC_t), col="darkorange")+
    geom_line(aes(date, D_t), col="darkred")+
    geom_line(aes(date, Hicu_t), col="darkblue")+
    geom_line(aes(date, H_t), col="blue")+
    facet_wrap(~country, scales = "free")
}


egg::ggarrange(plots = lapply(datasets[1:4], plot_dsets))





