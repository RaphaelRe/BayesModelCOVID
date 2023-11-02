###############################################################################
# This script simulates data according to the model
###############################################################################

PATH_DATA = "../data/"

N_CORES <- 50 # the number of used cores to generate data. Simulation time is not 
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
POPULATION <- 100000000 # arbitrary population size, can be the same over all countries


### load time distributions
gamma <- fread(paste0(PATH_DATA, "gamma_generation_time.csv"))[[3]]
Xi_C <- fread(paste0(PATH_DATA, "Xi_C_incubation_period.csv"))[[3]]
Xi_repC <- fread(paste0(PATH_DATA, "Xi_R_reporting_times_weekdays_estimated_lgl.csv"))[,-c(1:2)] 
Xi_D <- fread(paste0(PATH_DATA, "Xi_D_symptoms_to_death_weekdays.csv"))[[3]]
Xi_H <- fread(paste0(PATH_DATA, "XiH_all.csv"))[[1]]
Xi_Hicu <- fread(paste0(PATH_DATA, "XiH_all.csv"))[[3]]


### define parameters for tau, each country has its own tau deviating from the overall
mean_tau <- 10
sd_tau <- 2


### overdispersion negative binomial distributions, 1000 seens to give stability in the simulation
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
# p_reinf <- 1 # would imply a correction factor1 of 0, i.e. full reinfection


### Variants of concern
# over-contagiousness of VOCs
beta_alpha <- 1.4
beta_delta <- 2.4
# beta_alpha <- 1 # would imply same contagiousness as original type
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

# plot(1:N_DAYS, vaccinations1)
# points(1:N_DAYS,vaccinations2)

# effects of vaccinations on R0
vacc_on_R1 <- 0.6
vacc_on_R2 <- 0.3

# effects of vaccinations on ifr
vacc_on_ifr <- 0.8 # 80% reduced ifr


### define piD
# note: every country can use the same piD since it is deterministic
piD <- rep(0.007, N_DAYS)
# adapt it after 280+14, i.e. reduce it, to simulate vaccination
piD[(281+14):N_DAYS] <- piD[(281+14):N_DAYS] * c(seq(1 ,0.1, len=N_DAYS-280-40-14), rep(0.1, 40))

# add effect of variats of concern:
# assume effect of alpha to be 50% and delta 100% more severe
piD <- piD*proportion_original + piD*proportion_alpha*1.5 + piD*proportion_delta*2
# plot(piD)

### define piH & piHicu
# basic piH is 0.9
# basic piHicu is 0.2
# take ratio from piD to get actual piHs
piH <- 0.9
piHicu <- 0.2

piH <- sapply(piD, function(x) x/piD[1])*piH
piHicu <- sapply(piD, function(x) x/piD[1])*piHicu



### rhos aka reporting ratio
### make only 5 periods. This helps with the sampling time in the algorithm as there a less parameters to estimate
rhos <- 1/5:1



###############################################################################
# for the simulation function which generates one dataset we need some helpers

# short cpp function to make a faster convolution and sampling, saves a few seconds
# can be used for calculation of Dt, Ht, Hicut
# NOT for repCt, Ct and It, need some special tweaks, maybe I'll adapt that later
# old convolution is still as comment in the code
cppFunction(
  'NumericVector partial_conv_cpp(NumericVector x, 
                                NumericVector dist,
                                NumericVector p,
                                int phi) {
    int n = x.size();
    NumericVector res(n);
    for (int t = 0; t < n; ++t) {
      double sumut = 0.0;
      for(int u = 0; u < t; ++u) {
        sumut += (x[u] * dist[t-u]);
      }
      res[t] = rnbinom_mu(1, phi, sumut*p[t])[0];
      //res[t] = sumut;
    }
    return res;
  }')


#### test if the function works
# C_t <- rep(1:50, N_DAYS/50)
# piD <- rep(1, N_DAYS)
# phi_dea <- 10000000
# 
# D_t <- rep(-1, N_DAYS)
# for(t in 1:N_DAYS){
#  Sumut <- 0
#  for(u in 1:t){
#    Sumut <- Sumut + C_t[u] * Xi_D[t-u+1]
#  }
#  # D_t[t] <- rnbinom(1, mu = Sumut*piD[t], size = phi_dea)
#  D_t[t] <- Sumut
# }
# plot(D_t)
# partial_conv_cpp(C_t, Xi_D, piD, phi_dea) %>% points(1:500, ., col=2)




cppFunction(
  'int calc_one_t_cpp(int t, NumericVector x, 
                                NumericVector dist,
                                double p,
                                int phi) {
    int n = x.size();
    int res;
    t -= 1;
    double sumut = 0.0;
    for(int u = 0; u < t; ++u) {
      sumut += (x[u] * dist[t-u]);
    }
    res = rnbinom_mu(1, phi, sumut*p)[0];
    return res;
  }')

# r version
# foo <- function(t, x, dist, p, phi){
#   sumut <- 0
#   for (u in 1:(t-1)) {
#     sumut <- sumut + x[u] * dist[t-u+1]
#   }
#   rnbinom(1, mu = sumut*p, size = phi)
# }

#  compare whether cpp function is doning the right thing
# calc_one_t_cpp(100, 1000:1100, gamma, 2, 1000)
# foo(100, 1000:1100, gamma, 2, 1000)
# 
# system.time(replicate(10000, calc_one_t_cpp(100, 1000:1100, gamma, 2, 1000)))
# system.time(replicate(10000, foo(100, 1000:1100, gamma, 2, 1000)))
# 
# hist(replicate(100000, calc_one_t_cpp(100, 1000:1100, gamma, 2, 1000)))
# hist(replicate(100000, foo(100, 1000:1100, gamma, 2, 1000)), add = T, col=2)
# 



## define sigmoid function which defines the probability whether NPI k gets active
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

# this plot shows th eprobability that the k-th NPI is getting active as afunction of the ICU occupancy
grid <- 0:20000
plot(grid, sigmoid(grid, bb0, bb1), type = "l", ylim=c(0,1))
lines(grid, sigmoid(grid, bb0, bb2), type = "l", col=2)
lines(grid, sigmoid(grid, bb0, bb3), type = "l", col=3)
lines(grid, sigmoid(grid, bb0, bb4), type = "l", col=4)
lines(grid, sigmoid(grid, bb0, bb5), type = "l", col=5)



### define function to simulate data
simulate_data_country <- function(seed=321){
  # This function simulates the data for ONE country using the objects from above
  set.seed(seed)
  
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
  
  ### sample 1st value, i.e. index cases from nbinom with rate tau
  I_1 <- rnbinom(1, mu=tau, size=phi_inf)
  
  ### use generation time to generate the dynamics of the pandemic
  # calc container initial conditions
  I_t <- rep(-1, N_DAYS)
  R_t_vec <- rep(-1, N_DAYS) # just here to keep track and monitor it
  c1_t_vec <- rep(-1, N_DAYS) # just here to keep track and monitor it
  c2_t_vec <- rep(-1, N_DAYS) # just here to keep track and monitor it
  cfull_t_vec <- rep(-1, N_DAYS) # just here to keep track and monitor it
  C_t <- rep(-1, N_DAYS)
  Hicu_t <- rep(-1, N_DAYS)
  
  I_t[1] <- I_1
  R_t_vec[1] <- R0
  Hicu_t[1] <- 0
  c1_t_vec[1] <- 0
  c2_t_vec[1] <- 0
  cfull_t_vec[1] <- 1
  
  #### generate infections dynamically
  # at each time point t, calculate a probability whether NPI_k gets active
  
  # assign states of NPIs
  state_NPI1 <- rep(F, N_DAYS)
  state_NPI2 <- rep(F, N_DAYS)
  state_NPI3 <- rep(F, N_DAYS)
  state_NPI4 <- rep(F, N_DAYS)
  state_NPI5 <- rep(F, N_DAYS)
  
  # calculate the first Ct, loop starts at t=2, need decimals and Ct[1]
  Sumut <-  I_t[1] * Xi_C[1]
  decimals <- Sumut-floor(Sumut)
  if(decimals<0.5){
    decimals_last_day <- decimals
  } else{
    decimals_last_day <- decimals-1
  }
  C_t[1] <- round(Sumut)
  
  for (t in 2:N_DAYS){
    ### check whether NPI_k gets active at t
    # if NPI gets active it stays active for a random number of days between 60 to 120 days
    for (npi in 1:5) {
      NPI <- paste0("state_NPI",npi)
      if (!get(NPI)[t]) {
        # p_set_active <- sigmoid(Hicu_t[t-1], b0=bb0, b1=bb1)
        p_set_active <- sigmoid(Hicu_t[t-1], b0=bb0, b1=get(paste0("bb",npi)))
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
    c1_t <- sum(I_t[1:(t-1)] / POPULATION) * (1-p_reinf)
    
    c2_t <- (vaccinations1[t]) * vacc_on_R1 +
      (vaccinations2[t]) * vacc_on_R2
    
    cf_full <- (1 - c1_t - c2_t * (1-c1_t))
    
    Rt <- Rt * cf_full
    
    c1_t_vec[t] <- c1_t
    c2_t_vec[t] <- c2_t
    cfull_t_vec[t] <- cf_full
    R_t_vec[t] <- Rt # just here to keep track and monitor it
    
    # given final Rt, finally calculate I_t and Hicu_t
    
    ### calc I_t
    # slow version
    # sumut <- 0
    # for (u in 1:(t-1)) {
    #   sumut <- sumut + I_t[u] * gamma[t-u+1]
    # }
    # I_t[t] <- rnbinom(1, mu = sumut*Rt, size = phi_inf)
    I_t[t] <- calc_one_t_cpp(t, I_t[1:t], gamma, Rt, phi_inf)
    
    #### generate Ct
    
    Sumut <- decimals_last_day
    for(u in 1:t){
      Sumut <- Sumut + I_t[u] * Xi_C[t-u+1]
    }
    decimals <- Sumut-floor(Sumut)
    if(decimals<0.5){
      decimals_last_day <- decimals
    } else{
      decimals_last_day <- decimals-1
    }
    C_t[t] <- round(Sumut)
    
    ### calc Hicu_t
    #### generate Hicut
    # slow version
    # Sumut <- 0
    # for(u in 1:t){
    #   Sumut <- Sumut + C_t[u] * Xi_Hicu[t-u+1]
    # }
    # Hicu_t[t] <- rnbinom(1, mu = Sumut*piHicu[t], size = phi_icu)
    Hicu_t[t] <- calc_one_t_cpp(t, C_t, Xi_Hicu, piHicu[t], phi_icu) # Ct is longer, but only the first t values should be used, nice sanity check
  }
  
  data[, NPI1 := as.numeric(state_NPI1)]
  data[, NPI2 := as.numeric(state_NPI2)]
  data[, NPI3 := as.numeric(state_NPI3)]
  data[, NPI4 := as.numeric(state_NPI4)]
  data[, NPI5 := as.numeric(state_NPI5)]
  
  data[, I_t := I_t]
  data[, C_t := C_t]
  data[, Hicu_t := Hicu_t]
  
  data[, R_t := R_t_vec]
  data[, cf1 := c1_t_vec]
  data[, cf2 := c2_t_vec]
  data[, cf_full := cfull_t_vec]
  
  # just to check whether the simulation is doing some reasonable stuff
  # layout(1:2)
  # I_t %>% plot
  # R_t_vec %>% plot
  
  
  
  #### generate repCt
  
  # we need the weekdays since the model must know which column to pick from Xi_repC
#   data[, weekday := lubridate::wday(date, week_start = T)]
  data[, weekday := lubridate::wday(date, week_start = 1)]
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
      Sumut <- Sumut + data$C_t[u] * Xi_repC[[wd]][t-u+1]
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
  
  data[, D_t := D_t]
  
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
  data[, H_t := H_t]
  
  # add further required colmuns
  data[, ifr_t_m := piD]
  data[, first_vaccination := vaccinations1]
  data[, second_vaccination := vaccinations2]
  data[, alpha := proportion_alpha]
  data[, delta := proportion_delta]
  
  return(data[])
}



# test it for one country
data <- simulate_data_country(123)
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
  geom_line(aes(date, I_t))+
  geom_line(aes(date, C_t), col="darkgrey")+
  geom_line(aes(date, repC_t), col="darkorange")+
  geom_line(aes(date, D_t), col="darkred")+
  geom_line(aes(date, Hicu_t), col="darkblue")+
  geom_line(aes(date, H_t), col="blue")
g2 <- ggplot(data)+
  geom_line(aes(date, NPI1),col=2)+
  geom_line(aes(date, NPI2),col=3)+
  geom_line(aes(date, NPI3),col=4)+
  geom_line(aes(date, NPI4),col=5)+
  geom_line(aes(date, NPI5),col=6)
egg::ggarrange(g1,g2)



### use simulate_data_country function to generate N_COUNTRY countries
generate_data <- function(seed=321, write=F, path=paste0(PATH_DATA, "data_simulated.csv")){
  set.seed(seed)
  seeds <- sample(1:9999, N_COUNTRY)
  dds <- mapply(function(cc, seed){
    dd <- simulate_data_country(seed)
  } ,
  cc = paste0("country", LETTERS[1:N_COUNTRY]),
  seed = seeds,
  SIMPLIFY = F
  ) %>% 
    rbindlist(idcol = "country")
  if (write){
    fwrite(dds, path)
  }
  return(dds)
}


# test for many countries
d_tmp <- generate_data()
ggplot(d_tmp)+
  geom_line(aes(date, I_t))+
  geom_line(aes(date, C_t), col="darkgrey")+
  geom_line(aes(date, repC_t), col="darkorange")+
  geom_line(aes(date, D_t), col="darkred")+
  geom_line(aes(date, Hicu_t), col="darkblue")+
  geom_line(aes(date, H_t), col="blue")+
  facet_wrap(~country, scales = "free")


### generate N_DATASETS datasets
tt <- Sys.time()

set.seed(123)
seeds <- sample(1:9999, N_DATASETS)
datasets <- parallel::mcmapply(function(ds, seed){
  generate_data(seed, write=T, path=paste0(PATH_DATA, "/simulated_data/standard/",ds)) 
}, 
ds=paste0("data_sim_5NPIs_", 1:N_DATASETS, ".csv"),
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





