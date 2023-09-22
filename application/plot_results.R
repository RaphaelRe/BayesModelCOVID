library(data.table)
library(magrittr)
library(ggplot2)
library(egg)
library(scales)

CORES <- 50 

PATH_RESULTS = "results/main_results/"
PATH_DATA = "../data/real_data.csv"
PATH_PLOTS = "plots/main_results/"

if (!dir.exists(PATH_PLOTS)) dir.create(PATH_PLOTS, recursive = T)

# get all interventions
INTERVENTIONS = list.files(PATH_RESULTS) %>% grep("alpha",.) %>% 
  (function(x)(list.files(PATH_RESULTS)[x])) %>% 
  strsplit("_") %>% sapply(function(x) x[4]) %>% 
  unique() %>% setdiff("alpha.txt")

# get all country names
ccs = fread(PATH_DATA)$country %>% unique()

# reads for a parameter all chains apply optional thinning and returns a meldted data frame
get_chains = function(parameter, thin = 1){ # thinning here not necessary as chains are already thinned
  files <- list.files(PATH_RESULTS)
  l <- grepl(paste0(parameter), files)
  data = lapply(files[l], function(y){
    return(read.csv(paste0(PATH_RESULTS, y), header = F)[[1]])
  })
  if (length(data) == 0){
    paste0(parameter, " not available") %>% print
    return(data.table(chain=NA, value=NA, i=NA, parameter=parameter))
  }
  its = sapply(data, length)
  ids = lapply(its, function(x) seq(1, x, by = thin))
  data = mapply(function(d,ids) d[ids], d = data, ids = ids, SIMPLIFY = F) %>% 
    lapply(data.table) %>% 
    set_names(paste0("chain_", 1:length(its))) %>% rbindlist(idcol = "chain")
  
  iis <- sapply(ids, length)
  data[, i:=  sapply(iis, function(i) 1:i) %>% c %>%  unlist()]
  data[, parameter:= parameter]
  setnames(data, "V1", "value")
  #data = melt(data, id.vars = c("i", "parameter")) %>% set_names(c("i", "parameter","chain", "value"))
  return(data[])
}



### get all data

# make cluster for faster reading
cl = parallel::makeForkCluster(CORES)

# alpha
results_alpha <- expand.grid(c("mean","sd", ccs), INTERVENTIONS) %>% data.table
data.table::setorder(results_alpha,"Var1")
results_alpha <- Reduce(function(x,y) paste0(y,"_",x),results_alpha) %>% pbapply::pblapply(get_chains, cl = cl) %>% rbindlist()

# must use "by" and lapply since ':=' cannot be used here (probably because the same value appears too often at start) 
results_alpha[, intervention := lapply(parameter, function(p) strsplit(p, "_")[[1]][1]), by = 1:nrow(results_alpha)]
results_alpha[, country := lapply(parameter, function(p) strsplit(p, "_")[[1]][2]), by = 1:nrow(results_alpha)]


# rho
periods = fread(PATH_DATA) %>% extract2("rho_period") %>% unique

results_rho = lapply(ccs, function(x){paste0("rho_", x, periods, ".txt")}) %>% unlist() %>% 
  pbapply::pbsapply(get_chains, simplify = F,cl = cl) %>% rbindlist

results_rho[, period := lapply(parameter, function(p) gsub("[^[:digit:]]", "",  p)), by = 1:nrow(results_rho)]
results_rho[, country := lapply(parameter, function(p) gsub('[[:digit:]]+', '',gsub('.{4}$', '',strsplit(p, "_")[[1]][2]))), by = 1:nrow(results_rho)]

d_tmp <- fread(PATH_DATA)[,c("date", "rho_period")] %>% 
  set_names(c("date", "period")) %>% unique()
results_rho$date <- d_tmp[chmatch(results_rho$period, as.character(d_tmp$period)), date]

# beta_VOCs
results_beta_vocs = rbindlist(list(alpha = get_chains("beta_alpha"), delta = get_chains("beta_delta")), idcol="type")

parallel::stopCluster(cl)



# rename alphas to normal names
results_alpha[intervention == "closingSchools", intervention := "School closure"]
results_alpha[intervention == "restrictGatherings", intervention := "Restrictions on gatherings"]
results_alpha[intervention == "lockdown", intervention := "Lockdown"]
results_alpha[intervention == "subsequentLockdown", intervention := "Subsequent lockdown"]
results_alpha[intervention == "generalBehavioralChanges", intervention := "General behavioral changes"]

results_alpha[intervention == "seasonAutumn", intervention := "Autumn"]
results_alpha[intervention == "seasonWinter", intervention := "Winter"]
results_alpha[intervention == "seasonSpring", intervention := "Spring"]

results_alpha[, intervention_factor := factor(intervention, levels = c("School closure",
                                                                       "General behavioral changes",
                                                                       "Restrictions on gatherings",
                                                                       "Lockdown",
                                                                       "Subsequent lockdown",
                                                                       "Spring",
                                                                       "Autumn",
                                                                       "Winter"
))]


### results NPIs

## traceplots for alphas
(results_alpha[country == "mean"] %>%
    ggplot(aes(i, value, color = chain))+ 
    geom_line(alpha = 0.3)+
    facet_wrap(~intervention_factor,4,2,scales = "free")) %>% ggsave(paste0(PATH_PLOTS, "/trace_plots_alpha.pdf"),.,width = 7 , height = 10)


## NPI effects

# posterior
df_alphas <- results_alpha[country=="mean"]
# trafo = function(x) 1-exp(-x)
trafo = function(x) x
trafo_red = function(x) 1-exp(-x) # reduction in R0
trafo_inc = function(x) exp(-x)-1 # increase in R0

df_alphas[,value_t := trafo(value)]
m = df_alphas[, mean(value_t), by = .(intervention)] %>% set_names(c("parameter","mean"))
ql = df_alphas[, quantile(value_t, 0.25, na.rm = T), by = .(intervention)] %>% set_names(c("parameter","ql"))
qu = df_alphas[, quantile(value_t, 0.75, na.rm = T), by = .(intervention)]%>% set_names(c("parameter","qu"))
qll = df_alphas[, quantile(value_t, 0.025, na.rm = T), by = .(intervention)]%>% set_names(c("parameter","qll"))
quu = df_alphas[, quantile(value_t, 0.975, na.rm = T), by = .(intervention)]%>% set_names(c("parameter","quu"))


df_moments = Reduce(merge,list(m=m,ql=ql,qll=qll,qu=qu,quu=quu))
setorder(df_moments, -parameter)


# predictive effects
predictive_npis_mean <- results_alpha[country=="mean"]
predictive_npis_sd <- results_alpha[country=="sd"]
setnames(predictive_npis_sd, "value", "value_sd")
predictive_npis <- merge(predictive_npis_mean, predictive_npis_sd[, -c("parameter", "country")], by = c("chain", "i", "intervention"))
predictive_npis[, pred_val := rnorm(1, value, value_sd), by = 1:nrow(predictive_npis)]
setorder(predictive_npis, parameter, chain)


predictive_npis[,pred_val_trandsformed := trafo(pred_val)]
m = predictive_npis[, mean(pred_val_trandsformed), by = .(intervention)] %>% set_names(c("parameter","mean"))
ql = predictive_npis[, quantile(pred_val_trandsformed, 0.25, na.rm = T), by = .(intervention)] %>% set_names(c("parameter","ql"))
qu = predictive_npis[, quantile(pred_val_trandsformed, 0.75, na.rm = T), by = .(intervention)]%>% set_names(c("parameter","qu"))
qll = predictive_npis[, quantile(pred_val_trandsformed, 0.025, na.rm = T), by = .(intervention)]%>% set_names(c("parameter","qll"))
quu = predictive_npis[, quantile(pred_val_trandsformed, 0.975, na.rm = T), by = .(intervention)]%>% set_names(c("parameter","quu"))

df_moments_pred = Reduce(merge,list(m=m,ql=ql,qll=qll,qu=qu,quu=quu))
setorder(df_moments_pred, -parameter)


# combine plots from posterior mean effects and predictive ones

# setnames(df_moments_pred, names(df_moments_pred)[-1], paste0(names(df_moments_pred)[-1], "_pred"))
df_moments_both <- rbindlist(list(mean = df_moments, pred = df_moments_pred), idcol = "type") 
df_moments_both[, parameter := factor(parameter, 
                                      levels=c("School closure", "General behavioral changes",
                                               "Restrictions on gatherings", "Lockdown", "Subsequent lockdown",
                                               "Spring", "Autumn", "Winter"))]


# merge(df_moments, df_moments_pred) %>%
# npis
(g1 <- ggplot(df_moments_both[!(parameter %in% c("Winter", "Spring", "Autumn"))],
              aes(x=as.factor(parameter), color = type))+
    geom_point(aes(y=trafo_red(mean)),position=position_dodge(.3), size=4., shape=18)+
    geom_errorbar(aes(ymin=trafo_red(ql), ymax=trafo_red(qu)), linewidth=3, width=0, position=position_dodge(.3))+
    geom_errorbar(aes(ymin=trafo_red(qll), ymax=trafo_red(quu)), linewidth=1.6, width=0.0, position=position_dodge(.3))+
    geom_hline(yintercept = 0, linetype = "dotted")+
    scale_color_manual(values=c("orange2", "steelblue"))+
    ylim(-0.15, 0.8)+
    theme_bw()+
    theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, vjust=0.99, hjust=1))+
    theme(legend.position = "none", axis.title.x = element_blank(), legend.title = element_blank())+
    # ggtitle("Reduction in R (NPIs)")
    ylab("Reduction in R")
)  %>% ggsave(paste0(PATH_PLOTS, "/mean_effects_preds_est_quantiles.pdf"),., width = 7, height = 5)




# seasons
(g2 <- ggplot(df_moments_both[(parameter %in% c("Winter", "Spring", "Autumn"))],
              aes(x=as.factor(parameter), color = type))+
    geom_point(aes(y=trafo_inc(mean)),position=position_dodge(.3), size=4., shape=18)+
    geom_errorbar(aes(ymin=trafo_inc(ql), ymax=trafo_inc(qu)), linewidth=3, width=0, position=position_dodge(.3))+
    geom_errorbar(aes(ymin=trafo_inc(qll), ymax=trafo_inc(quu)), linewidth=1.6, width=0.0, position=position_dodge(.3))+
    geom_hline(yintercept = 0, linetype = "dotted")+
    scale_color_manual(values=c("orange2", "steelblue"))+
    ylim(-0.15, 0.8)+
    theme_bw()+
    theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, vjust=0.99, hjust=1))+
    theme(legend.position = "bone", axis.title.x = element_blank(), legend.title = element_blank())+
    # ggtitle("Increase in R")+
    ylab("Increase in R")
)  %>% ggsave(paste0(PATH_PLOTS, "/mean_effects_seasons_preds_est_quantiles.pdf"),., width = 7, height = 5)

# egg::ggarrange(g1,g2, ncol = 2)# %>% ggsave(paste0(PATH_PLOTS, "/mean_effects_npis_seasons_preds_est_quantiles_2in1.pdf"),., width = 12, height = 7)



## as table
# tmpdf1 <- df_moments_both[type == "mean", .(parameter, mean, qll, quu)][c(5,6,7,1,3,4,8,2)]
# tmpdf2 <- df_moments_both[type == "pred", .(parameter, mean, qll, quu)][c(5,6,7,1,3,4,8,2)]
# 
# # apply trafo
# tmpdf1[1:5, c("mean", "qll", "quu") := list(trafo_red(mean), trafo_red(qll), trafo_red(quu))]
# tmpdf1[6:8, c("mean", "qll", "quu") := list(trafo_inc(mean), trafo_red(qll), trafo_red(quu))]
# tmpdf2[1:5, c("mean", "qll", "quu") := list(trafo_red(mean), trafo_red(qll), trafo_red(quu))]
# tmpdf2[6:8, c("mean", "qll", "quu") := list(trafo_inc(mean), trafo_red(qll), trafo_red(quu))]
# 
#
# for latex
# data.table(
#   NPI = tmpdf1$parameter,
#   Mean = paste0(round(tmpdf1$mean, 3), " [",round(tmpdf1$ql, 3), ", ",round(tmpdf1$quu, 3), "]"),  
#   Pred = paste0(round(tmpdf2$mean, 3), " [",round(tmpdf2$ql, 3), ", ",round(tmpdf2$quu, 3), "]")  
# ) %>% 
#   xtable::xtable(digits=3) %>% 
#   print(include.rownames=FALSE)



# for all countries all effects
infos_all_alphas <- lapply(ccs, function(cc){
  df_alphas <- results_alpha[country==cc]
  trafo = function(x) 1-exp(-x)
  # trafo = function(x) x
  df_alphas[,value_t := trafo(value)]
  m = df_alphas[, mean(value_t), by = .(intervention)] %>% set_names(c("parameter","mean"))
  ql = df_alphas[, quantile(value_t, 0.25, na.rm = T), by = .(intervention)] %>% set_names(c("parameter","ql"))
  qu = df_alphas[, quantile(value_t, 0.75, na.rm = T), by = .(intervention)]%>% set_names(c("parameter","qu"))
  qll = df_alphas[, quantile(value_t, 0.025, na.rm = T), by = .(intervention)]%>% set_names(c("parameter","qll"))
  quu = df_alphas[, quantile(value_t, 0.975, na.rm = T), by = .(intervention)]%>% set_names(c("parameter","quu"))
  df_moments = Reduce(merge,list(m=m,ql=ql,qll=qll,qu=qu,quu=quu))
  setorder(df_moments, parameter)
}) %>%
  set_names(ccs) %>% 
  rbindlist(idcol = "country") %>% 
  .[, parameter_fac := factor(parameter, levels = c("School closure", "General behavioral changes",
                                                    "Restrictions on gatherings", "Lockdown", "Subsequent lockdown",
                                                    "Spring", "Autumn", "Winter"))]

(ggplot(infos_all_alphas) + 
    geom_errorbar(aes(y = parameter_fac, xmin=qll, xmax=quu, width = .15), linewidth=2 , alpha = 0.9,color = "steelblue3")+
    geom_errorbar(aes(y = parameter_fac, xmin=ql, xmax=qu, width = .05), linewidth=3.5, alpha = 0.9,color = "steelblue4")+
    geom_point(aes(y = parameter_fac, x = mean),fill = "white", shape = 21, size = 2)+
    xlab("")+ylab("") + 
    scale_x_continuous(labels = percent)+
    geom_vline(xintercept = 0, color = "darkred", linetype = "dotted")+
    ggplot2::ggtitle("")+theme_bw()+coord_flip()+
    theme(text = element_text(size = 10), axis.text.x = element_text(angle = 45, vjust=0.99, hjust=1))+
    facet_wrap(~country, scales = "free_y")
) %>% ggsave(paste0(PATH_PLOTS, "/effects_all_countries_quantiles.pdf"), ., width = 14, height = 8)





###  results rho
rho_means <- results_rho[, .(m = mean(value)), by = c("country", "date", "period")]
rho_sd <- results_rho[, .(StdDev = sd(value)), by = c("country", "date", "period")]
rho_m_sd <- merge(rho_means, rho_sd)

d_tmp <- fread(PATH_DATA)[, .(country, date, rho_period)]
setnames(d_tmp, "rho_period", "period")
rho_m_sd[, period := as.numeric(period)]
d_tmp <- merge(d_tmp, rho_m_sd, by = c("country", "period"), all=T)
d_tmp[, date.y := NULL]
setnames(d_tmp, "date.x", "date")


## all rhos
egg::ggarrange(plots = 
                 lapply(ccs, function(cc, d_tmp){
                   ggplot(d_tmp[country == cc]) +
                     geom_line(aes(date, m, group=as.factor(m), color = StdDev), linewidth=2, lineend="round")+
                     ylab("reporting_ratio")+
                     ggtitle(cc)+
                     theme_bw()
                 }, d_tmp=d_tmp)
) %>% ggsave(paste0(PATH_PLOTS, "/rho_period_all_countries.pdf"), ., width = 12, height = 12)

# info about testing
data_testing <- fread("https://opendata.ecdc.europa.eu/covid19/testing/csv/data.csv")
todate <- function(yw){
  ISOweek::ISOweek2date(paste0(yw, "-4"))
}
data_testing[, date := todate(year_week)]


##  rhos with testing info
# France
rhos_France <- d_tmp[country == "France"]
rhos_France[, date := as.Date(date)]
testing_France <- data_testing[country == "France" & region_name == "France", .(date, positivity_rate)]


g_fr <- (merge(rhos_France, testing_France, all.x = T, by = "date") %>% 
    .[positivity_rate == 0, positivity_rate := NA] %>% 
    ggplot() +
    geom_line(aes(date, m, group=as.factor(m), color = StdDev), linewidth=2, lineend="round")+
    geom_point(aes(date, 1/(positivity_rate/100)/100), color = "orange", alpha = 0.8)+
    ylab("Case detection ratio")+
    ggtitle("France")+
    theme_bw()+
    scale_y_continuous(sec.axis = sec_axis(~.*100, name = "Inverse test positivity rate")))# %>% ggsave(paste0(PATH_PLOTS, "/rho_period_inv_testpos_rate_France.pdf"), ., width = 8, height = 5)


# Hungary
rhos_Hungary <- d_tmp[country == "Hungary"]
rhos_Hungary[, date := as.Date(date)]
testing_Hungary <- data_testing[country == "Hungary" & region_name == "Hungary", .(date, positivity_rate)]


g_hun <- (merge(rhos_Hungary, testing_Hungary, all.x = T, by = "date") %>% 
    ggplot() +
    geom_line(aes(date, m, group=as.factor(m), color = StdDev), linewidth=2, lineend="round")+
    geom_point(aes(date, 1/(positivity_rate/100)/200), color = "orange", alpha = 0.8)+
    ylab("Case detection ratio")+
    ggtitle("Hungary")+
    theme_bw()+
    scale_y_continuous(sec.axis = sec_axis(~.*200, name = "Inverse test positivity rate")))# %>% ggsave(paste0(PATH_PLOTS, "/rho_period_inv_testpos_rate_Hungary.pdf"), ., width = 8, height = 5)

ggpubr::ggarrange(g_hun, g_fr, labels = c("A", "B")) %>% ggsave(paste0(PATH_PLOTS, "/rho_period_inv_testpos_rate_hungary_france.pdf"),., width = 14, height = 5)







### prior - posterior plots

plot_posterior <- function(d, add_prior = TRUE, generator=NULL, xlim=c(0.4,0.8),
                           main="Prior / Posterior"){
  g <- ggplot(d)+
    geom_density(aes(value), fill = "steelblue", color= "steelblue", alpha = 0.5)
  
  if (add_prior){
    dd <- generator(get("d"))
    g <- g + geom_ribbon(aes(dd$grid, ymax=dd$V2, ymin=0), fill = "steelblue2",
                         color = "steelblue2", alpha=0.4)
  }
  g+xlim(xlim)+
    ggtitle(main)
}

generator <- function(mu, sd, min, max, d){
  grid <- seq(min,max, len = length(d$value))
  data.table(grid, dnorm(grid, mu, sd))
}

plot_posterior(results_beta_vocs[type == "alpha"], TRUE, 
               generator = function(d) generator(0.6,0.01,0.5,0.65,d),
               xlim=c(0.5,0.65)) %>% ggsave(paste0(PATH_PLOTS, "/prior_posterior_VOC_alpha.pdf"), ., width = 7, height = 4)

plot_posterior(results_beta_vocs[type == "delta"], TRUE,
               generator = function(d) generator(1.5,0.02, 0.8,1.8,d), 
               xlim=c(1.4,1.6)) %>% ggsave(paste0(PATH_PLOTS, "/prior_posterior_VOC_delta.pdf"), ., width = 7, height = 4)




### posterior predictions time series

# takes a country and a series neame (e.g. reported_cases) and returns a list with all series (i.e. results from 4 chains)
get_all_predictions <- function(country, series){
  files <- list.files(PATH_RESULTS)
  l <- grepl(paste0("predictions_", series, "_", country), files)
  data = lapply(files[l], function(y){
    return(fread(paste0(PATH_RESULTS,y), header = F))
  })  
  
  return(data)
}


# assumes a time point in eah row, each column is a  sample point from its distribution
extract_info <- function(d, thin = T, thin_n = 2){
  if (thin){
    d <- as.matrix(d)[,seq(1, ncol(d), by = thin_n)]
  }
  m <- apply(d, 1, mean)
  q_ll <- apply(d, 1, quantile, probs = c(0.025))
  q_l <- apply(d, 1, quantile, probs = c(0.25))
  q_uu <- apply(d, 1, quantile, probs = c(0.975))
  q_u <- apply(d, 1, quantile, probs = c(0.75))
  df <- data.frame(m = m, q_ll = q_ll,q_l = q_l, q_uu = q_uu, q_u = q_u)
  return(df)
}


# calls get_all_predictions, extracts info and expands it with the original data
get_all_data_country <- function(country_sel, series = "reported_cases"){
  # print(country_sel)
  df <- get_all_predictions(country_sel, series) %>%
    Reduce(function(x,y) cbind(x,y),.) %>%
    as.matrix() %>% 
    extract_info() %>% 
    data.table()
  df[, country:=country_sel]
  df[,index := 1:nrow(df)]
  true_ts <- fread(PATH_DATA)[,c("country", "date","repC_t", "D_t", "H_t", "Hicu_t")]
  true_ts <- true_ts[country == country_sel,]
  true_ts[, index := 1:nrow(true_ts)]
  data <- merge(df, true_ts, all = T)
  start_date <- data$date %>% min(na.rm = T)
  dates <- seq(start_date, length.out = nrow(data), by = "day")
  data[, date := dates]
  return(data)
}



# data
cl = parallel::makeForkCluster(6)
all_reported_cases <- parallel::parLapply(cl, ccs, get_all_data_country) %>% rbindlist()
all_data_deaths <- parallel::parLapply(cl, ccs, get_all_data_country, series = "deaths") %>% rbindlist()
all_data_infections <- parallel::parLapply(cl, ccs, get_all_data_country, series = "infections") %>% rbindlist()
all_data_hospitalizations <- parallel::parLapply(cl, ccs, get_all_data_country, series = "hospitalizations") %>% rbindlist()
all_data_icu <- lapply(ccs, get_all_data_country, series = "intensiveCare") %>% rbindlist()
all_data_cases <- parallel::parLapply(cl, ccs, get_all_data_country, series = "cases") %>% rbindlist()
parallel::stopCluster(cl)



#plot function: plots the posterior predictions with CI and an original time series
plot_predictions <- function(dat, country_sel = "Germany", orig_series = c("cases", "deaths","hospitalizations","intensiveCare", "None"), x_lab = 'Day', y_lab = element_blank()){
  orig_series = match.arg(orig_series)
  if (country_sel == "all"){
    g <- ggplot(dat, aes(x = date)) + 
      geom_ribbon(aes(ymin = q_l, ymax = q_u), fill = "steelblue3", alpha = 0.3)+
      geom_ribbon(aes(ymin = q_ll, ymax = q_uu), fill = "steelblue2", alpha = 0.2)+
      geom_line(aes(y= m), color = "royalblue4") + 
      xlab("Days") + ylab(y_lab)+
      ggtitle(country_sel) + theme_bw()+
      facet_wrap(country)
  } else {
    g <- ggplot(dat[country ==country_sel], aes(x = date)) + 
      geom_ribbon(aes(ymin = q_l, ymax = q_u), fill = "steelblue3", alpha = 0.3)+
      geom_ribbon(aes(ymin = q_ll, ymax = q_uu), fill = "steelblue2", alpha = 0.2)+
      geom_line(aes(y= m), color = "royalblue4") + 
      xlab("Days") + ylab(y_lab)+
      ggtitle(country_sel) + theme_bw()
  }
  if (orig_series == "cases"){
    print("Using cases as original series")
    return(g + geom_line(aes(y = repC_t), alpha = 0.7, size = 0.3, alpha = 0.7))
  } else if (orig_series == "hospitalizations"){
    print("Using hospitalizations as original series")
    return(g + geom_line(aes(y = H_t), alpha = 0.7, size = 0.3, alpha = 0.7))
  }  else if (orig_series == "intensiveCare"){
    print("Using icu as original series")
    return(g + geom_line(aes(y = Hicu_t), alpha = 0.7, size = 0.3, alpha = 0.7))
  } else if (orig_series == "deaths"){
    print("Using deaths as original series")
    return(g + geom_line(aes(y = D_t), alpha = 0.7, size = 0.3, alpha = 0.7))
  } else {
    print("Using no original series")
    g
  }
}


## plot all time series in one plot

# icu
lapply(ccs, function(cc){
  plot_predictions(all_data_icu, country_sel = cc, orig_series = "intensiveCare")
}) %>% 
  egg::ggarrange(plots = .) %>% ggsave(paste0(PATH_PLOTS, "/icu_all.pdf"),., width = 12, height = 12)

# hospitalizations
lapply(ccs, function(cc){
  plot_predictions(all_data_hospitalizations, country_sel = cc, orig_series = "hospitalizations")
}) %>% 
  egg::ggarrange(plots = .) %>% ggsave(paste0(PATH_PLOTS, "/hosp_all.pdf"),., width = 12, height = 12)

# deaths
lapply(ccs, function(cc){
  plot_predictions(all_data_deaths, country_sel = cc, orig_series = "deaths")
}) %>% 
  egg::ggarrange(plots = .) %>% ggsave(paste0(PATH_PLOTS, "/deaths_all.pdf"),., width = 12, height = 12)

# reported cases
lapply(ccs, function(cc){
  plot_predictions(all_reported_cases, country_sel = cc, orig_series = "cases")
}) %>% 
  egg::ggarrange(plots = .) %>% ggsave(paste0(PATH_PLOTS, "/repCases_all.pdf"),., width = 12, height = 12)

# infections
lapply(ccs, function(cc){
  plot_predictions(all_data_infections, country_sel = cc, orig_series = "cases")
}) %>% 
  egg::ggarrange(plots = .) %>% ggsave(paste0(PATH_PLOTS, "/infections_all.pdf"),., width = 12, height = 12)


## plot 3 time series for given countries in main manuscript
# France
g1 <- plot_predictions(all_reported_cases, country_sel = "France", orig_series = "cases")+
  ggtitle("France", subtitle = "Reported cases")+ theme(plot.title = element_text(hjust = 0.5))
g2 <- plot_predictions(all_data_deaths, country_sel = "France", orig_series = "deaths")+
  ggtitle("", subtitle = "Deaths")
g3 <- plot_predictions(all_data_hospitalizations, country_sel = "France", orig_series = "hospitalizations")+
  ggtitle("", subtitle = "Hospital occupancy")
g_france <- egg::ggarrange(g1,g2,g3, ncol=1) #%>% ggsave(paste0(PATH_PLOTS, "/cases_hosp_deaths_France.pdf"),., width = 5, height = 12)

# Hungary
g1 <- plot_predictions(all_reported_cases, country_sel = "Hungary", orig_series = "cases")+
  ggtitle("Hungary", subtitle = "Reported cases") + theme(plot.title = element_text(hjust = 0.5))
g2 <- plot_predictions(all_data_deaths, country_sel = "Hungary", orig_series = "deaths")+
  ggtitle("", subtitle = "Deaths")
g3 <- plot_predictions(all_data_hospitalizations, country_sel = "Hungary", orig_series = "hospitalizations")+
  ggtitle("", subtitle = "Hospital occupancy")
g_hungary <- egg::ggarrange(g1,g2,g3, ncol=1) #%>% ggsave(paste0(PATH_PLOTS, "/cases_hosp_deaths_Hungary.pdf"),., width = 5, height = 12)

ggpubr::ggarrange(g_hungary, g_france, labels = c("A", "B")) %>% 
  ggsave(paste0(PATH_PLOTS, "/ts_hugary_france_both.pdf"),., width = 12, height = 14)





## for all countries in supplement
plot_all_series <- function(country_sel){
  g1 <- plot_predictions(all_reported_cases, country_sel = country_sel, orig_series = "cases")+
    ggtitle(country_sel, subtitle = "Reported cases")+ theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  g2 <- plot_predictions(all_data_deaths, country_sel = country_sel, orig_series = "deaths")+
    ggtitle("", subtitle = "Deaths")
  g3 <- plot_predictions(all_data_hospitalizations, country_sel = country_sel, orig_series = "hospitalizations")+
    ggtitle("", subtitle = "Hospital occupancy")
  g4 <- plot_predictions(all_data_icu, country_sel = country_sel, orig_series = "intensiveCare")+
    ggtitle("", subtitle = "ICU occupancy")
  egg::ggarrange(g1,g2,g3,g4, ncol=1)
}


plot_4_sereis <- lapply(ccs, plot_all_series)

ggpubr::ggarrange(plotlist = plot_4_sereis[1:8], labels = LETTERS[1:8], ncol=4, nrow=2) %>% ggsave(paste0(PATH_PLOTS, "/8_series_fig1.pdf"),., width = 12, height = 14)
ggpubr::ggarrange(plotlist = plot_4_sereis[9:16], labels = LETTERS[9:16], ncol=4, nrow=2) %>% ggsave(paste0(PATH_PLOTS, "/8_series_fig2.pdf"),., width = 12, height = 14)
ggpubr::ggarrange(plotlist = plot_4_sereis[17:20], labels = LETTERS[17:20], ncol=4) %>% ggsave(paste0(PATH_PLOTS, "/8_series_fig3.pdf"),., width = 12, height = 7)


### plot prior on I_1_m using n = 100k iid samples
n <- 1e5
phi_I <- truncnorm::rtruncnorm(n, a = 0, mean = 0, sd = 0.015)
sigma_tau <- truncnorm::rtruncnorm(n, a = 0, mean = 0, sd = 10)
tau <- rgamma(n, 10, scale = 1)
tau_m <- truncnorm::rtruncnorm(n, a = 0, mean = tau, sd = sigma_tau)
I_1_m <- rnbinom(n, size = 1/sqrt(phi_I), mu = tau_m)

plot_prior_I1 <- ggplot()+
  geom_histogram(aes(I_1_m, y=after_stat(density)),bins = 50, fill = "steelblue", alpha = 0.6, color = "steelblue4")+
  theme_bw()+
  xlab(expression(I["1,m"]))

ggsave(paste0(PATH_PLOTS, "../prior_I1.pdf"), plot_prior_I1, width = 8, height = 4.2)


