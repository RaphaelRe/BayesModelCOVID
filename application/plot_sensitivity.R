library(data.table)
library(magrittr)
library(ggplot2)
library(egg)
library(scales)

CORES <- 50 
PATH_RESULTS = "results/sensitivity_analysis/"
PATH_DATA = "../data/real_data.csv"
PATH_PLOTS = "plots/sensitivity_analysis/"  # assumes an existing directory


INTERVENTIONS = list.files(paste0(PATH_RESULTS, "high_cont/")) %>% grep("alpha",.) %>% 
  (function(x)list.files(paste0(PATH_RESULTS, "high_cont/"))[x]) %>% 
  strsplit("_") %>% sapply(function(x) x[4]) %>% 
  unique() %>% setdiff("alpha.txt")


# get all country names
ccs = fread(PATH_DATA)$country %>% unique()

# reads for a parameter all chains apply optional thinning and returns a meldted data frame
get_chains_setting = function(parameter, setting, thin = 1){ # thinning here not necessary as chains are already thinned
  files <- list.files(paste0(PATH_RESULTS, "/", setting))
  l <- grepl(paste0(parameter), files)
  data = lapply(files[l], function(y){
    return(read.csv(paste0(PATH_RESULTS, setting, "/", y), header = F)[[1]])
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


get_chains <- function(parameter){
  settings <- c("high_cont", "high_ifr", "high_vacc", "high_voc", "low_cont", "low_ifr", "low_vacc", "low_voc")
  sapply(settings, function(setting) get_chains_setting(parameter, setting), simplify = F, USE.NAMES = T) %>% 
    rbindlist(idcol = "setting")
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



## NPI effects

plot_effects_setting <- function(s){
  # posterior
  df_alphas <- results_alpha[get("setting") == s & country=="mean"]
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
  predictive_npis_mean <- results_alpha[get("setting") == s & country=="mean"]
  predictive_npis_sd <- results_alpha[get("setting") == s & country=="sd"]
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
  return(df_moments_both)
}

ests <- c("high_cont", "high_ifr", "high_vacc", "high_voc", "low_cont", "low_ifr", "low_vacc", "low_voc") %>% 
  sapply(plot_effects_setting, simplify = F, USE.NAMES = T) %>% 
  rbindlist(idcol = "setting")


ests[, setting_new := dplyr::case_when(
  ests$setting == "low_cont" ~ "Decreased contagiousness VOC",
  ests$setting == "high_cont" ~ "Increased contagiousness VOC",
  ests$setting == "low_ifr" ~ "Decreased IFR",
  ests$setting == "high_ifr" ~ "Increased IFR",
  ests$setting == "low_vacc" ~ "Decreased Effect vaccination",
  ests$setting == "high_vacc" ~ "Increased Effect vaccination",
  ests$setting == "low_voc" ~ "Decreased Effect VOC",
  ests$setting == "high_voc" ~ "Increased Effect VOC",
  TRUE ~ "-1"
)]


ests[, setting_fac := factor(setting_new, levels = c(
  "Decreased IFR", "Increased IFR",
  "Decreased Effect vaccination", "Increased Effect vaccination",
  "Decreased Effect VOC", "Increased Effect VOC",
  "Decreased contagiousness VOC", "Increased contagiousness VOC"
  ))]

trafo = function(x) x
trafo_red = function(x) 1-exp(-x) # reduction in R0
trafo_inc = function(x) exp(-x)-1 # increase in R0

(ggplot(ests[!(parameter %in% c("Winter", "Spring", "Autumn"))],
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
    ylab("Reduction in R")+
    facet_wrap(~setting_fac, ncol = 2)) %>% ggsave(paste0(PATH_PLOTS, "mean_effects_npis_sensitivity.pdf"), .,  width = 8, height = 10)



(ggplot(ests[(parameter %in% c("Winter", "Spring", "Autumn"))],
        aes(x=as.factor(parameter), color = type))+
    geom_point(aes(y=trafo_inc(mean)),position=position_dodge(.3), size=4., shape=18)+
    geom_errorbar(aes(ymin=trafo_inc(ql), ymax=trafo_inc(qu)), linewidth=3, width=0, position=position_dodge(.3))+
    geom_errorbar(aes(ymin=trafo_inc(qll), ymax=trafo_inc(quu)), linewidth=1.6, width=0.0, position=position_dodge(.3))+
    geom_hline(yintercept = 0, linetype = "dotted")+
    scale_color_manual(values=c("orange2", "steelblue"))+
    ylim(-0.15, 0.8)+
    theme_bw()+
    theme(text = element_text(size = 14), axis.text.x = element_text(angle = 45, vjust=0.99, hjust=1))+
    theme(legend.position = "none", axis.title.x = element_blank(), legend.title = element_blank())+
    # ggtitle("Reduction in R (NPIs)")
    ylab("Reduction in R")+
    facet_wrap(~setting_fac, ncol = 2)) %>% ggsave(paste0(PATH_PLOTS,"mean_effects_seasons_sensitivity.pdf"), ., width = 8, height = 10)

