library(data.table)
library(magrittr)
library(ggplot2)
library(egg)
library(scales)
library(cowplot)

CORES <- 50 # cores for parallelization

PATH_RESULTS = "results/standard/"
PATH_DATA_PART = "../data/simulated_data/standard/"
PATH_PLOTS = "plots/standard/"

if (!dir.exists(PATH_PLOTS)) dir.create(PATH_PLOTS, recursive = T)

N_datasets <- list.files(PATH_RESULTS) %>% length


# get all interventions
INTERVENTIONS = list.files(paste0(PATH_RESULTS, "res_dataset_1/")) %>% grep("alpha",.) %>% 
  (function(x)(list.files(paste0(PATH_RESULTS, "res_dataset_1/"))[x])) %>% 
  strsplit("_") %>% sapply(function(x) x[4]) %>% 
  unique() %>% setdiff("alpha.txt")

# get all countries
ccs = fread(paste0(PATH_DATA_PART, "/data_sim_5NPIs_1.csv"))$country %>% unique()


# reads for a parameter all chains apply thinning and returns data table with info
get_chains = function(parameter = "alpha_NPI1_mean", dataset="res_dataset_1",
                      start=1, thin = 1, path = getwd()){
  files <- list.files(paste0(PATH_RESULTS, dataset, "/"))
  l <- grepl(paste0(parameter), files)
  data = lapply(files[l], function(y){
    return(read.csv(paste0(paste0(PATH_RESULTS, dataset, "/"), y), header = F)[[1]])
  })
  if (length(data) == 0){
    paste0(parameter, " not available") %>% print
    return(data.table(chain=NA, value=NA, i=NA, parameter=parameter))
  }
  
  its = sapply(data, length)
  ids = lapply(its, function(x) seq(start, x, by = thin))
  data = mapply(function(d,ids) d[ids], d = data, ids = ids, SIMPLIFY = F) %>% 
    lapply(data.table) %>% 
    set_names(paste0("chain_", 1:length(its))) %>% rbindlist(idcol = "chain")
  
  iis <- sapply(ids, length)
  data[, i:=  sapply(iis, function(i) 1:i) %>% c %>%  unlist()]
  data[, parameter:= parameter]
  setnames(data, "V1", "value")
  return(data[])
}


get_chains_all_sets <- function(parameter = "alpha_NPI1_mean"){
  print(parameter)
  dsets <- paste0("res_dataset_", 1:N_datasets)
  chains <- lapply( dsets, function(ds) get_chains(parameter, ds))
  return(rbindlist(chains,  idcol = "dataset"))
}


# make cluster for faster reading
cl = parallel::makeForkCluster(CORES)

## alpha
results_alpha <- expand.grid(c("mean","sd", ccs), INTERVENTIONS) %>% data.table
data.table::setorder(results_alpha,"Var1")
results_alpha <- Reduce(function(x,y) paste0(y,"_",x),results_alpha) %>%
  pbapply::pblapply(get_chains_all_sets,cl=cl) %>% rbindlist()

# has to use "by" and lapply, cant use := directly 
results_alpha[, intervention := lapply(parameter, function(p) strsplit(p, "_")[[1]][1]), by = 1:nrow(results_alpha)]
results_alpha[, country := lapply(parameter, function(p) strsplit(p, "_")[[1]][2]), by = 1:nrow(results_alpha)]

parallel::stopCluster(cl)



####################################
# estimated effects and bias as table
estimated_mean_ds <- results_alpha[!(country %in% c("mean", "sd")), .(m = mean(value, na.rm = T)), by = c("dataset", "intervention", "country")]
estimated_ql_ds <- results_alpha[!(country %in% c("mean", "sd")), .(ql = quantile(value, 0.025, na.rm = T)), by = c("dataset", "intervention", "country")]
estimated_qu_ds <- results_alpha[!(country %in% c("mean", "sd")), .(qu = quantile(value, 0.975, na.rm = T)), by = c("dataset", "intervention", "country")]


calc_biases <- function(ds){
  # get true values
  true_vals_ds <- fread(paste0(PATH_DATA_PART, "/data_sim_5NPIs_",ds, ".csv"))[, c("country", paste0("true_val_NPI", 1:5))
  ][, lapply(.SD, unique), by=country] %>% 
    melt(id.vars="country", variable.factor = F) %>% 
    set_names(c("country", "intervention", "true_value"))
  true_vals_ds[, intervention := sapply(intervention, function(p) strsplit(p, "_")[[1]][3])]
  
  full_bias_ds <- merge(estimated_mean_ds[dataset==ds], true_vals_ds) %>%
    merge(estimated_ql_ds[dataset==ds]) %>% 
    merge(estimated_qu_ds[dataset==ds])
  
  full_bias_ds[,rel_bias := (m-true_value)/true_value, by = c("intervention", "country")]
  full_bias_ds[, inside := true_value > ql & true_value < qu]
  return(full_bias_ds)
}


cl <- parallel::makeForkCluster(CORES)
full_bias_ds <- pbapply::pblapply(1:N_datasets, calc_biases, cl = cl) %>% rbindlist
parallel::stopCluster(cl)


table_biases <- full_bias_ds[, .(true_mean = mean(true_value),
                                 mean_estimates = mean(m),
                                 mean_biasR = mean(rel_bias)*100,
                                 coverage = sum(inside / (N_datasets*10))*100
), by = c("intervention")]


# xtable::xtable(table_biases, digits=3) %>% 
#   print(include.rownames=FALSE)



## plot estimates as violin plot for 20 datasets
# this is the overall mean over all mean vlaues
# required for hline
overall_mean_true_values <- sapply(1:20, function(ds){
  d <- fread(paste0(PATH_DATA_PART, "/data_sim_5NPIs_",ds, ".csv"))
  d[, paste0("true_val_NPI", 1:5)] %>% sapply(mean) %>% t %>% data.table
}, simplify = F) %>% 
  rbindlist() %>% 
  sapply(mean)

hlinedf <- data.table(
  value = overall_mean_true_values[1:3],
  intervention = paste0("NPI", 1:3)
)

### plot as violin - only the first 20 datasets
(ggplot(results_alpha[country=="mean" & 
                        intervention %in% paste0("NPI", 1:3) &
                        dataset %in% 1:20
], 
aes(x=as.factor(dataset), y= value, fill = intervention)) + 
    geom_violin(trim=FALSE, alpha = 0.3)+
    labs(x="Dataset", y = "Posterior samples")+
    geom_boxplot(width=0.2)+
    facet_wrap(~intervention, scales = "free")+
    # add actual mean
    geom_hline(data=hlinedf, aes(yintercept=value), linetype="dashed")+
    theme_bw()+
    theme(legend.position = "none")
) %>%   ggsave(paste0(PATH_PLOTS, "/violin_mean_effects_standard.pdf"), ., width = 11, height = 3)



### plots posterior predictions


# assumes a time point in eah row, each column is a  sample point from its distribution
extract_info <- function(d, thin = T, thin_n = 1){
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


# takes a country and a series neame (e.g. reported_cases) and returns a list with all series (i.e. results from all chains)
get_all_predictions <- function(dataset, country, series = "reported_cases"){
  files <- list.files(paste0(PATH_RESULTS, "res_dataset_", dataset, "/"))
  l <- grepl(paste0("predictions_", series, "_", country), files)
  data = lapply(files[l], function(y){
    return(fread(paste0(PATH_RESULTS, "res_dataset_",dataset, "/", y), header = F))
  })  
  return(data)
}



# calls get_all_predictions, extracts info and appends the original data
get_all_data_country <- function(country_sel, dataset, series = "reported_cases"){
  print(country_sel)
  df <- get_all_predictions(dataset, country_sel, series) %>%
    Reduce(function(x,y) cbind(x,y),.) %>%
    as.matrix() %>% 
    extract_info() %>% 
    data.table()
  df[, country:=country_sel]
  df[,index := 1:nrow(df)]
  true_ts <- fread(paste0(PATH_DATA_PART, "data_sim_5NPIs_",dataset,".csv"))[,c("country", "date", "I_t", "C_t", "repC_t", "D_t", "H_t", "Hicu_t")]
  true_ts <- true_ts[country == country_sel,]
  true_ts[, index := 1:nrow(true_ts)]
  data <- merge(df, true_ts, all = T)
  start_date <- data$date %>% min(na.rm = T)
  dates <- seq(start_date, length.out = nrow(data), by = "day")
  data[, date := dates]
  return(data[])
}




# data
cl = parallel::makeForkCluster(CORES)
data_all_datasets <- pbapply::pblapply(paste0(1:10), function(dd){
  all_reported_cases <- lapply(ccs, get_all_data_country, dataset=dd) %>% rbindlist()
  all_data_deaths <- lapply(ccs, get_all_data_country, dataset=dd, series = "deaths") %>% rbindlist()
  all_data_infections <- lapply(ccs, get_all_data_country, dataset=dd, series = "infections") %>% rbindlist()
  all_data_hospitalizations <- lapply(ccs, get_all_data_country, dataset=dd, series = "hospitalizations") %>% rbindlist()
  all_data_icu <- lapply(ccs, get_all_data_country, dataset=dd, series = "intensiveCare") %>% rbindlist()
  all_data_cases <- lapply(ccs, get_all_data_country, dataset=dd, series = "cases") %>% rbindlist()
  d <- list(repCases = all_reported_cases,
            deaths = all_data_deaths,
            infections = all_data_infections,
            cases = all_data_cases,
            hospitalBeds = all_data_hospitalizations,
            hospitalICU = all_data_hospitalizations
  )
  return(rbindlist(d, idcol = "series"))
},cl=cl)
parallel::stopCluster(cl)

names(data_all_datasets) <- paste0("dataset_", 1:10)
data_all_datasets <- rbindlist(data_all_datasets, idcol = "dataset")



plot_ds3countryA_infs_rCases <- parallel::mcmapply(function(ds, ser, cc, title){
  g <- data_all_datasets[dataset == ds & series == ser & country == cc] %>% 
    ggplot(aes(x = date)) + 
    geom_ribbon(aes(ymin = q_l, ymax = q_u), fill = "steelblue3", alpha = 0.3)+
    geom_ribbon(aes(ymin = q_ll, ymax = q_uu), fill = "steelblue2", alpha = 0.2)+
    geom_line(aes(y= m), color = "royalblue4") + 
    xlab("date") + ylab(title)+
    ggtitle(paste("Posterior Predictions for", title, ":", ds, cc)) + theme_bw()
  
  if (ser == "cases"){
    return(g + geom_line(aes(y = C_t), alpha = 0.7, size = 0.3, alpha = 0.7))
  } else if (ser == "infections"){
    return(g + geom_line(aes(y = I_t), alpha = 0.7, size = 0.3, alpha = 0.7))
  } else if (ser == "deaths"){
    return(g + geom_line(aes(y = D_t), alpha = 0.7, size = 0.3, alpha = 0.7))
  } else if (ser == "repCases"){
    return(g + geom_line(aes(y = repC_t), alpha = 0.7, size = 0.3, alpha = 0.7))
  } else if (ser == "hospitalBeds"){
    return(g + geom_line(aes(y = H_t), alpha = 0.7, size = 0.3, alpha = 0.7))
  } else if(ser == "hospitalICU"){
    return(g + geom_line(aes(y = Hicu_t), alpha = 0.7, size = 0.3, alpha = 0.7))
  }
  
}, ds = c("dataset_3", "dataset_3"),
ser = c("infections", "repCases"),
cc = c("countryA", "countryA"),
title = c("infections", "reported cases"),
SIMPLIFY = F)

plot_ds3countryA_infs_rCases <- cowplot::plot_grid(plot_ds3countryA_infs_rCases[[1]],plot_ds3countryA_infs_rCases[[2]],
                                                      labels = c("A", "B"))

ggsave(paste0(PATH_PLOTS, "/posterior_predictions_countryA_ds3.pdf"),
       plot_ds3countryA_infs_rCases, width = 14, height = 5)

