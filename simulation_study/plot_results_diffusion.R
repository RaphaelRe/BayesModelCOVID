# this script is roughly similar to the standard evaluation, but for the diffusion model
library(data.table)
library(magrittr)
library(ggplot2)
library(egg)
library(scales)
library(cowplot)

CORES <- 50

PATH_RESULTS = "results/diffusion/"
PATH_DATA_PART = "../data/simulated_data/diffusion/"
PATH_PLOTS = "plots/diffusion/"

if (!dir.exists(PATH_PLOTS)) dir.create(PATH_PLOTS, recursive = T)


N_datasets <- list.files(PATH_RESULTS) %>% length


# get all interventions
INTERVENTIONS = list.files(paste0(PATH_RESULTS, "res_dataset_1/")) %>% grep("alpha",.) %>% 
  (function(x)(list.files(paste0(PATH_RESULTS, "res_dataset_1/"))[x])) %>% 
  strsplit("_") %>% sapply(function(x) x[4]) %>% 
  unique() %>% setdiff("alpha.txt")

# get all countries
ccs = fread(paste0(PATH_DATA_PART, "/data_sim_5NPIs_1_diffusion.csv"))$country %>% unique()



# reads for a parameter all chains apply thinning and returns a meldted data frame
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
  #data = melt(data, id.vars = c("i", "parameter")) %>% set_names(c("i", "parameter","chain", "value"))
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
  true_vals_ds <- fread(paste0(PATH_DATA_PART, "/data_sim_5NPIs_",ds, "_diffusion.csv"))[, c("country", paste0("true_val_NPI", 1:5))
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
  d <- fread(paste0(PATH_DATA_PART, "/data_sim_5NPIs_",ds, "_diffusion.csv"))
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

