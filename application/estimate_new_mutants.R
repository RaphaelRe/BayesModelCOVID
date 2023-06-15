################################################################################
### please note that due to changes in the data bases, the script does not work anymore
################################################################################

library(data.table)
library(ggplot2)
library(magrittr)

# weekly data on vaccinations
d <- fread("https://opendata.ecdc.europa.eu/covid19/virusvariant/csv/data.csv")

# add date: Assume always mid of the week
todate <- function(yw){
  a <- yw
  stringi::stri_sub(a, 6,5) <- "W"
  ISOweek::ISOweek2date(paste0(a, "-4"))
}

d[, date := todate(year_week)]
d$date


# plotly::ggplotly(
#   ggplot(d[country == "Germany"]) + 
#     geom_line(aes(date, percent_variant, color = variant))
# )

ggplot(d[country %in% c("Germany", "France", "Italy")]) + 
  geom_line(aes(date, percent_variant, color = variant)) + facet_wrap("country")




# for spain something is wrong in the data, set first 2 non zero points to zero
d[country == "Spain" & source == "GISAID" & variant %in% c("B.1.1.7", "B.1.617.2")] %>% 
  ggplot(aes(date, percent_variant, color = variant))+geom_point()
d[country == "Spain" & source == "GISAID" & variant == "B.1.1.7" & date < "2020-09-01", percent_variant := 0]





estimate_proportion <- function(country_sel, d, variant_sel = "B.1.1.7", plot = F){
  if (!(country_sel %in% d$country)) stop("Country not in dataset")
  y <- d[country == (country_sel) & variant ==  (variant_sel) & source == "GISAID", .(date,percent_variant, source)]
  if (country_sel == "Hungary"){
    # special case Hungary - only TESSy available
    y <- d[country == (country_sel) & variant ==  "B.1.617.2" & source == "TESSy", .(date,percent_variant, source)]
  }
  
  # get date of one week before first non zero value
  start_date <- y$date[1]
  start_date <- y[percent_variant > 0 ]$date[1]-7
  
  # let it start earlier (pad with zeros)
  y <- y[date >= (start_date)]$percent_variant %>% c(0, 0,0,.) %>% data.table::nafill(type = "locf")
  
  y <- y[1:which(y == max(y, na.rm = T))] /100
  x <- 0:(length(y)-1)
 
  # sigmoid function
  foo <- function(x, b = 0.4,c = 10){
    1/(1 + exp(-b * (x-c)))
  }
 
  # squared loss
  loss <- function(pars){
    sum(abs((y - foo(x, pars[1],pars[2])))^2)
  }
 
  # optimizer function
  optim_sigmoid <- function(pars){
    opt <- optim(pars, loss, method = "BFGS")
    return(list(par = opt$par, loss = opt$value))
  }
  
  
  n = 10
  b = runif(n, 0, 3)
  c = runif(n, 2, 30)
  pos = cbind(b,c)
  print(paste0("Start optimizer for: ", country_sel))
  cl <- parallel::makeForkCluster(7)
  results <- pbapply::pbapply(pos, 1,optim_sigmoid, cl = cl)
  parallel::stopCluster(cl)
  (pars <- sapply(results, function(l) l$loss) %>% which.min %>% results[.] %>% extract2(1) %>% extract2(1))
  
  days <- seq(0, 100, length.out = 100 * 7)
  pred <- foo(days, pars[1], pars[2])
  
  if (plot){
    # foo(x, pars[1], pars[2]) %>% plot(pch = 20, main = country_sel)
    # points(x,y, col=2, pch = 20)
    plot(x, y, pch = 20, xlim = c(0, 110), ylim = c(0,1), main = country_sel)
    points(days, pred, pch = 20, col = 2)
  }
  # check if first value is close enough to zero
  #if (pred[1] > 0.05) stop("first value should be smaller 5%")
  if (pred[1] > 0.05) warning(paste0("first value for ", country_sel, " was ", pred[1], "Is rather high, should be near to 0"))
  return(list(pred=pred ,start = start_date-21, pars = pars)) # add 2 week earlier start
}

ccs <- c("Austria", "Czechia","Germany", "Italy","France","Spain","Denmark",
               "Netherlands","Portugal", "United Kingdom", "Norway", "Switzerland",
               "Poland", "Sweden", "Finland", "Ireland",
         "Belgium", "Greece", "Hungary", "Slovenia")



# alpha aka B.1.1.7
voc_alpha <- sapply(ccs, function(cc, d) tryCatch(estimate_proportion(cc, d), error = function(e) e), d = d, simplify = F, USE.NAMES = T)


# add UK (other data source)
UK_alpha <- data.table(country = "United Kingdom",
                 variant = c(rep("B.1.1.7", 17)),
                 percent_variant = c(0, 0.05, 0.35, 1,2.7,6.3,10.2,10.2,13.9,
                                       32.9,45.7,51.3,70.6,74.1,78.6,86,89.1),
                 date = seq.Date(as.Date("2020-10-15"), len = 17, by = "week"),
                 source = "GISAID" # not correct but the necessary to use the function above
                 ) %>% 
  merge(data.table(date = seq.Date(as.Date("2020-10-01"), d$date %>% max, by = "week")),., all=T)

UK_delta <- data.table(country = "United Kingdom",
                 variant = c(rep("B.1.617.2", 22)),
                 percent_variant = c(0.0833,0.1428,0.0084,0,0,0,0,0,0,0,0.0014,
                                     0,0,0.0074,0.0867,0.3950,0.6287,0.8394,0.9163,
                                     0.9442,0.9668,0.9688)*100,
                 date = as.Date(c("31.01.2021","07.02.2021","14.02.2021","21.02.2021",
                                  "28.02.2021","07.03.2021","14.03.2021", "21.03.2021",
                                  "28.03.2021","04.04.2021","11.04.2021","18.04.2021",
                                  "25.04.2021","02.05.2021", "09.05.2021","16.05.2021",
                                  "23.05.2021","30.05.2021","06.06.2021","13.06.2021",
                                  "20.06.2021","27.06.2021"), format = "%d.%m.%Y")-3,# to get it on mif of the week
                 source = "GISAID" # not correct but the necessary to use the function above
                 ) %>% 
  merge(data.table(date = seq.Date(as.Date("2020-10-01"), d$date %>% max, by = "week")),., all=T)




# add Switzerland (other data source https://www.covid19.admin.ch/api/data/20210806-h6laxk40/sources/COVID19Variants_wgs.csv)

d_CH <- fread("https://www.covid19.admin.ch/api/data/20210806-h6laxk40/sources/COVID19Variants_wgs.csv")

d_CH_alpha <- d_CH[which(!is.na(d_CH$freq) & (d_CH$variant_type == "B.1.1.7")),]

thursdays = seq.Date(as.Date("2020-11-19"), len = 30, by = "week")
percent <- c()
for(day in thursdays){
  percent <- c(percent, colMeans(d_CH_alpha[c((which(d_CH_alpha$date==day)-3):(which(d_CH_alpha$date==day)+3)), "prct"]))
}

CH_alpha <- data.table(country = "Switzerland",
                       variant = c(rep("B.1.1.7", 30)),
                       percent_variant = percent,
                       date = thursdays,
                       source = "GISAID" # not correct but the necessary to use the function above
) %>% 
  merge(data.table(date = seq.Date(as.Date("2020-10-01"), d$date %>% max, by = "week")),., all=T)


d_CH_delta <- d_CH[which(!is.na(d_CH$freq) & (d_CH$variant_type == "B.1.617.2")),]

thursdays = seq.Date(as.Date("2021-04-01"), len = 17, by = "week")
percent <- c()
for(day in thursdays){
  percent <- c(percent, colMeans(d_CH_delta[c((which(d_CH_delta$date==day)-3):(which(d_CH_delta$date==day)+3)), "prct"]))
}

CH_delta <- data.table(country = "Switzerland",
                       variant = c(rep("B.1.617.2", 17)),
                       percent_variant = percent,
                       date = thursdays,
                       source = "GISAID" # not correct but the necessary to use the function above
) %>% 
  merge(data.table(date = seq.Date(as.Date("2020-10-01"), d$date %>% max, by = "week")),., all=T)



voc_alpha[["United Kingdom"]] <- estimate_proportion("United Kingdom", d = UK_alpha)
voc_alpha[["Switzerland"]] <- estimate_proportion("Switzerland", d = CH_alpha)

#ccs <- setdiff(ccs, c("Switzerland"))



# delta aka B.1.617.2

voc_delta <- sapply(ccs, function(cc, d) tryCatch(estimate_proportion(cc, d, variant_sel = "B.1.617.2"), error = function(e) e), d = d, simplify = F, USE.NAMES = T)

voc_delta[["United Kingdom"]] <- estimate_proportion("United Kingdom", d = UK_delta,variant_sel = "B.1.617.2")
voc_delta[["Switzerland"]] <- estimate_proportion("Switzerland", d = CH_delta,variant_sel = "B.1.617.2")




# plot this for the countries
plot_transition <- function(country_sel, voc_alpha, d, variant_sel = "B.1.1.7"){
  print(country_sel)
  d_sel <- d[country == country_sel & variant == (variant_sel) & source == "GISAID",.(date, percent_variant)]
  
  start <- voc_alpha[[country_sel]]$start
  end <- start + 300
  data <- data.table(date = seq.Date(start, end, by = "day")) %>% merge(d_sel, all.x = T)
  est <- voc_alpha[[country_sel]]$pred
  data[, estimate := est[1:nrow(data)]]
  ggplot(data, aes(x=date)) + 
    geom_point(aes(y=percent_variant/100), color = "steelblue1")+
    geom_line(aes(y=estimate),color = "steelblue4")+ ggtitle(paste(variant_sel, country_sel))
}

# lapply(ccs, plot_transition, voc_alpha=voc_delta, d=d, variant_sel = "B.1.617.2")

(g1 <- plot_transition("Germany", voc_alpha, d))
(g2 <- plot_transition("Germany", voc_delta, d, variant_sel = "B.1.617.2"))

ggsave("transition_alpha_ger.pdf", g1, height = 5, width = 4)
ggsave("transition_delta_ger.pdf", g2, height = 5, width = 4)




# get data together
combine_data <- function(country_sel, d1 = voc_alpha, d2 = voc_delta){
  print(country_sel)
  d1_sel <- d1[[country_sel]]
  d2_sel <- d2[[country_sel]]
  
 
  alpha_ts <- data.table(
    date = seq.Date(d1_sel$start, length.out = length(d1_sel$pred), by = "day"),
    proportion_alpha = d1[[country_sel]]$pred
    )
  delta_ts <- data.table(
    date = seq.Date(d2_sel$start, length.out = length(d2_sel$pred), by = "day"),
    proportion_delta = d2[[country_sel]]$pred
    )
  
  data.table(date = seq.Date(as.Date("2020-01-01"), 
                             as.Date("2021-08-31"), 
                             by = "day")) %>% 
    merge(., alpha_ts, all = T) %>%
    merge(., delta_ts, all = T) %>% return()
}

#melt(combine_data("Germany"), id.vars = "date") %>% ggplot(aes(date, value, color = variable)) + geom_line()

data_final <- sapply(ccs, combine_data, simplify = F, USE.NAMES = T) %>% rbindlist(idcol = "country")


data_final[, proportion_alpha := nafill(proportion_alpha, "locf"), by = country]
data_final[, proportion_delta := nafill(proportion_delta, "locf"), by = country]

data_final[, proportion_alpha := nafill(proportion_alpha, fill=0)]
data_final[, proportion_delta := nafill(proportion_delta, fill=0)]


data_final[, proportion_alpha := ifelse(proportion_alpha < 1-proportion_delta, proportion_alpha, 1-proportion_delta)]
data_final[, proportion_original := 1-proportion_alpha-proportion_delta]

ggplot(data_final) + geom_line(aes(date, proportion_alpha),color = "darkblue")+
  geom_line(aes(date, proportion_delta),color = "lightblue")+
  geom_line(aes(date, proportion_original),color = "orange")+
  facet_wrap("country", nrow = 3)

(g <- ggplot(data_final[country == "Germany"]) + geom_line(aes(date, proportion_alpha),color = "darkblue")+
  geom_line(aes(date, proportion_delta),color = "lightblue")+
  geom_line(aes(date, proportion_original),color = "orange")+
  xlim(as.Date("2020-08-01"), as.Date("2022-01-01"))+
  ggtitle("Combination of variants Germany"))

ggsave("transition_VOCs_ger.pdf", g, height = 5, width = 9)



data_final[country == "United Kingdom", country := "UnitedKingdom"]

# check if there are some NAs 
data_final %>% sapply(function(x) sum(is.na(x)))


# (path <- paste0("data/new_mutants_",
#                Sys.Date() %>% stringr::str_remove_all("-"), ".csv"))
# fwrite(data_final, path)
(path <- "../data/variants_of_concern.csv")
fwrite(data_final, path)


