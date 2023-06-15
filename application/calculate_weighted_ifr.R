############################################################################
#### calculates time dependent ifr which gets corrected for vaccinated people
############################################################################
library(magrittr)
library(data.table)
library(ggplot2)

# age specific ifrs from O'Driscol
ifrs <- c(0.002,0, 0.009, 0.122, 0.992, 7.274)/100

countries <- c("Austria", "Czechia","Germany", "Italy","France","Spain","Denmark",
               "Netherlands","Portugal", "United Kingdom", "Norway", "Switzerland",
               "Poland", "Sweden", "Finland", "Ireland",
               "Belgium", "Greece", "Hungary", "Slovenia")

# vaccination data from france
# d <- fread("https://www.data.gouv.fr/fr/datasets/r/54dd5f8d-1e2e-4ccb-8fb8-eac68245befd")
d <- fread("../data/data_vaccination/vacsi-a-fra-2023-02-14-19h00.csv")


d[, fra:= NULL][, couv_dose1 := NULL][,couv_complet := NULL]
d[, age_group := as.factor(clage_vacsi)]
d <- d[clage_vacsi != "0", ]
# d[, lapply(.SD, sum), by = c("jour"), .SDcols = c("n_dose1", "n_complet", "n_cum_dose1", "n_cum_complet")]
setnames(d, c("jour", "n_dose1", "n_complet", "n_cum_dose1", "n_cum_complet"), 
         c("date", "first_vacc", "second_vacc", "first_vacc_cum", "second_vacc_cum"))


vaccinations_france <- d[, .(date, age_group, first_vacc, second_vacc)]
# vaccinations_france <- d[, .(date, age_group, first_vacc_cum, second_vacc_cum)]
vaccinations_france[, sum_first_vacc := cumsum(first_vacc), by = age_group]
vaccinations_france[, sum_second_vacc := cumsum(second_vacc), by = age_group]

ggplot(vaccinations_france) + 
  geom_line(aes(date, sum_first_vacc, color = age_group)) + 
  geom_line(aes(date, sum_second_vacc, color = age_group)) 



# marginal vaccinations over all strata?
vaccinations_france[, sum_first_vacc_marginal := sum(sum_first_vacc), by = date]
vaccinations_france[, sum_second_vacc_marginal := sum(sum_second_vacc), by = date]

ggplot(vaccinations_france) + 
  geom_line(aes(date, sum_first_vacc_marginal))+
  geom_line(aes(date, sum_second_vacc_marginal), col = "darkred")



###############################
# adapt this to other countries
###############################

# calc proportions
vaccinations_france[, proportion_strat_first := sum_first_vacc/sum_first_vacc_marginal]
vaccinations_france[, proportion_strat_second := nafill(sum_second_vacc/sum_second_vacc_marginal, fill = 0)]


#### now use france to calculate vaccinations to other countries

# need nb of vaccinations in each country to each time point
# data_vaccs_all_countries <- fread("https://covid.ourworldindata.org/data/owid-covid-data.csv")
data_vaccs_all_countries <- fread("../data/data_vaccination/owid-covid-data.csv")
data_vaccs_all_countries <- data_vaccs_all_countries[location %in% countries & date >= "2020-12-27",
                                 .(date, location, people_vaccinated, people_fully_vaccinated)]

# check where first vaccination is zero
data_vaccs_all_countries[, lapply(.SD, function(x) return(x[1:3])), .SDcols = c("people_vaccinated", "people_fully_vaccinated", "date"), by = location] 

foo <- function(people_vaccinated){
  print(people_vaccinated[1])
  if (is.na(people_vaccinated[1])){
    people_vaccinated[1] <- 0
  }
  return(people_vaccinated)
}
data_vaccs_all_countries[, people_vaccinated := foo(people_vaccinated), by = location]
data_vaccs_all_countries[, people_fully_vaccinated := foo(people_fully_vaccinated), by = location]


data_vaccs_all_countries[, c("people_vaccinated","people_fully_vaccinated") := lapply(.SD, nafill, type = "locf"), .SDcols = c("people_vaccinated", "people_fully_vaccinated")]

ggplot(data_vaccs_all_countries) +
    geom_line(aes(date, people_vaccinated, color = location))



calc_vaccinations_country <- function(country, vaccinations_france, data_vaccs_all_countries){
  vaccs_single_country <- data_vaccs_all_countries[location == country]
  merged_data <- merge(vaccinations_france[, .(date, age_group, proportion_strat_first,proportion_strat_second)],
                       vaccs_single_country, by = "date", all.x = T)
  merged_data[, sum_first_vacc := proportion_strat_first * people_vaccinated]
  merged_data[, sum_second_vacc := proportion_strat_second * people_fully_vaccinated]
  
  return(merged_data[, .(date, age_group, sum_first_vacc, sum_second_vacc)])
}

vaccinations_age_specific <- lapply(countries, calc_vaccinations_country, 
       vaccinations_france=vaccinations_france, data_vaccs_all_countries=data_vaccs_all_countries) %>% 
  set_names(countries)


ggplot(rbindlist(vaccinations_age_specific, idcol = "country")) + geom_line(aes(date, sum_first_vacc, color = age_group))+
  geom_line(aes(date, sum_second_vacc, color = age_group), linetype = "dashed")+ 
  facet_wrap("country", scales = "free")



# redefine age groups
vaccinations_age_specific[[1]]$age_group %>% unique()


tmp <- vaccinations_age_specific[[1]]
tmp[, age_num := as.numeric(age_group)-1]
# tmp <- tmp[-nrow(tmp)]
foo <- function(z){
  if (z <= 3) return("0-34")
  if (z %in% c(4,5)) return("35-59")
  if (z %in% c(6,7,8,9)) return("60-79")
  if (z == 10) return("80+")
}

foo <- function(z){
  if (z <= 6) return("0-34")
  if (z %in% c(7,8,9)) return("35-59")
  if (z %in% c(10:13)) return("60-79")
  if (z == 14) return("80+")
}

tmp[, age_group_new := sapply(age_num, foo)]
tmp[, .(sum_first_vacc = sum(sum_first_vacc), sum_second_vacc = sum(sum_second_vacc)), by = c("date","age_group_new")]


ggplot(tmp[, sum(sum_first_vacc), by = c("date","age_group_new")]) + geom_line(aes(date, V1, color = age_group_new))



vaccinations_age_specific <- lapply(vaccinations_age_specific, function(tmp){
  tmp[, age_num := as.numeric(age_group)-1]
  # tmp <- tmp[-nrow(tmp)]
  foo <- function(z){
    if (z <= 3) return("0-34")
    if (z %in% c(4,5)) return("35-59")
    if (z %in% c(6,7,8,9)) return("60-79")
    if (z == 10) return("80+")
  }
  foo <- function(z){
    if (z <= 6) return("0-34")
    if (z %in% c(7,8,9)) return("35-59")
    if (z %in% c(10:13)) return("60-79")
    if (z == 14) return("80+")
  }
  tmp[, age_group_new := sapply(age_num, foo)]
  tmp[, .(sum_first_vacc = sum(sum_first_vacc), sum_second_vacc = sum(sum_second_vacc)), by = c("date","age_group_new")]
})

ggplot(vaccinations_age_specific[[1]]) + geom_line(aes(date, sum_first_vacc, color = age_group_new))
ggplot(vaccinations_age_specific[[3]]) + geom_line(aes(date, sum_first_vacc, color = age_group_new))
ggplot(vaccinations_age_specific[[5]]) + geom_line(aes(date, sum_first_vacc, color = age_group_new))
ggplot(vaccinations_age_specific[[6]]) + geom_line(aes(date, sum_first_vacc, color = age_group_new))



# check if every date is existing
date_seq <- seq.Date(from = as.Date("2020-12-27"), to = as.Date("2021-05-05"), by = "day")
lapply(vaccinations_age_specific, function(x){
  date_seq %in% x$date
}) %>% sapply(all) %>% all

# check if every age group is existing at each date
date_seq <- seq.Date(from = as.Date("2020-12-27"), to = as.Date("2021-05-05"), by = "day")
age_groups <- vaccinations_age_specific %>% rbindlist() %>% extract2("age_group_new") %>% unique()
lapply(vaccinations_age_specific, function(x){
  x[, length(age_groups %in% age_group_new) == 4, by = date]$V1 %>% all
}) %>% unlist %>% all



# assign ifr in each group
vaccinations_age_specific_new <- lapply(vaccinations_age_specific, function(tmp){
  ifrs <- data.table(
    age_group_new = c("0-34", "35-59", "60-79", "80+"),
    ifr_age = c(0.008, 0.122, 0.992, 7.274)/100
  )
  merge(tmp, ifrs, by = "age_group_new", all.x = T)
}) %>% rbindlist(idcol = "country")


# show vaccinations in each category in each country
ggplot(vaccinations_age_specific_new) + geom_line(aes(date, sum_first_vacc, color = age_group_new)) +
  facet_wrap("country", scales = "free")


# should be 4 horizontal lines, one for each age group
vaccinations_age_specific_new[country == "Germany"] %>% 
  ggplot()+
  geom_line(aes(date, ifr_age, color = age_group_new))



## add info about variant of concern to inflate ifr 


vocs <- fread("../data/variants_of_concern.csv")
vocs[country == "UnitedKingdom", country := "United Kingdom"]


# voc time series starts much earlier, we need all infos about alpha and delta
# therefore expand vaccinations_age_specific_new back to a date where alpha is 0
vocs[proportion_alpha != 0, date] %>% min
vaccinations_age_specific_new[, lapply(.SD, first), by = c("country", "age_group_new")]

vaccinations_age_specific_new[, date := as.Date(date)]

tmp_data_seq <- seq.Date(from = vocs[proportion_alpha != 0, date] %>% min,
                         to = vaccinations_age_specific_new$date %>% max,
                         by = "day")

vaccinations_age_specific_new <- vaccinations_age_specific_new %>% split(vaccinations_age_specific_new$country) %>% 
  lapply(function(df){
    split(df, df$age_group_new)
  }) %>% lapply(function(ll){
    lapply(ll, function(df){
      df_res <- merge(data.table(date=tmp_data_seq), df, all = T) %>% 
        tidyr::fill(., country, country, age_group_new, ifr_age, .direction = "up") %>% 
        tidyr::replace_na(list(sum_first_vacc=0, sum_second_vacc=0))
      return(df_res)
    }) %>% rbindlist()
  }) %>% rbindlist



vaccinations_age_specific_new <- merge(vaccinations_age_specific_new, vocs, all.x = T, by = c("country","date"))

# source for effect of alpha and delta: doi 10.1503/cmaj.211248;
vaccinations_age_specific_new[, ifr_age_voc_adjusted := ifr_age*proportion_original +
                     ifr_age*proportion_alpha*1.51 +
                     ifr_age*proportion_delta*2.08
    ]

# blocks with two lines indicate additional code for alternative IFRs used in the sensitivity analysis
# further IFRs for sensitivity analysis
# ______________________________________________________________________________
# ______________________________________________________________________________
vaccinations_age_specific_new_voc_low <- data.table::copy(vaccinations_age_specific_new)
vaccinations_age_specific_new_voc_high <- data.table::copy(vaccinations_age_specific_new)

vaccinations_age_specific_new_voc_low[, ifr_age_voc_adjusted := ifr_age*proportion_original +
                                ifr_age*proportion_alpha*(1 + 0.51*0.75)  +
                                ifr_age*proportion_delta*(1 + 1.08*1.25)
]
vaccinations_age_specific_new_voc_high[, ifr_age_voc_adjusted := ifr_age*proportion_original +
                                ifr_age*proportion_alpha*(1 + 0.51*1.25)  +
                                ifr_age*proportion_delta*(1 + 1.08*1.25)
]

# vaccinations_age_specific_new_voc_high[country == "Germany"] %>% ggplot(aes(date))+
#   geom_line(aes(y=ifr_age),color="black")+
#   geom_line(aes(y=ifr_age_voc_adjusted),color="blue")+
#   facet_wrap(~age_group_new, scales = "free")

# ______________________________________________________________________________
# ______________________________________________________________________________



vaccinations_age_specific_new[country == "Germany"] %>% ggplot(aes(date))+
  geom_line(aes(y=ifr_age),color="black")+
  geom_line(aes(y=ifr_age_voc_adjusted),color="blue")+
  facet_wrap(~age_group_new, scales = "free")



## get info about age distribution in each country
## source for population stems from oDriscoll

# function to define age groups
aggregate_groups <- function(d){
  d[, full := M+F]
  d[, age_group := 1:21]  
  foo <- function(x){
    if (x <= 7) return("0-34")
    if (x %in% c(8:12)) return("35-59")
    if (x %in% c(13:16)) return("60-79")
    if (x >16) return("80+")
  }
  d[, age_group := sapply(age_group, foo)]
  d[, sum(full), by = age_group]
}

# need specific vector:
# problems occur at uk (splitted in sub-locations) and czechia (name)
countries_tmp <- c(countries, c("England", "Northern%20Ireland", "Scotland", "Wales", "Czech%20Republic"))
countries_tmp <- countries_tmp[! countries_tmp == "Czechia"]
countries_tmp <- countries_tmp[! countries_tmp == "United Kingdom"]


populations <- lapply(countries_tmp, function(cc){
  print(paste("Downloading", cc))
  path <- paste0("https://raw.githubusercontent.com/meganodris/International-COVID-IFR/master/data/population/", cc, "-2019.csv")
  fread(path)
}) %>% lapply(aggregate_groups) %>% 
  set_names(countries_tmp) %>% 
  rbindlist(idcol = T) %>% 
  set_names(c("country", "age_group_new", "population_age_strata"))

populations[, population := sum(population_age_strata), by = country]

# sum UK
pop_UK <- populations[country %in% c("England", "Northern%20Ireland", "Scotland", "Wales"), lapply(.SD, sum), by = age_group_new, .SDcols = c("population_age_strata", "population")]
pop_UK[, country := "United Kingdom"]
populations <- rbindlist(list(populations, pop_UK), use.names = T)
populations <- populations[!(country %in% c("England", "Northern%20Ireland", "Scotland", "Wales"))]

# rename Czechia
populations[country == "Czech%20Republic", country := "Czechia"]
populations


##############################
# now reweight the ifr with the vaccinations in each age strata
##############################

# bring population information together with vaccination infos
vaccinations_age_specific_new <- merge(vaccinations_age_specific_new, populations, all.x = T, by = c("country", "age_group_new"))

# control for age population: cant get more people vaccinated than existing in its strata (problem can occur since we use approximation from france)
vaccinations_age_specific_new[sum_first_vacc > population_age_strata, sum_first_vacc := population_age_strata]

vaccinations_age_specific_new %>% ggplot() + geom_line(aes(date, sum_first_vacc, color = age_group_new))+
  facet_wrap("country", scales = "free")


# ______________________________________________________________________________
# ______________________________________________________________________________
vaccinations_age_specific_new_voc_low <- merge(vaccinations_age_specific_new_voc_low, populations, all.x = T, by = c("country", "age_group_new"))
vaccinations_age_specific_new_voc_low[sum_first_vacc > population_age_strata, sum_first_vacc := population_age_strata]

vaccinations_age_specific_new_voc_high <- merge(vaccinations_age_specific_new_voc_high, populations, all.x = T, by = c("country", "age_group_new"))
vaccinations_age_specific_new_voc_high[sum_first_vacc > population_age_strata, sum_first_vacc := population_age_strata]
# ______________________________________________________________________________
# ______________________________________________________________________________





#  calculate ifr, assume 80% safety against dying
# this line is without voc adjustment
vaccinations_age_specific_new[, ifr_t_m := ((population_age_strata-sum_first_vacc*0.8)*ifr_age)/population]
# with voc adjustment
vaccinations_age_specific_new[, ifr_t_m_voc_adjusted := ((population_age_strata-sum_first_vacc*0.8)*ifr_age_voc_adjusted)/population]



# ______________________________________________________________________________
# ______________________________________________________________________________
vaccinations_age_specific_new_voc_low[, ifr_t_m_voc_adjusted := ((population_age_strata-sum_first_vacc*0.8)*ifr_age_voc_adjusted)/population]
vaccinations_age_specific_new_voc_high[, ifr_t_m_voc_adjusted := ((population_age_strata-sum_first_vacc*0.8)*ifr_age_voc_adjusted)/population]

vaccinations_age_specific_new_vacc_low <- data.table::copy(vaccinations_age_specific_new)
vaccinations_age_specific_new_vacc_high <- data.table::copy(vaccinations_age_specific_new)
vaccinations_age_specific_new_vacc_low[, ifr_t_m_voc_adjusted := ((population_age_strata-sum_first_vacc*(1-0.2*0.75))*ifr_age_voc_adjusted)/population]
vaccinations_age_specific_new_vacc_high[, ifr_t_m_voc_adjusted := ((population_age_strata-sum_first_vacc*(1-0.2*1.25))*ifr_age_voc_adjusted)/population]
# ______________________________________________________________________________
# ______________________________________________________________________________


ggplot(vaccinations_age_specific_new)+
  geom_line(aes(date, ifr_age_voc_adjusted, color = age_group_new)) + # adjusted
  geom_line(aes(date, ifr_age, color = age_group_new))+ # without adjustment
  facet_wrap(~country)


# sum ifrtm together to get the true ifrtm
ifr_t_m_all_countries <- vaccinations_age_specific_new[, .(ifr_t_m = sum(ifr_t_m)), by = c("country", "date")]
ifr_t_m_all_countries_new <- vaccinations_age_specific_new[, .(ifr_t_m_voc_adjusted = sum(ifr_t_m_voc_adjusted)), by = c("country", "date")]

# ______________________________________________________________________________
# ______________________________________________________________________________
ifr_t_m_all_countries_new_voc_low <- vaccinations_age_specific_new_voc_low[, .(ifr_t_m_voc_adjusted = sum(ifr_t_m_voc_adjusted)), by = c("country", "date")]
ifr_t_m_all_countries_new_voc_high <- vaccinations_age_specific_new_voc_high[, .(ifr_t_m_voc_adjusted = sum(ifr_t_m_voc_adjusted)), by = c("country", "date")]

ifr_t_m_all_countries_new_vacc_low <- vaccinations_age_specific_new_vacc_low[, .(ifr_t_m_voc_adjusted = sum(ifr_t_m_voc_adjusted)), by = c("country", "date")]
ifr_t_m_all_countries_new_vacc_high <- vaccinations_age_specific_new_vacc_high[, .(ifr_t_m_voc_adjusted = sum(ifr_t_m_voc_adjusted)), by = c("country", "date")]
# ______________________________________________________________________________
# ______________________________________________________________________________



g1 <- ggplot(ifr_t_m_all_countries) + geom_line(aes(date,ifr_t_m,color = country))  
g2 <- ggplot(ifr_t_m_all_countries_new) + geom_line(aes(date,ifr_t_m_voc_adjusted,color = country))  
g3 <- ggplot(merge(ifr_t_m_all_countries_new, ifr_t_m_all_countries)) + geom_line(aes(date,ifr_t_m_voc_adjusted-ifr_t_m,color = country))  
egg::ggarrange(g1,g2,g3,ncol = 1)
# plotly::ggplotly(g2)


# expand this back to start of pandemic 
ifr_t_m_all_countries
ifr_t_m_all_countries_new
ifr_t_m_all_countries_long <- lapply(unique(ifr_t_m_all_countries$country), function(cc){
  tmp <- data.table(country = cc, date = seq.Date(as.Date("2020-01-01"), ifr_t_m_all_countries$date %>% max, by = "day"))
  tmp <- merge(tmp, ifr_t_m_all_countries[country == cc], all.x = T)
  tmp[, ifr_t_m := nafill(ifr_t_m, type = "nocb")]
}) %>% rbindlist()

ifr_t_m_all_countries_long_new <- lapply(unique(ifr_t_m_all_countries_new$country), function(cc){
  tmp <- data.table(country = cc, date = seq.Date(as.Date("2020-01-01"), ifr_t_m_all_countries_new$date %>% max, by = "day"))
  tmp <- merge(tmp, ifr_t_m_all_countries_new[country == cc], all.x = T)
  tmp[, ifr_t_m_voc_adjusted := nafill(ifr_t_m_voc_adjusted, type = "nocb")]
}) %>% rbindlist()


# ______________________________________________________________________________
# ______________________________________________________________________________
ifr_t_m_all_countries_long_new_voc_low <- lapply(unique(ifr_t_m_all_countries_new_voc_low$country), function(cc){
  tmp <- data.table(country = cc, date = seq.Date(as.Date("2020-01-01"), ifr_t_m_all_countries_new_voc_low$date %>% max, by = "day"))
  tmp <- merge(tmp, ifr_t_m_all_countries_new_voc_low[country == cc], all.x = T)
  tmp[, ifr_t_m_voc_adjusted := nafill(ifr_t_m_voc_adjusted, type = "nocb")]
}) %>% rbindlist()

ifr_t_m_all_countries_long_new_voc_high <- lapply(unique(ifr_t_m_all_countries_new_voc_high$country), function(cc){
  tmp <- data.table(country = cc, date = seq.Date(as.Date("2020-01-01"), ifr_t_m_all_countries_new_voc_high$date %>% max, by = "day"))
  tmp <- merge(tmp, ifr_t_m_all_countries_new_voc_high[country == cc], all.x = T)
  tmp[, ifr_t_m_voc_adjusted := nafill(ifr_t_m_voc_adjusted, type = "nocb")]
}) %>% rbindlist()

ifr_t_m_all_countries_long_new_vacc_low <- lapply(unique(ifr_t_m_all_countries_new_vacc_low$country), function(cc){
  tmp <- data.table(country = cc, date = seq.Date(as.Date("2020-01-01"), ifr_t_m_all_countries_new_vacc_low$date %>% max, by = "day"))
  tmp <- merge(tmp, ifr_t_m_all_countries_new_vacc_low[country == cc], all.x = T)
  tmp[, ifr_t_m_voc_adjusted := nafill(ifr_t_m_voc_adjusted, type = "nocb")]
}) %>% rbindlist()

ifr_t_m_all_countries_long_new_vacc_high <- lapply(unique(ifr_t_m_all_countries_new_vacc_high$country), function(cc){
  tmp <- data.table(country = cc, date = seq.Date(as.Date("2020-01-01"), ifr_t_m_all_countries_new_vacc_high$date %>% max, by = "day"))
  tmp <- merge(tmp, ifr_t_m_all_countries_new_vacc_high[country == cc], all.x = T)
  tmp[, ifr_t_m_voc_adjusted := nafill(ifr_t_m_voc_adjusted, type = "nocb")]
}) %>% rbindlist()

# ______________________________________________________________________________
# ______________________________________________________________________________



# delay ifr with 14 days
ifr_t_m_all_countries_long[, ifr_t_m := nafill(shift(ifr_t_m, 14), "nocb"), by = country]
ifr_t_m_all_countries_long_new[, ifr_t_m := nafill(shift(ifr_t_m_voc_adjusted, 14), "nocb"), by = country]

ggplot(ifr_t_m_all_countries_long) + geom_line(aes(date, ifr_t_m, color = country))
ggplot(ifr_t_m_all_countries_long_new) + geom_line(aes(date, ifr_t_m, color = country))

# ______________________________________________________________________________
# ______________________________________________________________________________
ifr_t_m_all_countries_long_new_voc_low[, ifr_t_m := nafill(shift(ifr_t_m_voc_adjusted, 14), "nocb"), by = country]
ifr_t_m_all_countries_long_new_voc_high[, ifr_t_m := nafill(shift(ifr_t_m_voc_adjusted, 14), "nocb"), by = country]
ifr_t_m_all_countries_long_new_vacc_low[, ifr_t_m := nafill(shift(ifr_t_m_voc_adjusted, 14), "nocb"), by = country]
ifr_t_m_all_countries_long_new_vacc_high[, ifr_t_m := nafill(shift(ifr_t_m_voc_adjusted, 14), "nocb"), by = country]
g1 <- ggplot(ifr_t_m_all_countries_long_new_voc_low) + geom_line(aes(date, ifr_t_m, color = country)) + ylim(0, 0.013) + ggtitle("voc_low")
g2 <- ggplot(ifr_t_m_all_countries_long_new_voc_high) + geom_line(aes(date, ifr_t_m, color = country)) + ylim(0, 0.013) + ggtitle("voc_high")
g3 <- ggplot(ifr_t_m_all_countries_long_new_vacc_low) + geom_line(aes(date, ifr_t_m, color = country)) + ylim(0, 0.013) + ggtitle("vacc_low")
g4 <- ggplot(ifr_t_m_all_countries_long_new_vacc_high) + geom_line(aes(date, ifr_t_m, color = country)) + ylim(0, 0.013) + ggtitle("vacc_high")

egg::ggarrange(g1,g2,g3,g4)
# ______________________________________________________________________________
# ______________________________________________________________________________


fwrite(ifr_t_m_all_countries_long_new_voc_low, "../data/ifr/ifr_t_m_shifted_voc_low.csv")
fwrite(ifr_t_m_all_countries_long_new_voc_high, "../data/ifr/ifr_t_m_shifted_voc_high.csv")
fwrite(ifr_t_m_all_countries_long_new_vacc_low, "../data/ifr/ifr_t_m_shifted_vacc_low.csv")
fwrite(ifr_t_m_all_countries_long_new_vacc_high, "../data/ifr/ifr_t_m_shifted_vacc_high.csv")



####### 
# write data
#fwrite(ifr_t_m_all_countries_long, "../data/ifr_t_m.csv")
fwrite(ifr_t_m_all_countries_long_new, "../data/ifr/ifr_t_m_shifted.csv")

