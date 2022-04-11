################################################################################
# This script prepares the data to fit the model
################################################################################

library(dplyr)
library(magrittr)
library(ggplot2)
library(data.table)
library(lubridate)


PATH = getwd()
setwd(PATH)

################
# get basic data

### for cases
d <- fread("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv") %>% 
  filter(`Province/State` == "") %>% 
  .[,c("Province/State", "Lat", "Long") := NULL] %>% 
  melt("Country/Region") %>% 
  .[, variable := as.Date(variable, format = "%m/%d/%y")]

setnames(d, c("Country/Region", "variable", "value"),
         c("country", "date", "cum_cases"))
setorder(d, country, date)

# select countries
countries <- c("Austria", "Belgium", "Czechia", "Denmark", "Finland", "France",
               "Germany", "Greece", "Hungary", "Ireland", "Italy", "Netherlands",  
               "Norway", "Poland", "Portugal", "Slovenia", "Spain", "Sweden",
               "Switzerland", "United Kingdom")

d <- d[country %in% countries]

# NAs in data?
d$cases %>% is.na() %>% sum

# plot
ggplot(d) + geom_line(aes(date, cum_cases, color = country))

# add differences
d[, cases := c(NA, diff(cum_cases)), by = country]
ggplot(d) + geom_line(aes(date, cases, color = country))

# set negative values to 0 firstly
d[cases <= 0, cases := 0]


### now for deaths
dd <- fread("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv") %>% 
  filter(`Province/State` == "") %>% 
  .[,c("Province/State", "Lat", "Long") := NULL] %>% 
  melt("Country/Region") %>% 
  .[, variable := as.Date(variable, format = "%m/%d/%y")]

setnames(dd, c("Country/Region", "variable", "value"),
         c("country", "date", "cum_deaths"))
setorder(dd, country, date)


dd <- dd[country %in% countries]

# plot
ggplot(dd) + geom_line(aes(date, cum_deaths, color = country))

# add differences
dd[, deaths := c(NA, diff(cum_deaths)), by = country]
ggplot(dd) + geom_line(aes(date, deaths, color = country))


# remove 0
dd[deaths <= 0, deaths := 0]
ggplot(dd) +geom_line(aes(x = date, y = deaths, col = country))


####  bring cases and deaths together
data <- merge(d,dd, all = T)

ggplot(data) + 
  geom_line(aes(date, deaths, col = country), alpha = 0.6, linetype = "dashed") + 
  geom_line(aes(date, cases, col = country), alpha = 0.6)

data[, c("cum_cases", "cum_deaths") := NULL]
setnames(data, c("cases", "deaths"), c("repC_t", "D_t"))

# set name of UK to camelcase (algorithm needs this)
data[country == "United Kingdom", country := "UnitedKingdom"]



################################################
# define period changes of rho (under reporting)
periods <- c("2020-03-15", "2020-03-26", "2020-04-22", "2020-05-14",
             "2020-06-09", "2020-06-28", "2020-08-08", "2020-08-27", "2020-09-22",
             "2020-10-07", "2020-11-03", "2020-12-02", "2020-12-16", "2021-02-03",
             "2021-03-01", "2021-04-01", "2021-04-21", "2021-06-15",
             "2021-07-15", "2021-08-15", "2021-09-15", "2021-10-15","2021-11-15")


periods <- replicate(length(countries), periods, simplify = F) %>% 
  set_names(unique(data$country))


# function to generate periods from the cutoffs
calc_rho_period <- function(country_sel, periods_list, data){
  cutpoints <- periods_list[[country_sel]] %>% as.Date %>% sort
  cutpoints <- c(cutpoints, max(data[country == country_sel, date])+1)
  nb_periods <- length(cutpoints)
  rho_period <- data.table(date = data[country == country_sel, date], rho_period = -999) # last date is -999, but will be deleted anyway
  for (i in nb_periods:1) {
    rho_period[date < cutpoints[i], rho_period:= i]
  }
  return(rho_period)
}

rho_period <- sapply(names(periods), calc_rho_period, periods_list = periods, 
                     data = data, simplify = F) %>% 
  rbindlist(idcol = T) %>% 
  set_names(c("country", "date", "rho_period"))


data <- merge(data, rho_period, all.x = T)


#################################
# add weekdays (start of the week is encoded on monday with 1 and runs up to 7)
data[, weekday := lubridate::wday(date, week_start = T)]


##########
## Add VOC
data_voc <- fread("../data/variants_of_concern.csv")
setnames(data_voc, c("proportion_alpha", "proportion_delta"), c("alpha", "delta"))

data <- merge(data, data_voc, all.x = T)



##################
# Add vaccinations
vaccs <- fread("https://covid.ourworldindata.org/data/owid-covid-data.csv") %>%
  select(c("location","date","population", "people_vaccinated", "people_fully_vaccinated"))

vaccs[,date := as.Date(date)]

vaccs <- vaccs[location %in% countries]
vaccs[, first_vaccination := nafill(people_vaccinated/population, "locf"), by = location]
vaccs[, second_vaccination := nafill(people_fully_vaccinated/population, "locf"), by = location]
# viz
melt(vaccs, id.vars = c("location", "date"), measure.vars = c("first_vaccination", "second_vaccination")) %>% 
  ggplot(.) + geom_line(aes(date, value, color = location))+facet_grid(rows="variable")

setnames(vaccs, "location", "country")
vaccs[country == "United Kingdom", country:= "UnitedKingdom"]

# plot vaccination for one country
# vaccs_country <- vaccs[country == "Germany"]
# cols <- c("first_vaccination", "second_vaccination")
# vaccs_country[, (cols) := lapply(.SD, nafill, fill=0), .SDcols = cols]
# melt(vaccs_country, id.vars = c("country", "date"), measure.vars = c("first_vaccination", "second_vaccination")) %>%
#   ggplot(.) + geom_line(aes(date, value, color = variable))+scale_color_manual(values=c("steelblue4", "steelblue2"))+
#   theme(legend.position = "bottom", legend.title = element_blank())+ylab("vaccinated in %")
# ggsave("~/plot_vaccinations_ger.pdf", width = 7, height = 4)


# delay vaccinations 14 days
vaccs[, first_vaccination := shift(first_vaccination, 14), by = country]
vaccs[, second_vaccination := shift(second_vaccination, 14), by = country]

data <- merge(data, vaccs[, .(country,date, first_vaccination, second_vaccination)], all = T, by = c("country", "date"))
data[, first_vaccination := nafill(first_vaccination, fill = 0)]
data[, second_vaccination := nafill(second_vaccination, fill = 0)]



######################
# Add Hospitalizations

hh <- fread("https://opendata.ecdc.europa.eu/covid19/hospitalicuadmissionrates/csv/data.csv")
hh <- hh[country %in% countries][indicator %in% c("Daily hospital occupancy", "Daily ICU occupancy")][, .(country, indicator, date, value)]

ggplot(hh) + geom_line(aes(date, value, color=indicator))+facet_wrap(~country, scales = "free")

hh[, indicator := ifelse(indicator == "Daily hospital occupancy", "H_t", "Hicu_t")]

hh_wards = hh[indicator == "H_t"]
setnames(hh_wards, "value", "H_t")

hh_icu = hh[indicator == "Hicu_t"]
setnames(hh_icu, "value", "Hicu_t")

hh <- merge(hh_wards[, -c("indicator")], hh_icu[, -c("indicator")], by = c("country", "date"), all = T) 


# data from UK from other source
h_UK <- fread("https://api.coronavirus.data.gov.uk/v2/data?areaType=overview&metric=covidOccupiedMVBeds&metric=hospitalCases&format=csv")[, .(date, covidOccupiedMVBeds, hospitalCases)]

setnames(h_UK, c("covidOccupiedMVBeds", "hospitalCases"), c("Hicu_t", "H_t"))
h_UK[, country := "UnitedKingdom"]
# ggplot(h_UK) + 
#   geom_line(aes(date, Hicu_t), color = "darkblue")+
#   geom_line(aes(date, H_t), color = "steelblue3")
hh <- rbindlist(list(hh, h_UK), use.names = T)

# Switzerland: not yet used
# h_CH <- fread("https://www.covid19.admin.ch/api/data/20211006-cnzsfrjp/downloads/sources-csv.zip")

ggplot(hh) + 
  geom_line(aes(date, H_t),color = "blue")+
  geom_line(aes(date, Hicu_t), color = "red")+
  facet_wrap(~country, scales = "free")


data <- merge(data, hh, all = T)



# function to show all the data with 2 axis and rescaling option
plot_comparison <- function(cc, sD=20, sH=1){
    ggplot(data[country == cc]) + 
    geom_line(aes(date, repC_t, color = "repC_t")) + 
    geom_line(aes(date, D_t*sD, color="D_t"), alpha = 0.9) + 
    geom_line(aes(date, H_t*sH, color= "H_t"),alpha = 0.7) + 
    scale_y_continuous(sec.axis = sec_axis(~./sD, name = "deaths")) + 
    ggplot2::ggtitle(cc)+
    
    scale_colour_manual(values = c("steelblue4","steelblue","steelblue1"),name = "",
                        labels = c(bquote(paste(deaths%*%20)),
                                   bquote(paste("hospital ", occupancy)),
                                bquote(paste("reported ", cases)))
                        )+
    labs(y = "reported cases / hospital occupancy",  x = "Date",  colour = "Series")+
    theme(legend.position = 'bottom')+
    theme(axis.text.y.left = element_text(color = "steelblue1"))+
    theme(axis.title.y.left = element_text(color = "steelblue1"))+
    theme(axis.ticks.y.left = element_line(color = "steelblue1"))+
    theme(axis.text.y.right = element_text(color = "steelblue4"))+
    theme(axis.title.y.right = element_text(color = "steelblue4"))+
    theme(axis.ticks.y.right = element_line(color = "steelblue4")) #+ facet_wrap(facets = "country", scales = "free_y", ncol = 1)
}



plot_comparison("Spain")

# save plots
# ccs <- data$country %>% unique()
# pdf("~/Desktop/comparison_all_countries.pdf", width = 7, height = 4)
# lapply(ccs, plot_comparison)
# dev.off()



################################
# cut down to end of October
data <- data[date <= "2021-10-31"]


##################################
# delete outlier

# thresholds
# 2.5 
# 300 for cases
# 30 for deaths

# Specific case: peak in Spain v 2020-06-18
data[country == "Spain" & date == "2020-06-19", D_t := 0]

# smooth all other data points which meet the criteria with 25 50 25 %
data[, smooth_cases := as.vector(stats::filter(repC_t, rep(1/7, 7))), by = country]
data[, smooth_deaths := as.vector(stats::filter(D_t, rep(1/7, 7))), by = country]

data[, abs_deviation_cases := abs(repC_t - smooth_cases)]
data[, abs_deviation_deaths := abs(D_t - smooth_deaths)]
data[, rel_deviation_cases := abs(repC_t - smooth_cases)/smooth_cases]
data[, rel_deviation_deaths := abs(D_t - smooth_deaths)/smooth_deaths]

### cases
data[rel_deviation_cases > 2.5 & abs_deviation_cases > 300]
# if NaN: 0/0 - in this case rel deviation is 1
data[is.nan(rel_deviation_cases), rel_deviation_cases := 1]
data[is.nan(rel_deviation_deaths), rel_deviation_cases := 1]

ii <- (data$rel_deviation_cases > 2.5 & data$abs_deviation_cases > 300) %>% 
  as.numeric() %>% 
  nafill(., fill = 0) %>% as.logical()
all(data[ii] == data[rel_deviation_cases > 2.5 & abs_deviation_cases > 300], na.rm = T)

ii <- which(ii)


# function to reduce extreme outlieres
# # gives 50% of its value to its neighbors (i.e. 25% each)
smooth_value_cases <- function(index, d = data){
  value <- d[index, repC_t]
  value_i <- round(value * 0.5)
  value_neigh <- round(value * 0.25)
  d[index, repC_t := value_i]
  d[index-1, repC_t:= repC_t+value_neigh]
  d[index+1, repC_t:= repC_t+value_neigh]
  return(d)
}

for (i in ii) {
  data <- smooth_value_cases(i)
}


### deaths
data[rel_deviation_deaths > 2.5 & abs_deviation_deaths > 30]
ii <- (data$rel_deviation_deaths > 2.5 & data$abs_deviation_deaths > 30) %>% 
  as.numeric() %>% 
  nafill(., fill = 0) %>% as.logical()
all(data[ii] == data[rel_deviation_deaths > 2.5 & abs_deviation_deaths > 30], na.rm = T)

ii <- which(ii)

smooth_value_deaths <- function(index,d = data){
  value <- d[index, D_t]
  value_i <- round(value * 0.5)
  value_neigh <- round(value * 0.25)
  d[index, D_t := value_i]
  d[index-1, D_t:= D_t+value_neigh]
  d[index+1, D_t:= D_t+value_neigh]
  return(d)
}
for (i in ii) {
  data <- smooth_value_deaths(i, data)
}


data[rel_deviation_deaths > 2.5 & abs_deviation_deaths > 30]
data[rel_deviation_cases > 2.5 & abs_deviation_cases > 300]


data[, c("smooth_cases", "smooth_deaths", "abs_deviation_cases", 
         "abs_deviation_deaths", "rel_deviation_cases",
         "rel_deviation_deaths") := NULL]



ggplot(data)+
  geom_line(aes(date, D_t*20), color = "darkred")+
  geom_line(aes(date, repC_t),color = "darkblue")+
  facet_wrap(~country, scales = "free")




##############################################
# define start of the epidemic in each country

# delte rows before 23-01. NAs
data <- data[date > "2020-01-22"]

# start is 30 days before we observe a sum at least 10 deaths
data[, cum_deaths := cumsum(D_t), by = country]

start_dates <- data[cum_deaths >= 10, head(date, 1)-30,by = country]


# could do this line son data as well but I am superstitious o_O
data_tmp <- copy(data)
for (cc in start_dates$country){
  # print(cc)
  start_tmp <- start_dates[country == cc, V1]  
  data_tmp <- data_tmp[!(country == cc & date < start_tmp)]
}


data <- data_tmp
data[,cum_deaths := NULL]


ggplot(data) + geom_line(aes(date, repC_t, color = country))+
  geom_line(aes(date, alpha*100000, color = country))+
  geom_line(aes(date, delta*100000, color = country), linetype = "dashed")


data[, first(date), by = country]


####################################################
############# assign piD
####################################################

pi_D <- fread("../data/ifr_t_m_shifted.csv")
pi_D[country == "United Kingdom", country := "UnitedKingdom"]
pi_D[, date := as.Date(date)]
pi_D[, ifr_t_m_voc_adjusted := NULL]


data <- merge(data, pi_D, all.x = T, by = c("country","date"))



##################################
################# add NPIs
##################################

# load interventions data
ints <- fread("https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/OxCGRT_latest.csv")

# rename Czech Republic
ints[CountryName == "Czech Republic", CountryName := "Czechia"]
ints <- ints[CountryName %in% countries]
ints$CountryName %>% unique %>% length
# delete all where RegionName is empty
ints <- ints[RegionName == ""]


# setting 1
# workplacesClosed, smallGatherings, mediumGatherings, contactsTraced, 
# lockdown, secondLockdown, schoolsClosed, secondSchoolsclosed


# smallGatherings, mediumGatherings, contactsTraced, lockdown, secondLockdown,
# schoolsClosed, secondSchoolsclosed)

# setting 3
# generalBehavioralChanges, smallGatherings, mediumGatherings, contactsTraced, lockdown, secondLockdown, schoolsClosed




################ define interventions
ints[, date := as.Date(as.character(Date), format = "%Y%m%d")]
cols_setting <- names(ints)[c(1,62,7:20,33)]
ints <- ints[, ..cols_setting]


interventions_data <- copy(ints)

################################################################################
# closingSchools

interventions_data[, closingSchools := nafill(as.numeric(`C1_School closing` >= 2 & C1_Flag==1), "locf")]


ggplot(interventions_data)+
  geom_line(aes(date, `C1_School closing`))+
  geom_line(aes(date, C1_Flag), color=2)+
  geom_line(aes(date, closingSchools+0.1), color=3)+
  facet_wrap(~CountryName)


################################################################################
# closingWorkplaces

# no workplace variable, but could be something like this:
# interventions_data[, closingWorkplaces := as.numeric(`C2_Workplace closing` >= 2 & C2_Flag==1)]



################################################################################
# cancelPublicEvents
# should be incorporated in restricted gatherings
# interventions_data[, cancelPublicEvents := nafill(as.numeric(`C3_Cancel public events` == 2 & C3_Flag==1), "locf")]
# 
# 
# ggplot(interventions_data)+
#   geom_line(aes(date, `C3_Cancel public events`), color="darkblue")+
#   geom_line(aes(date, C3_Flag), color="steelblue")+
#   geom_line(aes(date, cancelPublicEvents+0.1), color=1)+
#   facet_wrap(~CountryName)



################################################################################
# restrictGatherings
interventions_data[, restrictGatherings := nafill(as.numeric(`C4_Restrictions on gatherings` == 4 & C4_Flag==1),"locf")]


ggplot(interventions_data)+
  geom_line(aes(date, `C4_Restrictions on gatherings`), color="darkblue")+
  geom_line(aes(date, C4_Flag), color="steelblue")+
  geom_line(aes(date, restrictGatherings+0.1), color=1)+
  facet_wrap(~CountryName)


################################################################################
# restrictPublicTransport

# not in the model, could be like this:
# interventions_data[, restrictPublicTransport := as.numeric(`C5_Close public transport` >= 1 & C5_Flag==1)]
# 
# 
# ggplot(interventions_data)+
#   geom_line(aes(date, `C5_Close public transport`), color="darkblue")+
#   geom_line(aes(date, C5_Flag), color="steelblue")+
#   geom_line(aes(date, restrictPublicTransport+0.1), color=1)+
#   facet_wrap(~CountryName)



################################################################################
# lockdown


interventions_data[, lockdown := nafill(as.numeric(`C6_Stay at home requirements` >= 2 & C6_Flag==1), "locf")]
# adpat lockdown for Germany 
interventions_data[CountryName == "Germany", lockdown := nafill(as.numeric(`C6_Stay at home requirements` >= 2), "locf")]

ggplot(interventions_data)+
  geom_line(aes(date, `C6_Stay at home requirements`), color="darkblue")+
  geom_line(aes(date, C6_Flag), color="steelblue")+
  geom_line(aes(date, lockdown+0.1), color=1)+
  facet_wrap(~CountryName)




################################################################################
# generalBehavioralChanges

# define starts
start_dates <- interventions_data[
  (closingSchools > 0 | restrictGatherings > 0| lockdown > 0), # | cancelPublicEvents > 0
  date, by = CountryName
  ][, first(date),  by = CountryName]

# assign variable
interventions_data[, generalBehavioralChanges := 0]
for (cc in start_dates$CountryName) {
  interventions_data[CountryName == cc & date > start_dates[CountryName == cc, V1],
          generalBehavioralChanges := 1]
}


ggplot(interventions_data)+
  geom_line(aes(date, generalBehavioralChanges))+
  facet_wrap(~CountryName)


################################################################################
# subsequentLockdown
  


define_subsequent_intervention <- function(df, country_sel, intervention){
  sec_int <- paste0("subsequent",stringr::str_to_title(intervention))
  if (sum((df[CountryName == country_sel][[intervention]] %>% diff()) != 0) >= 3){
    break_date <- df[CountryName == country_sel][which(c(FALSE,(df[CountryName == country_sel][[intervention]] %>% diff()) != 0))[3]-1]$date
    break_date <- df[CountryName == country_sel][which(c(FALSE,(df[CountryName == country_sel][[intervention]] %>% diff()) != 0))[3]-1]$date
    vals <- df[CountryName == country_sel & date >=break_date][[intervention]]
    df[CountryName == country_sel & date >=break_date, (sec_int) := vals]
    df[CountryName == country_sel & date >=break_date, (intervention) := 0]
  } else {
    df[CountryName == country_sel, (sec_int) := 0]
  }
  # note that this only implements second interventions IF the there IS an intervention
  # for example Finland had no lockdown at all
  return(df)
}

for (cc in unique(interventions_data$CountryName)){
  print(cc)
  interventions_data <- define_subsequent_intervention(interventions_data, cc, "lockdown")
}

interventions_data[is.na(subsequentLockdown), subsequentLockdown := 0]


interventions_data %>% ggplot()+
  geom_line(aes(date, subsequentLockdown+0.1), color=2)+ 
  geom_line(aes(date, lockdown), color = 1)+ 
  facet_wrap(~CountryName)





# plot NPIs
ggplot(interventions_data)+ 
  geom_line(aes(date, closingSchools+0.01), color="red")+
  # geom_line(aes(date, cancelPublicEvents+0.02), color="blue")+
  geom_line(aes(date, restrictGatherings+0.03), color="green")+
  geom_line(aes(date, lockdown+0.04), color="yellow")+
  geom_line(aes(date, subsequentLockdown+0.06), color="orange")+
  geom_line(aes(date, generalBehavioralChanges+0.05), color="black")+
  facet_wrap(~CountryName)


# delete old names
cols_keep <- c("CountryName", "date","closingSchools","restrictGatherings", "lockdown", "subsequentLockdown", "generalBehavioralChanges") # "cancelPublicEvents" 

interventions_data <- interventions_data[, ..cols_keep]

# NAs?
sapply(interventions_data, function(x) sum(is.na(x)))
# rename UK and columnname
interventions_data[CountryName == "United Kingdom", CountryName := "UnitedKingdom"]
setnames(interventions_data, "CountryName", "country")


tmp_npis <- melt(interventions_data, id.vars = c("country", "date"))
match_table <- tmp_npis$variable %>% unique()
tmp_npis[, value := value*match(variable, match_table)]
tmp_npis[, value := ifelse(value == 0, NA, value)]

ggplot(tmp_npis, aes(date, value, color = variable))+
  geom_line()+
  facet_wrap(~country) +
  theme(axis.text.x = element_text(angle=45, vjust=0.5))+
  scale_y_continuous(labels = rep("", 5))
  # + scale_y_continuous(labels = match_table)
  



####################################################
############# merge with rest
####################################################

data <- merge(data, interventions_data, all.x = T)



####################################################
############# assign seasonality
####################################################


define_dummy <- function(x, months){
   return(as.numeric(x %in% months))
}

define_variable <- function(d, name, months){
  d[,(name) := define_dummy(month(date), months)]
  return(d)
}

define_variable(data, "winter", c(12,1,2))
define_variable(data, "spring", c(3,4,5))
define_variable(data, "summer", c(6,7,8))
define_variable(data, "autumn", c(9,10,11))




source("plot_time_series.R")


plot_time_series(
  dataset = data[country == "Germany"],
  interventions = interventions_data %>% names %>% extract(-(1:2)),
  scale_deaths = 20,
  shifttext = 2000
)
plot_time_series(
  dataset = data[country == "France"],
  interventions = interventions_data %>% names %>% extract(-(1:2)),
  scale_deaths = 20,
  shifttext = 3000,
  shiftstart = 8000
)




# check data, NAs only for H(icu)_t valid
sapply(data, function(x) sum(is.na(x)))

data[, head(.SD, 1), by = country]
data[, tail(.SD, 1), by = country]


# date <- Sys.Date() %>% stringr::str_remove_all("-")
# (write_path = paste0("data/final_data",date,".csv"))


(write_path = "../data/real_data.csv")

fwrite(x = data, file = write_path)
