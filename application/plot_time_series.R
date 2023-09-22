# This script implements the plot functions to viz the time series
# 1. function to plot everything in one plot
# 2. the same as 1. but without orhwe stuff - not in the manuscript, only used in presentations for didactical reasons
# 3. script to plot the NPIs  
require(data.table)
require(ggplot2)
require(magrittr)

plot_time_series <- function(dataset, interventions, scale_deaths=20, 
                             shifttext = 2000, shiftstart=1000, hide_legend=T, xpos_shift=25){
  #' Plot function to plot various effects
  #' 
  #' @param dataset A data.table obejct containing all informations
  #' @param interventions A character vector with names of the NPIs
  #' @param shifttext Defines each offset for the y distance between NPIs
  #' @param shiftstart The offset to start from the time series
  #' 
  #' @return A ggplot object
  
  
  shift = shiftstart
  dataset <- copy(dataset) # just to be sure that nothing outside the function gets overwritten
  
  for(intervention in interventions){
    shift = shift + shifttext
    id <- dataset[[intervention]]  == 0
    if (sum(id) == length(id)) {
      id[1] <- FALSE
      dataset[1, (intervention)] <- 1
    }
    dataset[id, (intervention) := NA]
    dataset[[intervention]] <- dataset[[intervention]] + max(dataset$repC_t)+shift
  }
  
  dataset[, seasons_pos := max(dataset$repC_t) + shiftstart]
  
  # only required if you want to plot a ribbon
  for(season in c("summer", "autumn","winter","spring")){
    id <- dataset[[season]]  == 0
    dataset[id, (season) := NA]
  }
  dataset[, winter := winter * seasons_pos]
  dataset[, spring := spring * seasons_pos]
  dataset[, summer := summer * seasons_pos]
  dataset[, autumn := autumn * seasons_pos]
  
  
  
  g <- ggplot(dataset) + 
    geom_line(aes(date, repC_t, color = "cases"),alpha = 0.9)+
    geom_line(aes(date, D_t*scale_deaths, color = "deaths"), alpha = 0.9)+
    geom_line(aes(date, H_t, color = "hospital"),alpha = 0.5)+
    geom_line(aes(date, Hicu_t, color = "ICU"),alpha = 0.9)+
    geom_line(aes(date, (1-alpha-delta)*max(repC_t), color = "original"),alpha = 0.2)+
    geom_line(aes(date, (alpha)*max(repC_t),color = "alpha"),alpha = 0.2)+
    geom_line(aes(date, (delta)*max(repC_t),color = "delta"),alpha = 0.2)+
    geom_line(aes(date, first_vaccination*max(repC_t), color="Vaccination 1"),alpha = 0.4, linetype="dashed")+
    geom_line(aes(date, second_vaccination*max(repC_t), color="Vaccination 2"),alpha = 0.4, linetype="dashed")+
    ylab("reported time series:\n cases, deaths, hospital/icu occupancy")+
    theme_bw()
  
  
  xpos <- dataset$date %>% min %>% subtract(xpos_shift)
  
  for(intervention in interventions){
    pos <- dataset[[intervention]] %>% na.omit() %>% extract(1)
    
    g <- g + geom_line(aes_string("date", intervention), size = 2.5, color = "steelblue", lineend="round", alpha=0.5)+
      annotate("text", x = xpos, y = pos, label = intervention, color = "steelblue", hjust = 1, size = 3.7)
  }
  
  # lines
  g <- g +
    geom_line(aes(date, winter,color="winter"), size = 2.5,lineend="round", alpha = .3)+
    geom_line(aes(date, spring,color="spring"), size = 2.5,lineend="round", alpha = .3)+
    geom_line(aes(date, summer,color="summer"), size = 2.5,lineend="round", alpha = .3)+
    geom_line(aes(date, autumn,color="autumn"), size = 2.5,lineend="round", alpha = .3)+
    annotate("text", x = xpos, y = dataset$seasons_pos[1], label = "Seasons", color = "steelblue", hjust = 1, size = 3.7)
  
  
  
  max_val <- dataset[[intervention]] %>% max(na.rm=T)
  g <- g + coord_cartesian(ylim = c(0, max_val+1), 
                           xlim = c(min(dataset$date)+5, max(dataset$date)+5), 
                           clip = "off")+
    scale_y_continuous(breaks = dataset$repC_t %>% max %>% seq(0,.,by=10000),
                       sec.axis = sec_axis(~./max(dataset$repC_t),breaks=c(0,0.5,1),
                                           name = "proportion vaccinations & variants of concern"))
  
  g <- g+
    theme(legend.position = 'bottom')+
    scale_color_manual(name="Series", values=c("original"="black",
                                               "alpha"="darkgreen",
                                               "delta"="darkred",
                                               "Vaccination 1" = "deeppink1",
                                               "Vaccination 2" = "deeppink3",
                                               "cases" = "steelblue1",
                                               "deaths" = "grey30",
                                               "hospital" = "blue3",
                                               "ICU" = "darkblue",
                                               "autumn" = "orange", "winter"="royalblue3",
                                               "spring"="green3", "summer"="gold"
    ))
  if (hide_legend){
    g <- g + theme(legend.position="none") + theme(plot.margin = unit(c(0.3,.3,.3, 3), "cm"))+
      theme(axis.title.y = element_text(vjust=5))
    
  }
  return(g)
}


d <- fread("../data/real_data.csv")
ints <- names(d)[15:19]
g <- plot_time_series(d[country == "UnitedKingdom"], ints, xpos_shift = 30)


# please note, that the legend in the manuscript is created by hand as ggplot cannot arrange it in the expected way
ggsave("plot_all_series_UK.pdf",width = 15, height = 8.5)






#################################################################################
# same plot but with invisible lines and no sec axis
#(not used in the paper, but for didactical reasons usefull e.g. in presentations)
#################################################################################

plot_time_series_less_stuff <- function(dataset, interventions, scale_deaths=20, 
                             shifttext = 2000, shiftstart=1000, hide_legend=T, xpos_shift=25){
  shift = shiftstart
  dataset <- copy(dataset) # just to be sure that nothing outside the function gets overwritten
  
  for(intervention in interventions){
    shift = shift + shifttext
    id <- dataset[[intervention]]  == 0
    if (sum(id) == length(id)) {
      id[1] <- FALSE
      dataset[1, (intervention)] <- 1
    }
    dataset[id, (intervention) := NA]
    dataset[[intervention]] <- dataset[[intervention]] + max(dataset$repC_t)+shift
  }
  
  dataset[, seasons_pos := max(dataset$repC_t) + shiftstart]
  # for(ss in c("summer", "autumn","winter","spring")){
  #   id <- dataset[[ss]]  == 0
  #   dataset[id, seasons := ss]
  # }
  
  # only required if you want to plot a ribbon
  for(season in c("summer", "autumn","winter","spring")){
    id <- dataset[[season]]  == 0
    dataset[id, (season) := NA]
  }
  dataset[, winter := winter * seasons_pos]
  dataset[, spring := spring * seasons_pos]
  dataset[, summer := summer * seasons_pos]
  dataset[, autumn := autumn * seasons_pos]
  
  
  
  g <- ggplot(dataset) + 
    geom_line(aes(date, repC_t, color = "cases"),alpha = 0.9)+
    geom_line(aes(date, D_t*scale_deaths, color = "deaths"), alpha = 0.9)+
    geom_line(aes(date, H_t, color = "hospital"),alpha = 0.7)+
    geom_line(aes(date, H_t, color = "ICU"),alpha = 0.7)+
    geom_line(aes(date, (1-alpha-delta)*max(repC_t), color = "original"),alpha = 0)+
    geom_line(aes(date, (alpha)*max(repC_t),color = "alpha"),alpha = 0)+
    geom_line(aes(date, (delta)*max(repC_t),color = "delta"),alpha = 0)+
    geom_line(aes(date, first_vaccination*max(repC_t), color="Vaccination 1"),alpha = 0, linetype="dashed")+
    geom_line(aes(date, second_vaccination*max(repC_t), color="Vaccination 2"),alpha = 0, linetype="dashed")+
    ylab("reported time series:\n cases, deaths, hospital/icu occupancy")
  
  
  xpos <- dataset$date %>% min %>% subtract(xpos_shift)
  
  for(intervention in interventions){
    pos <- dataset[[intervention]] %>% na.omit() %>% extract(1)
    
    g <- g + geom_line(aes_string("date", intervention), size = 2.5, color = "steelblue", lineend="round", alpha=0)+
      annotate("text", x = xpos, y = pos, label = intervention, color = "steelblue", alpha=0)
  }
  g <- g +
    geom_line(aes(date, winter,color="winter"), size = 2.5,lineend="round", alpha = 0)+
    geom_line(aes(date, spring,color="spring"), size = 2.5,lineend="round", alpha = 0)+
    geom_line(aes(date, summer,color="summer"), size = 2.5,lineend="round", alpha = 0)+
    geom_line(aes(date, autumn,color="autumn"), size = 2.5,lineend="round", alpha = 0)+
    annotate("text", x = xpos, y = dataset$seasons_pos[1], label = "Seasons", color = "steelblue", alpha=0)
  
  
  
  max_val <- dataset[[intervention]] %>% max(na.rm=T)
  g <- g + coord_cartesian(ylim = c(0, max_val+1), 
                           xlim = c(min(dataset$date)-5, max(dataset$date)+5), 
                           clip = "off")+
    scale_y_continuous(breaks = dataset$repC_t %>% max %>% seq(0,.,by=10000),
    )
  
  
  g <- g+
    theme(legend.position = 'bottom')+
    scale_color_manual(name="Series", values=c("original"="black",
                                               "alpha"="darkgreen",
                                               "delta"="darkred",
                                               "Vaccination 1" = "deeppink1",
                                               "Vaccination 2" = "deeppink3",
                                               "cases" = "steelblue1",
                                               "deaths" = "steelblue4",
                                               "hospital" = "steelblue2",
                                               "ICU" = "steelblue3",
                                               "autumn" = "orange", "winter"="royalblue3",
                                               "spring"="green3", "summer"="gold"
    ))
  if (hide_legend){
    g <- g + theme(legend.position="none")
    
  }
  return(g)
}







#################################################################################
# code to plot the interventions as timeline
#################################################################################



ddd <- fread("../data/real_data.csv")[,.(country, date, 
                                         closingSchools,
                                         restrictGatherings,
                                         lockdown,
                                         subsequentLockdown,
                                         generalBehavioralChanges
                                         #winter,
                                         #spring,
                                         #summer,
                                         #autumn
                                         )]

cols <- colnames(ddd)[3:7]

i <- 1
for(npi in cols){
  ddd[[npi]] <- ddd[[npi]] * 0.1*i
  i <- i+1
}

ddd <- melt(ddd, id.vars = c("country", "date")) 
ddd[value == 0, value := NA]

npi_names <- c("School closure", "Gatherings", "Lockdown", "Subsequent lockdown","General behavioral changes")

g <- ggplot(ddd)+
  geom_line(aes(date, value, color = variable),size=2, lineend = "round")+
  facet_wrap(~country)+
  scale_y_continuous(labels=npi_names)+
  theme(legend.position = "none")+
  ylab(element_blank())+
  xlab(element_blank())+
  theme(axis.text.x = element_text(angle=45, hjust = 0.8))

  
ggsave("plots_NPIs_europe.pdf", g, width = 10, height = 8)
