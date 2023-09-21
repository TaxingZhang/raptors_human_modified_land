library(auk)
library(dplyr)

raptor<-get_ebird_taxonomy()
raptor<-raptor[grepl('Cathartiformes|Accipitriformes|Strigiformes|Cariamiformes|Falconiformes',raptor$order) & 
                 raptor$category=='species',]
raptor<-unique(raptor$scientific_name)

# f_in_ebd <- file.path("F:/Faster study/raptor/ebd_relMar-2021.txt")
# f_in_sampling <- file.path("F:/Faster study/raptor/ebd_sampling_relNov-2021.txt")
# 
# #run filters using auk packages
# ebd_filters <- auk_ebd(f_in_ebd, f_in_sampling) %>%
#   auk_species(raptor[251:571]) %>%
#   auk_protocol(protocol = c("Stationary", "Traveling")) %>%
#   auk_date(c("2001-01-01", "2020-12-31")) %>%
#   auk_complete()
# 
# # check filters
# ebd_filters
# 
# # specify output location and perform filter
# f_out_ebd <- "F:/Faster study/raptor/ebd_raptor_2001-2020_251-571.txt"
# f_out_sampling <- "F:/Faster study/raptor/sampling_raptor_2001-2020.txt"
# 
# auk_filter(ebd_filters,f_out_ebd,f_out_sampling,overwrite = T,filter_sampling=F)

library(data.table)
library(dplyr)
columnsOfInterest <- c(
  "global_unique_identifier", "scientific_name", "observation_count",
  "locality", "locality_id", "locality_type", "latitude",
  "longitude", "observation_date",
  "time_observations_started", "observer_id",
  "sampling_event_identifier", "protocol_type",
  "duration_minutes", "effort_distance_km", "effort_area_ha",
  "number_observers") %>% toupper() 
columnsOfInterest<-gsub('_',' ',columnsOfInterest)

# ebd_1<-fread("F:/Faster study/raptor/ebd_raptor_2001-2020_1-250.txt",select=columnsOfInterest)
# ebd_2<-fread("F:/Faster study/raptor/ebd_raptor_2001-2020_251-571.txt")

ebd<-c("F:/Faster study/raptor/ebd_raptor_2001-2020_1-250.txt",
  "F:/Faster study/raptor/ebd_raptor_2001-2020_251-571.txt") %>%
  lapply(function(x){fread(x,select=columnsOfInterest)}) %>%
  rbindlist()    

names <- names(ebd) %>%
  stringr::str_to_lower() %>%
  stringr::str_replace_all(" ", "_")

setnames(ebd, names)


ebdSpSum <- ebd[, .(
  nSp = .N,
  totSoiSeen = length(intersect(scientific_name, raptor))
),
by = list(sampling_event_identifier)
]

library(lubridate)
time_to_decimal <- function(x) {
  x <- lubridate::hms(x, quiet = TRUE)
  lubridate::hour(x) + lubridate::minute(x) / 60 + lubridate::second(x) / 3600
}

ebd[, `:=`(
  decimalTime = time_to_decimal(time_observations_started),
  julianDate = yday(as.POSIXct(observation_date))
)]

ebdEffChk <- setDF(ebd) %>%
  mutate(year = year(observation_date)) %>%
  distinct(
    sampling_event_identifier, observer_id,
    year,protocol_type,
    duration_minutes, effort_distance_km, effort_area_ha,
    longitude, latitude,
    locality, locality_id,
    decimalTime, julianDate, number_observers
  ) %>%
  # drop rows with NAs in cols used in the model
  tidyr::drop_na(
    sampling_event_identifier, observer_id,
    duration_minutes, decimalTime, julianDate
  )

# remove ebird data
rm(ebd)
gc()

# 3. join to covariates and remove large groups (> 10)
ebdChkSummary <- inner_join(ebdEffChk, ebdSpSum)
fwrite(ebdChkSummary, file = "F:/Faster study/raptor/ebdChkSummary.csv")


library(sf)
library(data.table)
library(dggridR)

ebdChkSummary<-fread('F:/Faster study/raptor/ebdChkSummary.csv')

dggs<-dgconstruct(area = 2500)

system.time({ebdChkSummary[, `:=`(
  cell = dgGEO_to_SEQNUM(dggs,longitude, latitude)$seqnum
)]})

list_year_max_landcover<-list()

for(y in 1:20) {
  
  # y=1
  load(paste0("D:/work/QTP/study_of_raptor/hexs/hex_landcover/",y+2000,".RData"))
  
  year_landcover<-bind_rows(list_year_landcover)
  
  year_landcover$max_class<-apply(year_landcover[,!(colnames(year_landcover)=='year'|colnames(year_landcover)=='cell')],
                            1,function(x){which.max(x) %>% names()})
  
  list_year_max_landcover[[y]]<-select(year_landcover,year,cell,max_class)
  
  
  print(y)
  
}
year_max_landcover<-bind_rows(list_year_max_landcover)

ebdChkSummary<-left_join(ebdChkSummary,year_max_landcover,by=c('year','cell'))

# change names for easy handling
setnames(ebdChkSummary, c(
  "locality","locality_id", "latitude", "longitude", "observer", "sei",
  "protocol","duration", "distance", "area", "nObs", "decimalTime", "julianDate",
  "year", "nSp", "nSoi",'cell','landcover'
))

# count data points per observer
# obscount <- count(ebdChkSummary, observer) %>%
#   filter(n >= 3)

# make factor variables and remove obs not in obscount
# also remove 0 durations
ebdChkSummary <- ebdChkSummary %>%
  mutate(
    distance = ifelse(is.na(distance), 0, distance),
    duration = if_else(is.na(duration), 0.0, as.double(duration))
  ) %>%
  # filter(
  #   observer %in% obscount$observer,
  #   duration > 0,
  #   duration <= 300,
  #   nSoi >= 0,
  #   distance <= 5,
  #   !is.na(nSoi)
  # ) %>%
  mutate(
    landcover = as.factor(landcover),
    observer = as.factor(observer)
  ) %>%
  tidyr::drop_na(landcover)

# editing julian date to model it in a linear fashion
unique(ebdChkSummary$julianDate)

ebdChkSummary <- ebdChkSummary %>%
  mutate(
    newjulianDate =
      case_when(
        julianDate >= 334 & julianDate <= 365 ~ (julianDate - 333),
        julianDate >= 1 & julianDate <= 152 ~ (julianDate + 31)
      )
  ) %>%
  tidyr::drop_na(newjulianDate)


fwrite(ebdChkSummary, file = "F:/Faster study/raptor/ebdChkSummary_landcover_without_filter.csv")


library(data.table)

ebdChkSummary<-fread('F:/Faster study/raptor/ebdChkSummary_landcover_without_filter.csv')


# ebdChkSummary$landcover<-as.character(ebdChkSummary$landcover)
library(lmerTest)
library(dplyr)
system.time({
  modObsExp <- glmer(nSoi ~ duration + sqrt(duration) +
                       # landcover + 
                       protocol+ year +
                       sqrt(decimalTime) +
                       I((sqrt(decimalTime))^2) +
                       log(newjulianDate) +
                       I((log(newjulianDate)^2)) +
                       (1 | observer) + (0 + duration | observer),
                     data = ebdChkSummary,#[sample(1:12125730,10000000)], 
                     family = "poisson",nAGQ=0,verbose=2)
})

# make df with means
observer <- unique(ebdChkSummary$observer)

# 求众数
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

getmode(ebdChkSummary$landcover)


load("D:/work/QTP/study_of_raptor/Manuscript/CB/observer_model.RData")
# predict at 60 mins on the most common landcover (deciduous forests)
dfPredict_sta <- ebdChkSummary %>%
  summarise_at(vars(duration, decimalTime, newjulianDate), list(~ mean(.))) %>%
  mutate(duration = 60, protocol="Stationary", year=2011,landcover = as.factor(210)) %>%
  tidyr::crossing(observer)

# run predict from model on it
dfPredict_sta <- mutate(dfPredict_sta,
                    score = predict(modObsExp,
                                    newdata = dfPredict_sta,
                                    type = "response",
                                    allow.new.levels = TRUE
                    )
) %>%
  mutate(score = scales::rescale(score))


# predict at 60 mins on the most common landcover (deciduous forests)
dfPredict_tra <- ebdChkSummary %>%
  summarise_at(vars(duration, decimalTime, newjulianDate), list(~ mean(.))) %>%
  mutate(duration = 60, protocol="Traveling", year=2011,landcover = as.factor(210)) %>%
  tidyr::crossing(observer)

# run predict from model on it
dfPredict_tra <- mutate(dfPredict_tra,
                        score = predict(modObsExp,
                                        newdata = dfPredict_tra,
                                        type = "response",
                                        allow.new.levels = TRUE
                        )
) %>%
  mutate(score = scales::rescale(score))

mean_observer_score<-data.frame(observer=dfPredict_sta$observer,
                                mean_score=cbind(dfPredict_sta$score,dfPredict_tra$score) %>% apply(1,mean))




