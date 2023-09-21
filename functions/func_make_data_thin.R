make_data_thin<-function (s,out_class='spocc',out_path) {

  time_to_decimal <- function(x) {
    x <- lubridate::hms(x, quiet = TRUE)
    lubridate::hour(x) + lubridate::minute(x) / 60 + lubridate::second(x) / 3600
  }
  
  
  if(out_class=='spocc') {
    dggs_thin=dgconstruct(area = 0.5*0.5)
  }
  
  if(out_class=='unmarked') {
    dggs_thin=dgconstruct(area = 1*1)
  }
  
  #####读取数据####
  sp=all_sp[s] %>% gsub('.csv','',.)
  
  if(!dir.exists(out_path)) {
    dir.create(out_path)
  }
    
  target_file<-paste0(out_path,'/',sp,'.csv')
  
  if(file.exists(target_file)){
    return()
  }
  
  ebird_zf<-fread(all_sp[s],
                      select = c("latitude","longitude","site",'year','cell','sampling_event_identifier','locality_id',
                                 "observation_date","time_observations_started","observer_id",
                                 "duration_minutes", "effort_distance_km", 'protocol_type',
                                 "number_observers",'species_observed','observation_count'))

  if(nrow(ebird_zf)==0) return() 
  
  ebird_zf <- ebird_zf %>%
    mutate(observation_count = ifelse(observation_count == "X",
                                      "1", observation_count
    ))
  
  # remove records where duration is 0
  ebird_zf <- filter(ebird_zf, duration_minutes > 0)
  
  # group data by site and sampling event identifier
  # then, summarise relevant variables as the sum

  # bind rows combining data frames, and filter
  ebird_zf  <- ebird_zf  %>%
    filter(
      duration_minutes <= 300,
      effort_distance_km <= 5,
      number_observers <= 10
    )
  
  ebird_zf <- mutate(ebird_zf,
                        pres_abs = observation_count >= 1,
                        decimalTime = time_to_decimal(time_observations_started)
  )
  
  
  ebird_zf <- ebird_zf %>%
    mutate(effort_distance_km = if_else(
      effort_distance_km == 0 &
        protocol_type == "Stationary",
      0.1, effort_distance_km
    ))
  
  ebird_zf <-
    ebird_zf %>%
    mutate(
      min_obs_started = as.integer(
        as.difftime(
          time_observations_started,
          format = "%H:%M:%S", units = "mins"
        )
      )
    )
  
  # Adding the julian date to the dataframe
  ebird_zf  <- ebird_zf  %>%
    mutate(julian_date = lubridate::yday(observation_date))
  
  # recode julian date to model it as a linear predictor
  ebird_zf  <- ebird_zf  %>%
    mutate(
      newjulianDate =
        case_when(
          (julian_date >= 334 & julian_date) <= 365 ~
            (julian_date - 333),
          (julian_date >= 1 & julian_date) <= 152 ~
            (julian_date + 31)
        )
    ) %>%
    drop_na(newjulianDate)
  
  # recode time observations started to model it as a linear predictor
  ebird_zf  <- ebird_zf  %>%
    mutate(
      newmin_obs_started = case_when(
        min_obs_started >= 300 & min_obs_started <= 720 ~
          abs(min_obs_started - 720),
        min_obs_started >= 720 & min_obs_started <= 1140 ~
          abs(720 - min_obs_started)
      )
    ) %>%
    drop_na(newmin_obs_started)
  
  
  ####thin data####
  
  ebird_zf <- group_by(ebird_zf, locality_id) %>%
    mutate(tot_effort = length(sampling_event_identifier)) %>%
    ungroup()

  
  ebird_zf$dgg_thin<-dgGEO_to_SEQNUM(dggs_thin,ebird_zf$longitude %>% as.numeric()
                                           ,ebird_zf$latitude %>% as.numeric())[[1]] #%>%
  # paste0("_",ebird_zf$year)

  ebird_zf <- ebird_zf %>%
    group_by(dgg_thin) %>%
    dplyr::filter(tot_effort == max(tot_effort)) %>%
    # there may be multiple rows with max effort, select first
    dplyr::filter(row_number() == 1)
  
  data_absences<-ebird_zf %>% dplyr::filter(species_observed==0)

  data_prsences<-ebird_zf %>% dplyr::filter(species_observed==1)

  data_thin_absences <- map_if(
    .x = split(x = ungroup(data_absences), f = ungroup(data_absences)$site),
    .p = function(x) {
      nrow(x) > 10
    },
    .f = function(x) sample_n(x, 10, replace = FALSE)
  ) %>% bind_rows()

  thin_ebird_zf<-bind_rows(data_thin_absences,data_prsences)

  
  rm(list = c('ebird_zf','data_prsences','data_absences'))
  
  gc()
  
  thin_ebird_zf <- dplyr::group_by_at(thin_ebird_zf, 'site') %>%
    dplyr::mutate(.obs_id = dplyr::row_number()) %>%
    dplyr::ungroup()
  

  # all_year_site<-select(thin_ebird_zf,year,cell,site)

  if(out_class=='unmarked') {
    
    all_cell_year<-crossing(year=2001:2020,
                            cell=thin_ebird_zf$cell %>% unique(),
                            .obs_id= 1:max(thin_ebird_zf$.obs_id))
    all_cell_year$site<-paste0(all_cell_year$cell,'_',all_cell_year$year)
    
    
    thin_ebird_zf<-left_join(all_cell_year,thin_ebird_zf,by=c('year','cell','.obs_id'))
  }
  
  
  fwrite(thin_ebird_zf,target_file)
  rm(thin_ebird_zf)
  gc()
}


