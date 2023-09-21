*This file contains all the processes from data preparation to results
output*

# 1 Filtering Species data

    library(auk)
    library(dplyr)
    library(rgdal)
    library(dggridR)
    library(lubridate)
    library(plyr)
    library(doParallel)
    library(foreach)
    library(data.table)
    library(purrr)
    library(stringr)
    library(tidyr)

## 1.1 Get Latin name of all raptors

    birds<-get_ebird_taxonomy()
    raptors<-birds[grepl('Cathartiformes|Accipitriformes|Strigiformes|Cariamiformes|Falconiformes',birds$order) & 
                     birds$category=='species',]
    raptors<-unique(raptors$scientific_name)

## 1.2 filter data during 2001 to 2020

    #ebd_relMar-2021.txt and sampling_relMar-2021 can obtain from https://ebird.org/

    ebd <- auk_ebd('.../ebd_relMar-2021.txt',
                   '.../sampling_relMar-2021.txt')

    ##set filter criteria
    ebd_filters <- ebd %>%
      auk_species(raptor) %>%
      auk_protocol(protocol = c("Stationary", "Traveling")) %>%
      auk_complete() %>%
      auk_date(date = c("2001-1-1", "2020-12-31"))

    ##output path
    f_ebd <- file.path('.../filtered_raptor_ebd_2001-2020.txt')
    f_sampling <- file.path('.../filtered_raptor_sampling_2001-2020.txt')

    ##start filter
    auk_filter(ebd_filters,f_ebd,f_sampling,overwrite = T,filter_sampling=T)

## 1.3 filter data of each raptors

    ##input path
    ebd <- auk_ebd('.../filtered_raptor_ebd_2001-2020.txt',
                   '.../filtered_raptor_sampling_2001-2020.txt')

    ##species loop: this step will take a lot of time
    ##if your hardware is good enough
    ##parallel computing is recommended

    for (s in 1:length(raptor)) {
      
      ##set filter criteria
      ebd_filters <- ebd %>%
        auk_species(raptor[s]) %>%
        auk_complete()
      
      ##output path
      f_ebd <- paste0('.../filtered_raptor_ebd/filtered_raptor_ebd_',raptor[s],'.txt')
      f_sampling <- paste0('.../filtered_raptor_sampling/filtered_raptor_sampling_',raptor[s],'.txt')
      
      ##start filter
      auk_filter(ebd_filters,f_ebd,f_sampling,overwrite = T,filter_sampling=T)
      
      print(s)
      
    } ##species loop

## 1.4 get zero filled data for modeling

    ## Using parallel computing, import functions from the functions folder
    source('.../functions/func_zero_filling.R')
    ## Set 10 parallel computations
    cl<- makeCluster(detectCores()-2)
    registerDoParallel(cl)

    foreach(
      s=1:length(raptors),
      .packages = c("auk","dggridR","dplyr",'lubridate'),
      .verbose=F
    ) %dopar% func_zero_filling(s,
                                f_ebd_path='filtered_raptor_ebd',
                                f_sampling_path='filtered_raptor_sampling',
                                expert_map_path='...',
                                output_path'.../raptor_zerofill_breeding')

    ## close parallel
    stopCluster(cl)

# 2 Calculate observer score (CCI)

The code of section 2-3 refer from *RAMESH V, GUPTE P R, TINGLEY M W, et
al. 2022. Using citizen science to parse climatic and land cover
influences on bird occupancy in a tropical biodiversity hotspot.
Ecography \[J\], 2022: e06075.*

## 2.1 Prepare species data

    library(auk)
    library(dplyr)
    library(data.table)

    columnsOfInterest <- c(
      "global_unique_identifier", "scientific_name", "observation_count",
      "locality", "locality_id", "locality_type", "latitude",
      "longitude", "observation_date",
      "time_observations_started", "observer_id",
      "sampling_event_identifier", "protocol_type",
      "duration_minutes", "effort_distance_km", "effort_area_ha",
      "number_observers") %>% toupper() 

    columnsOfInterest<-gsub('_',' ',columnsOfInterest)

    ebd<-fread(".../filtered_raptor_ebd_2001-2020.txt",select=columnsOfInterest)

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

## 2.2 Prepare checklists for observer score

    library(lubridate)

    ## add new columns of decimal time and julian date
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
      ## drop rows with NAs in cols used in the model
      tidyr::drop_na(
        sampling_event_identifier, observer_id,
        duration_minutes, decimalTime, julianDate
      )

    ## remove ebird data
    rm(ebd)
    gc()

    ## join to covariates
    ebdChkSummary <- inner_join(ebdEffChk, ebdSpSum)

## 2.3 Filter checklist data

    ## change names for easy handling
    setnames(ebdChkSummary, c(
      "locality","locality_id", "latitude", "longitude", "observer", "sei",
      "protocol","duration", "distance", "area", "nObs", "decimalTime", "julianDate",
      "year", "nSp", "nSoi",'cell','landcover'
    ))

    ## make factor variables and remove obs not in obscount
    ## also remove 0 durations
    ebdChkSummary <- ebdChkSummary %>%
      mutate(
        distance = ifelse(is.na(distance), 0, distance),
        duration = if_else(is.na(duration), 0.0, as.double(duration))
      )

    ## editing julian date to model it in a linear fashion
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

## 2.4 Model observer expertise

    library(lmerTest)
    library(dplyr)

    modObsExp <- glmer(nSoi ~ duration + sqrt(duration) +
                           protocol+ year +
                           sqrt(decimalTime) +
                           I((sqrt(decimalTime))^2) +
                           log(newjulianDate) +
                           I((log(newjulianDate)^2)) +
                           (1 | observer) + (0 + duration | observer),
                         data = ebdChkSummary, 
                         family = "poisson",nAGQ=0,verbose=2)

    ## make df with means
    observer <- unique(ebdChkSummary$observer)

    ## predict stationary
    dfPredict_sta <- ebdChkSummary %>%
      summarise_at(vars(duration, decimalTime, newjulianDate), list(~ mean(.))) %>%
      mutate(duration = 60, protocol="Stationary", year=2011) %>%
      tidyr::crossing(observer)

    ## run predict from model on it
    dfPredict_sta <- mutate(dfPredict_sta,
                        score = predict(modObsExp,
                                        newdata = dfPredict_sta,
                                        type = "response",
                                        allow.new.levels = TRUE
                        )
    ) %>%
      mutate(score = scales::rescale(score))


    ## predict traveling
    dfPredict_tra <- ebdChkSummary %>%
      summarise_at(vars(duration, decimalTime, newjulianDate), list(~ mean(.))) %>%
      mutate(duration = 60, protocol="Traveling", year=2011,landcover = as.factor(210)) %>%
      tidyr::crossing(observer)

    ## run predict from model on it
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

    fwrite(mean_observer_score,'.../mean_observer_score.csv)

# 3 Make data thinning

    library(data.table)
    source('.../functions/func_make_data_thin.R')

    setwd('.../raptor_zerofill_breeding')
    all_sp<-Sys.glob('*.csv')
    set.seed(1)

    ## set site grid
    dggs=dgconstruct(area = 55*55)

    ## Make data thinning for spOccupancy model
    ## Use parallel computing
    cl<- makeCluster(detectCores()-2)
    registerDoParallel(cl)

    foreach(
      s=1:length(all_sp), 
      .errorhandling = c("stop"),
      .packages = c('dggridR','dplyr','purrr','plyr','stringr','tidyr','auk','lubridate','data.table'),
      .verbose=T
    ) %dopar% make_data_thin(s,10,'spOcc','.../raptor_zerofill_breeding_thin_0.25km_spOcc')

    stopCluster(cl)

    ## Make data thinning for dynamic model
    ## Use parallel computing
    cl<- makeCluster(detectCores()-2)
    registerDoParallel(cl)

    foreach(
      s=1:length(all_sp), 
      .errorhandling = c("stop"),
      .packages = c('dggridR','dplyr','purrr','plyr','stringr','tidyr','auk','lubridate','data.table'),
      .verbose=T
    ) %dopar% make_data_thin(s,10,'unmarked','.../raptor_zerofill_breeding_thin_1km_unmarked')

    stopCluster(cl)

# 4 Link species data and occupancy covariates

We extract mean\_observer\_score covariates from the global data for the
period 2001-2020 The global data has been processed into RData files in
the folder: …/hex

    ##
    load('.../hexs/hex_pop/all_year_pop.RData',.GlobalEnv)
    all_year_pop<-all_year_v
    ##
    load('.../hexs/hex_npp/all_year_npp.RData'),.GlobalEnv)
    all_year_npp<-all_year_v
    ##
    load(paste0('.../hexs/hex_shdi/all_year_shdi.RData'),.GlobalEnv)
    all_year_shdi<-all_year_v
    ##
    rm(all_year_v)
    ##
    mean_observer_score<-fread(".../mean_observer_score.csv")
    mean_observer_score<-mutate(mean_observer_score,numObserver = str_extract(observer, "\\d+")) %>%
        dplyr::select(-observer)

    source('.../functions/func_link_data_cov.R')

    ## link data for spOccupancy
    func_link_cov(input_path='.../raptor_zerofill_breeding_thin_0.25km_spOcc',
                  output_path='.../raptor_zerofill_breeding_thin_0.25km_spOcc_cov')
    ## link data  for dynamic model
    func_link_cov(input_path='.../raptor_zerofill_breeding_thin_1km_unmarked',
                  output_path='.../raptor_zerofill_breeding_thin_1km_unmarked_cov')

# 5 Formatting data for model

    library(doParallel)
    library(parallel)
    source('.../functions/func_make_data_for_spOcc_unmarked.R')

    ## set site grid
    dggs=dgconstruct(area = 55*55)

    ## Formatting data for spOccupancy
    setwd('.../raptor_zerofill_breeding_thin_0.25km_spOcc_cov/')
    all_sp<-Sys.glob('*.csv')

    cl<- makeCluster(detectCores()-2)
    registerDoParallel(cl)

    foreach(
      s=1:length(all_sp), 
      .errorhandling = c("stop"),
      .packages = c('dggridR','dplyr','purrr','plyr','stringr','tidyr','auk','lubridate','data.table'),
      .verbose=T
    ) %dopar% make_data(s,'.../data_for_spOcc')

    stopCluster(cl)

    ## Formatting data for dynamic mod
    setwd('.../raptor_zerofill_breeding_thin_1km_unmarked_cov/')
    all_sp<-Sys.glob('*.csv')

    cl<- makeCluster(detectCores()-2)
    registerDoParallel(cl)

    foreach(
      s=1:length(all_sp), 
      .errorhandling = c("stop"),
      .packages = c('dggridR','dplyr','purrr','plyr','stringr','tidyr','auk','lubridate','data.table'),
      .verbose=T
    ) %dopar% make_data(s,'.../data_for_unmarked')

    stopCluster(cl)

# 6 Fit multi-season single-species occupancy models (MSSOMs)

    library(doParallel)
    library(parallel)
    source('.../functions/func_MSSOMs.R')

    setwd('.../data_for_spOcc')
    all_sp<-Sys.glob('*.RData')


    cl<- makeCluster(detectCores()-2)
    registerDoParallel(cl)

    foreach(
      s=1:length(all_sp), 
      .errorhandling = c("stop"),
      .packages = c('dplyr','spOccupancy'),
      .verbose=T
    ) %dopar% func_MSSOMs(s,output_path='.../spOcc_out')

    stopCluster(cl)

# 7 Fit Multi-season single-species dynamic occupancy models (MSSDOMs)

    library(doParallel)
    library(foreach)
    source('.../functions/func_MSSDOMs.R')

    setwd('.../data_for_unmarked')
    all_sp<-Sys.glob('*.RData')

    cl<- makeCluster(detectCores()-2)
    registerDoParallel(cl)

    foreach(
      s=1:length(all_sp), 
      .errorhandling = c("stop"),
      .packages = c('dplyr','unmarked'),
      .verbose=T
    ) %dopar% func_occ_model(s,'.../unmarked_out')

    stopCluster(cl)

# 8 Extract estimated coefficients from MSSOMs

    library(doParallel)
    library(foreach)
    library(data.table)
    source('.../functions/func_get_coef_from_MSSOMs')

    setwd('.../spOcc_out')
    all_sp<-Sys.glob('*.RData')

    ##
    cl<- makeCluster(detectCores()-2)
    registerDoParallel(cl)
    coefs<-foreach(
      s=1:length(all_sp),
      .errorhandling = 'remove',
      .packages = c('stats','dplyr','spOccupancy')
    ) %dopar% get_MSSOMs_coef(s)

    stopCluster(cl)

    coefs<-Reduce(function(x,y) {bind_rows(x,y)},coefs)

    fwrite(coefs,'.../coefs_MSSOMs.csv')

# 9 Calculate Bayesian p-value of MSSOMs

    library(spOccupancy)
    library(dplyr)

    setwd('.../spOcc_out')
    all_sp<-Sys.glob('*.RData')
    all_bayes_p<-list()

    ## Calculating Bayes p values takes up a lot of memory, so we use a loop instead of parallel computation
    for (s in length(all_sp)) {
      load(all_sp[s])
      ppc.out<-ppcOcc(out, fit.stat = 'freeman-tukey', group = 1)
      ppc.df<-ppc.out$fit.y.rep>ppc.out$fit.y
      all_bayes_p[[s]]<-mean(apply(ppc.df, 2, sum)/nrow(ppc.df))
      rm(ppc.out)
      gc()
      print(s)
    }

    mean(unlist(all_bayes_p))
    sd(unlist(all_bayes_p))/sqrt(length(all_bayes_p))

# 10 Extract estimated coefficients from MSSDOMs and calculate sensitivity of raptor occupancy probability to covariates

## 10.1 Extract estimated coefficients from MSSDOMs

    library(doParallel)
    library(foreach)
    library(data.table)
    source('.../functions/func_get_coef_from_MSSDOMs')

    setwd('.../unmarked_out')
    all_sp<-Sys.glob('*.RData')

    ##
    cl<- makeCluster(detectCores()-2)
    registerDoParallel(cl)
    foreach(
      s=1:length(all_sp),
      .errorhandling = 'remove',
      .packages = c('stats','dplyr','unmarked')
    ) %dopar% get_MSSDOMs_coef(s,'.../unmarked_coefs')

    stopCluster(cl)

## 10.2 Calculate sensitivity of raptor occupancy probability to covariates

    source('.../functions/func_calculate_sensitivity.R')

    setwd('.../unmarked_out')
    all_sp<-Sys.glob('*.RData')

    cl<- makeCluster(detectCores()-2)
    registerDoParallel(cl)
    senss<-foreach(
      s=1:length(all_sp),
      .packages = c('data.table','dplyr','plyr')
    ) %dopar% calculate_sens(s)

    stopCluster(cl)

    senss<-Reduce(function(x,y) {bind_rows(x,y)},senss)

    fwrite(senss,'.../senss_MSSDOMs.csv')

# 11 Fit PGLSs and plot figures

## 11.1 Fit PGLSs

    library(dplyr)
    library(plyr)
    library(doParallel)
    library(foreach)
    library(data.table)
    library(caper)
    source('.../functions/func_pgls.R')

    ## read data
    coefs<-fread('.../coefs_MSSOMs.csv')
    senss<-fread('.../sensss_MSSDOMs.csv')

    diet_raptor<-read.csv('.../data/AVONET1_BirdLife.csv')
    crosswalk<-read.csv('.../data/BirdLife-BirdTree crosswalk.csv')
    tree_sys<-read.csv('.../data/BLIOCPhyloMasterTax.csv')
    tree<-read.nexus('.../data/9993.nex')

    ## Combine data
    all_sp_coef_plot<-left_join(coefs,diet_raptor,by=c('species'='Species1')) %>%
      left_join(crosswalk,by=c('species'='Species1')) %>%
      left_join(tree_sys,by=c('Species3'='Scientific'))
    all_sp_coef_plot<-all_sp_coef_plot[all_sp_coef_plot$coefs_name!='(Intercept)',]
    all_sp_coef_plot$signf<-(all_sp_coef_plot$L_coef) * (all_sp_coef_plot$U_coef)
    all_sp_coef_plot$signf<-ifelse(all_sp_coef_plot$signf>0,'S','NS')
    all_sp_coef_plot$TipLabel[all_sp_coef_plot$signf>0] %>% unique() %>% length()

    ## Pruning phylogenetic tree
    group_species<-dplyr::select(all_sp_coef_plot,TipLabel,Order1) %>% unique()
    tree_dropped_tip<-tidytree::drop.tip(tree,tree$tip.label[!tree$tip.label %in% group_species$TipLabel])

    ## fit pgls

    mod_list1<-batch_pgls('LC_sett','coef')
    mod_list2<-batch_pgls('LC_agri','coef')
    mod_list3<-batch_pgls('pasture','coef')

    mod_list4<-batch_pgls('LC_sett','sens')
    mod_list5<-batch_pgls('LC_agri','sens')
    mod_list6<-batch_pgls('pasture','sens')

    all_mod_list<-list(mod_list1,mod_list2,mod_list3,mod_list4,mod_list5,mod_list6)

    save(all_mod_list,file = '.../all_mod_list.RData')

    ## print model output
    all_mod_list[[1]][['Mass']] %>% summary()
    all_mod_list[[2]][['Mass']] %>% summary()
    all_mod_list[[3]][['Mass']] %>% summary()

    all_mod_list[[4]][['Mass']] %>% summary()
    all_mod_list[[5]][['Mass']] %>% summary()
    all_mod_list[[6]][['Mass']] %>% summary()

    all_mod_list[[1]][['Range.size']] %>% summary()
    all_mod_list[[2]][['Range.size']] %>% summary()
    all_mod_list[[3]][['Range.size']] %>% summary()

    all_mod_list[[4]][['Range.size']] %>% summary()
    all_mod_list[[5]][['Range.size']] %>% summary()
    all_mod_list[[6]][['Range.size']] %>% summary()

    all_mod_list[[1]][['Trophic Niche']] %>% anova()
    all_mod_list[[2]][['Trophic Niche']] %>% anova()
    all_mod_list[[3]][['Trophic Niche']] %>% anova()

    all_mod_list[[4]][['Trophic Niche']] %>% anova()
    all_mod_list[[5]][['Trophic Niche']] %>% anova()
    all_mod_list[[6]][['Trophic Niche']] %>% anova()

    all_mod_list[[1]][['Habitat']] %>% anova()
    all_mod_list[[2]][['Habitat']] %>% anova()
    all_mod_list[[3]][['Habitat']] %>% anova()

    all_mod_list[[4]][['Habitat']] %>% anova()
    all_mod_list[[5]][['Habitat']] %>% anova()
    all_mod_list[[6]][['Habitat']] %>% anova()

    all_mod_list[[1]][['Migration']] %>% anova()
    all_mod_list[[2]][['Migration']] %>% anova()
    all_mod_list[[3]][['Migration']] %>% anova()

    all_mod_list[[4]][['Migration']] %>% anova()
    all_mod_list[[5]][['Migration']] %>% anova()
    all_mod_list[[6]][['Migration']] %>% anova()

## 11.2 Plot figures

    library(ggsci)
    library(ggplot2)
    library(cowplot)

    load('.../all_mod_list.RData')
    source('.../functions/func_plot.RData')
    col_list <- pal_npg(palette = c("nrc"), alpha = 1)(5)

    ps_sett<-batch_plot('LC_sett','coef')
    ps_agri<-batch_plot('LC_agri','coef')
    ps_pasture<-batch_plot('pasture','coef')

    ## Or
    # ps_sett<-func_plot('LC_sett','sens')
    # ps_agri<-func_plot('LC_agri','sens')
    # ps_pasture<-func_plot('pasture','sens')
    # 

    ## plot Fig. 6
    pp<-plot_grid(ps_sett[[1]]+theme(legend.position = 'none')+xlab(NULL)+ylab(NULL)+theme(axis.text.x = element_blank(),
                                                                                      axis.text.y = element_blank()),
                  ps_sett[[2]]+theme(legend.position = 'none')+xlab(NULL)+ylab(NULL)+theme(axis.text.x = element_blank(),
                                                                                      axis.text.y = element_blank()),
                  ps_agri[[1]]+theme(legend.position = 'none')+xlab(NULL)+ylab(NULL)+theme(axis.text.x = element_blank(),
                                                                                      axis.text.y = element_blank()),
                  ps_agri[[2]]+theme(legend.position = 'none')+xlab(NULL)+ylab(NULL)+theme(axis.text.x = element_blank(),
                                                                                      axis.text.y = element_blank()),
                  ps_pasture[[1]]+theme(legend.position = 'none')+xlab(NULL)+ylab(NULL)+theme(axis.text.x = element_blank(),
                                                                                       axis.text.y = element_blank()),
                  ps_pasture[[2]]+theme(legend.position = 'none')+xlab(NULL)+ylab(NULL)+theme(axis.text.x = element_blank(),
                                                                                       axis.text.y = element_blank()),
                  nrow =3,byrow = T)
    plot(pp)

    ## plot Fig. 7
    pp<-plot_grid(ps_sett[[5]]+theme(legend.position = 'none')+xlab(NULL)+ylab(NULL)+theme(axis.text.x = element_blank(),
                                                                                      axis.text.y = element_blank()),
                  ps_agri[[5]]+theme(legend.position = 'none')+xlab(NULL)+ylab(NULL)+theme(axis.text.x = element_blank(),
                                                                                      axis.text.y = element_blank()),
                  ps_pasture[[5]]+theme(legend.position = 'none')+xlab(NULL)+ylab(NULL)+theme(axis.text.x = element_blank(),
                                                                                       axis.text.y = element_blank()),
                  nrow =3,byrow = T)
    plot(pp)

# 12 Calculate temporal trend of mean occupancy probility, and plot map

## 12.1 Calculate temporal trend of mean occupancy probility

    library(unmarked)
    library(data.table)
    library(dplyr)
    library(doParallel)
    library(foreach)

    setwd('.../unmarked_out')
    all_sp<-Sys.glob('*.RData')
    source('.../functions/func_get_species_richness_map.R')
    source('.../functions/func_get_psi_change_map.R')
    map_raster<-raster(".../data/world_raster.tif")


    cl<- makeCluster(detectCores()-2)
    registerDoParallel(cl)

    list_psi_change<-foreach(
      s=1:length(all_sp),
      .packages = c('unmarked','data.table','dplyr'),
      .verbose=F 
    ) %dopar% get_psi_change(s)

    stopCluster(cl)

    mean(unlist(list_psi_change))
    sd(unlist(list_psi_change))/sqrt(425)

## 12.2 Calculate species richness map

    library(rgdal)
    library(raster)

    setwd('.../unmarked_mean_yearly_psi')
    all_sp<-Sys.glob('*.RData')

    ## calculate species richness map
    cl<- makeCluster(detcetCores()-2)
    registerDoParallel(cl)

    list_sp_raster<-foreach(
      s=1:length(all_sp),
      .packages = c('dplyr','raster','rgdal'),
      .verbose=F
    ) %dopar% get_map_SR(s,'.../data/split_raptor_shp_breeding')

    stopCluster(cl)
    SR_raster<-Reduce(function(x,y){sum(x,y,na.rm = T)},list_sp_raster)
    values(SR_sp_raster)[values(SR_sp_raster)==0]<-NA
    plot(SR_sp_raster)

    writeRaster(SR_sp_raster,
                '.../raptor_SR.tif',overwrite=T)

## 12.3 Calculate map of temporal trend of mean occupancy probility

    cl<- makeCluster(detcetCores()-2)
    registerDoParallel(cl)

    list_sp_raster<-foreach(
      s=1:length(all_sp),
      .packages = c('dplyr','raster','rgdal','data.table'),
      .verbose=F
    ) %dopar% get_psi_change_map(s,'.../data/split_raptor_shp_breeding')

    stopCluster(cl)

    mean_psi_change_sp_raster<-Reduce(function(x,y){sum(x,y,na.rm = T)},list_sp_raster)/SR_raster

    writeRaster(mean_psi_change_sp_raster,
                '.../mean_psi_change_sp_raster.tif',overwrite=T)

# 13 Calculate map of spatial patterns of mean spatiotemporal correlation and mean sensitivity

    library(dplyr)
    library(plyr)
    library(doParallel)
    library(foreach)
    library(data.table)
    library(raster)
    library(rgdal)

    ## read data
    coefs<-fread('.../coefs_MSSOMs.csv')
    senss<-fread('.../sensss_MSSDOMs.csv')

    diet_raptor<-read.csv('.../data/AVONET1_BirdLife.csv')
    crosswalk<-read.csv('.../data/BirdLife-BirdTree crosswalk.csv')
    tree_sys<-read.csv('.../data/BLIOCPhyloMasterTax.csv')

    ## combine data

    all_sp_coef_plot<-left_join(coefs,diet_raptor,by=c('species'='Species1')) %>%
      left_join(crosswalk,by=c('species'='Species1')) %>%
      left_join(tree_sys,by=c('Species3'='Scientific'))

    all_sp_coef_plot<-all_sp_coef_plot[all_sp_coef_plot$coefs_name!='(Intercept)',]
    all_sp_coef_plot$signf<-(all_sp_coef_plot$L_coef) * (all_sp_coef_plot$U_coef)
    all_sp_coef_plot$signf<-ifelse(all_sp_coef_plot$signf>0,'S','NS')
    all_sp_coef_plot$TipLabel[all_sp_coef_plot$signf>0] %>% unique() %>% length()


    all_sp_one_coef_plot1<-all_sp_coef_plot[c(grepl('LC_sett',all_sp_coef_plot$coefs_name)),]
    all_sp_one_sens_plot1<-senss[c(grepl('LC_sett',senss$sens_name)),]
    all_sp_one_coef_sens_plot1<-left_join(all_sp_one_coef_plot1, all_sp_one_sens_plot1,by='species')

    all_sp_one_coef_plot2<-all_sp_coef_plot[c(grepl('LC_agri',all_sp_coef_plot$coefs_name)),]
    all_sp_one_sens_plot2<-senss[c(grepl('LC_agri',senss$sens_name)),]
    all_sp_one_coef_sens_plot2<-left_join(all_sp_one_coef_plot2, all_sp_one_sens_plot2,by='species')

    all_sp_one_coef_plot3<-all_sp_coef_plot[c(grepl('pasture',all_sp_coef_plot$coefs_name)),]
    all_sp_one_sens_plot3<-senss[c(grepl('pasture',senss$sens_name)),]
    all_sp_one_coef_sens_plot3<-left_join(all_sp_one_coef_plot3, all_sp_one_sens_plot3,by='species')

    all_sp_one_coef_plot4<-all_sp_coef_plot[c(grepl('shdi',all_sp_coef_plot$coefs_name)),]
    all_sp_one_sens_plot4<-senss[c(grepl('shdi',senss$sens_name)),]
    all_sp_one_coef_sens_plot4<-left_join(all_sp_one_coef_plot4, all_sp_one_sens_plot4,by='species')


    map_raster<-raster(".../data/world_raster.tif")
    values(map_raster)<-0


    SR_sp_raster<-raster('.../raptor_SR.tif')


    ## for example: calculate map of raptor' sensitivity to settlement
    cl<- makeCluster(detectCores()-2)
    registerDoParallel(cl)

    list_sp_raster<-foreach(
      s=1:nrow(all_sp_one_coef_sens_plot1),
      .packages = c('dplyr','raster','rgdal'),
      .verbose=F
    ) %dopar% get_map_mean(s,.../data/split_raptor_shp_breeding,all_sp_one_coef_sens_plot1,'sens')

    stopCluster(cl)
    mean_sp_raster<-Reduce(function(x,y){sum(x,y,na.rm = T)},list_sp_raster)/SR_sp_raster

    writeRaster(mean_sp_raster,
                '.../sp_mean_sett_sens_map.tif',overwrite=T)
