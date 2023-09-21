library(auk)
library(data.table)
library(dplyr)
library(dggridR)
library(doParallel)
library(foreach)
library(purrr)
library(plyr)
library(stringr)
library(tidyr)
library(unmarked)

#####
func_get_pasture<-function (year) {
  # year=2001
  load(paste0('.../hexs/hex_pasture/',year,'.RData'),.GlobalEnv)
  cell_pasture<-list()
  for (c in 1:nrow(all_cells)) {
    # c=1
    if(all_cells$year[c]!=year) next
    cell_pasture[[c]]<-list_year_NEW_land[[all_cells$cell[c]]]
  }
  cell_pasture<-Reduce(function(x,y) {bind_rows(x,y)},cell_pasture)
  return(cell_pasture)
}
#####
func_get_landcover<-function (year) {
  # year=2001
  load(paste0('.../hexs/hex_landcover/',year,'.RData'),.GlobalEnv)
  cell_landcover<-list()
  for (c in 1:nrow(all_cells)) {
    if(all_cells$year[c]!=year) next
    cell_landcover[[c]]<-list_year_landcover[[all_cells$cell[c]]]
  }
  cell_landcover<-Reduce(function(x,y) {bind_rows(x,y)},cell_landcover)
  return(cell_landcover)
}
#####
func_link_cov<-function(input_path,output_path) {
  
  setwd(input_path)
  all_sp<-Sys.glob('*.csv')
  dggs=dgconstruct(area = 50*50)
  
  for(s in 1:length(all_sp)) {

    sp<-all_sp[s] %>% gsub('.csv','',.)
    print(sp)
    target_file<-paste0(output_path,'/',sp,'.csv')
    if(file.exists(target_file)){
      next
    }
    

    ebird_zf<-fread(all_sp[s])
    
    if(length(unique(ebird_zf$cell))<10){
      next
    }
    
    ##extract occupancy covariates
    all_cells<-select(ebird_zf,year,cell) %>% unique()
    
    cl<- makeCluster(30)
    registerDoParallel(cl)
    
    print("extracting pasture")
    cell_pasture<- foreach(
      year=2001:2020, 
      .combine=rbind,
      .packages = c('dplyr'),
      .verbose=F 
    ) %dopar% func_get_pasture(year)
    
    stopCluster(cl)
    
    colnames(cell_pasture)[1]<-'pasture'
    ########
    cl<- makeCluster(30)      
    registerDoParallel(cl)
    
    print("extracting landcover")
    cell_landcover<- foreach(
      year=2001:2020,
      .combine=dplyr::bind_rows,
      .packages = c('dplyr'),
      .verbose=F
    ) %dopar% func_get_landcover(year)
    
    stopCluster(cl)
    
    ebird_zf_cov<-left_join(ebird_zf,all_year_pop,by=c('cell','year')) %>%
      left_join(all_year_npp,by=c('cell','year')) %>%
      left_join(all_year_shdi,by=c('cell','year')) %>%
      left_join(cell_landcover,by=c('cell','year')) %>%
      left_join(cell_pasture,by=c('cell','year'))
    
    rm(ebird_zf)
    gc()
    
    landcover_ID<-colnames(ebird_zf_cov) %>% as.numeric()
    landcover_ID<-landcover_ID[!is.na(landcover_ID)]
    
    agri_ID<-landcover_ID[landcover_ID<=40] %>% as.character() %>% as.vector()
    sett_ID<-landcover_ID[landcover_ID==190] %>% as.character() %>% as.vector()
    
    
    ebird_zf_cov$LC_agri<-ebird_zf_cov[,..agri_ID] %>% apply(1,function(x){sum(x,na.rm=T)})
    ebird_zf_cov$LC_sett<-ebird_zf_cov[,..sett_ID] %>% apply(1,function(x){sum(x,na.rm=T)})
    
    landcover_ID<-as.character(landcover_ID)
    
    ebird_zf_cov<-ebird_zf_cov[,!..landcover_ID]
    
    gc()
    
    occ_cov<-c('pop','npp','shdi')
    
    if('LC_sett' %in% colnames(ebird_zf_cov)){
      occ_cov<-c(occ_cov,'LC_sett')
    }
    
    if('pasture' %in% colnames(ebird_zf_cov)){
      occ_cov<-c(occ_cov,'pasture')
    }
    
    if('LC_agri' %in% colnames(ebird_zf_cov)){
      occ_cov<-c(occ_cov,'LC_agri')
    }
    
    ebird_zf_cov<-drop_na(ebird_zf_cov,all_of(occ_cov))
    
    ## extract observer score
    
    # get unique observers per sei
    dataSeiScore <- distinct(
      ebird_zf_cov, sampling_event_identifier,
      observer_id
    ) %>%
      # make list column of observers
      mutate(observers = str_split(observer_id, ",")) %>%
      unnest(cols = c(observers)) %>%
      # add numeric observer id
      mutate(numObserver = str_extract(observers, "\\d+")) %>%
      # now get distinct sei and observer id numeric
      distinct(sampling_event_identifier, numObserver)
    
    dataSeiScore <- left_join(dataSeiScore, mean_observer_score,
                              by = "numObserver") %>%
      # get max expertise score per sei
      ddply('sampling_event_identifier',summarise,expertise = max(mean_score))
    
    ebird_zf_cov<-left_join(ebird_zf_cov,dataSeiScore,by=c('sampling_event_identifier')) %>%
      filter(!is.na(expertise)) %>%setDT()
    
    fwrite(ebird_zf_cov,target_file)
    rm(ebird_zf_cov)
    gc()
  }
}




