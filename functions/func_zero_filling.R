func_zero_filling<-function (s,f_ebd_path,f_sampling_path,expert_map_path,output_path) {
  
  if(!dir.exists(output_path)) {
    dir.create(output_path)
  }
  ## test if target file exist
  target_file<-paste0(output_path,'/',raptor[s],'.csv')
  if (file.exists(target_file)) {return()}
  
  ## input path
  f_ebd<-paste0(f_ebd_path,'/filted_raptor_ebd_',raptor[s],'.txt')
  f_sampling<-read_sampling(f_sampling_path,'/filtered_raptor_sampling_',raptor[s],'.txt')
  
  ## Some species have zero rows data after being filtered
  if("try-error" %in% class(try(read_ebd(f_ebd),silent = T))) {return()}
  
  ## input expert map
  
  ## The layer was obtained from Birdlife International
  ## we split the layer to single species file in Arcgis
  ## and removed the non-breeding range.
  ## To avoid taxonomic errors
  ## we analyzed only species with consistent Latin names 
  ## in the taxonomic systems of Ebird and Birdlife
  expert_map<- try(readOGR(paste0(expert_map_path,'/',gsub(' ','_',raptor[s]),'.shp')),
                   silent = T)
  if("try-error" %in% class(expert_map)) {
    return()
  }

  
  filter_sampling<-read_sampling(f_sampling)
  filter_ebd<-read_ebd(f_ebd)
  
  ##Make sure checklists in the two files match
  filter_ebd<-filter_ebd[filter_ebd$checklist_id %in% filter_sampling$checklist_id,]
  
  ##get zero filling data
  zf<-auk_zerofill(filter_ebd,filter_sampling,collapse = T) %>%
    mutate(
      # convert X to NA
      observation_count = if_else(observation_count == "X",
                                  NA_character_, observation_count),
      observation_count = as.integer(observation_count),
      # effort_distance_km to 0 for non-travelling counts
      effort_distance_km = if_else(protocol_type != "Traveling",
                                   0, effort_distance_km),
      # split date into year and day of year
      year = year(observation_date),
      day_of_year = yday(observation_date),
      # occupancy modeling requires an integer response
      species_observed = as.integer(species_observed),
      # Assign each record to a hexagonal grid of responses
      cell = dgGEO_to_SEQNUM(dggs,longitude, latitude)$seqnum
    ) %>%# Convert the data into the form required by the model
    filter_repeat_visits(min_obs =5,
                         max_obs = 999999,
                         annual_closure = T,
                         date_var = "observation_date",
                         site_vars = c("cell"))
  
  site_ID<-unique(select(zf,cell))
  
  site_ID<-site_ID %>%
    cbind(dgSEQNUM_to_GEO(dggs,site_ID$cell))
  
  cell_we_need<-site_ID[,2:3] %>%
    as.data.frame() %>%
    SpatialPoints(proj4string=expert_map@proj4string) %>%
    over(expert_map) %>%
    select(PRESENCE)
  
  cell_we_need<-!is.na(point_need)
  
  site_ID<-site_ID[point_need,]
  
  zf<-zf[zf$cell %in% site_ID$cell,]
  
  ##Some species have very little data 
  ##and are left with an empty table 
  ##after a series of filtering
  if(nrow(zf)==0) return()
  
  fread(zf, target_file)
}


