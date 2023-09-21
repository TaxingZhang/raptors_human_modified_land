get_map_mean<-function (s,shp_path,data,x) {
  sp_shp<-readOGR(paste0(shp_path,'/',data$species[s] %>% 
                           gsub(' ','_',.),'.shp'))
  
  sp_raster<-mask(map_raster,sp_shp)
  if(x=='coef') {
    values(sp_raster)[!is.na(values(sp_raster))]<-data$mean_coef[s]
  } else {
    values(sp_raster)[!is.na(values(sp_raster))]<-data$mean_sens[s]
  }
  return(sp_raster)
}