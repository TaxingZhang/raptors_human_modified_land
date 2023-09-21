get_psi_change_map<-function (s,shp_path) {
  sp_shp<-readOGR(paste0(shp_path,'/',
                         all_sp[s] %>% 
                           gsub(' ','_',.) %>% 
                           gsub('.RData','.shp',.)))
  
  sp_raster<-mask(map_raster,sp_shp)
  values(sp_raster)[!is.na(values(sp_raster))]<-list_psi_change[[s]]
  return(sp_raster)
}