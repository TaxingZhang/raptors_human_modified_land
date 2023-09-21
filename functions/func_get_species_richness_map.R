get_map_SR<-function (s,shp_path) {
  sp_shp<-readOGR(paste0(shp_path,'/',
                         all_sp[s] %>% 
                           gsub(' ','_',.) %>% 
                           gsub('.RData','.shp',.)))
  
  sp_raster<-mask(map_raster,sp_shp)
  values(sp_raste)[!is.na(values(sp_raster))]<-1
  values(sp_raster)[is.na(values(sp_raster))]<-0
  return(sp_raster)
}