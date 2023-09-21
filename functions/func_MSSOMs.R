func_occ_model<- function (s,output_path) {
  
  set.seed(1)
  target_file<-paste0(output_path,all_sp[s])
  
  if(file.exists(target_file)) {return()}
  
  load(all_sp[s])
  
  occ_covs<-names(data.list$occ.covs)
  occ_covs<-occ_covs[occ_covs!='site']
  
  occ_covs<-paste("scale(",occ_covs,sep = '') %>%
    paste(')',sep = '')
  
  det_covs<-names(data.list$det.covs)
  det_covs<-paste("scale(",det_covs,sep = '') %>%
    paste(')',sep = '')
  
  out<-stPGOcc(occ.formula = as.formula(paste0('~',paste(occ_covs,collapse  = '+'))), 
               det.formula = as.formula(paste0('~',paste(det_covs,collapse  = '+'))), 
               data = data.list,
               n.batch = 1000, 
               batch.length = 25,
               ar1 = T,
               n.neighbors = 5,
               n.report = 100,
               n.thin = 25,
               n.chains = 1)
  save(out,file = target_file)
  rm(out)
  gc()
}








