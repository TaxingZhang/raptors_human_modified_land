func_occ_model<- function (s,output_path) {

  set.seed(1)
  target_file<-paste0(output_path,all_sp[s])

  if(file.exists(target_file)) {return()}
  
  load(all_sp[s])
  
  if(dim(occ_multi_frame@y)[1]<10){return()}
  
  year_site_cov<-names(occ_multi_frame@yearlySiteCovs)
  
  year_site_cov<-year_site_cov[!grepl('year|cell|site',year_site_cov)]
  
  site_cov<-year_site_cov[apply(occ_multi_frame@siteCovs %>% select(all_of(year_site_cov)), 2, sum)!=0]
  
  site_cov<-paste("scale(",site_cov,sep = '') %>%
    paste(')',sep = '')
  
  year_site_cov<-paste("scale(",year_site_cov,sep = '') %>%
    paste(')',sep = '')

  obs_cov<-names(occ_multi_frame@obsCovs)
  obs_cov<-obs_cov[!grepl('.obs_id|protocol_type',obs_cov)]

  obs_cov<-paste("scale(",obs_cov,sep = '') %>%
    paste(')',sep = '')
    
  fm<-colext(
    psiformula = as.formula(paste0('~',paste(site_cov,collapse  = '+'))),
    gammaformula = as.formula(paste0('~',paste(year_site_cov,collapse  = '+'))),
    epsilonformula = as.formula(paste0('~',paste(year_site_cov,collapse  = '+'))),
    pformula =  as.formula(paste0('~',paste(obs_cov,collapse  = '+'))),
    occ_multi_frame,
    control = list(trace=1,maxit=1e4))
  
  
  save(fm,file = target_file)
  rm(fm,occ_multi_frame)
  gc()
}


