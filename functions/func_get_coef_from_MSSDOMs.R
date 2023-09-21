get_MSSDOMs_coef<-function(s,output_path) { 
  load(all_sp[sp])
  
  target_file<-paste0(output_path,'/',
                      gsub('.RData','.csv',all_sp[sp]))
  if(file.exists(target_file)) {return()}
  
  coefs<-coef(fm) %>% as.data.frame()
  
  psi_conf<-confint(fm,type='psi')
  col_conf<-confint(fm,type='col')
  ext_conf<-confint(fm,type='ext')
  p_conf<-confint(fm,type='det')
  
  confs<-rbind(psi_conf,col_conf,ext_conf,p_conf)
  coefs_confs<-data.frame(coefs,confs)
  colnames(coefs_confs)<-c('Est','L','U')
  write.csv(coefs_confs,target_file)
}