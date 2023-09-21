get_MSSOMs_coef<-function(s) {
  load(all_sp[s])
  coefs_mean<-out$beta.samples %>% apply(2,mean) %>% as.data.frame()
  coefs_CI<-out$beta.samples %>% apply(2,function(x){quantile(x,c(0.025,0.975))}) %>% t() %>% as.data.frame()
  coefs_mean_CI<-cbind(coefs_CI,coefs_mean)
  
  colnames(coefs_mean_CI)<-c('L_coef','U_coef','mean_coef')
  coefs_mean_CI$coefs_name<-rownames(coefs_mean_CI)
  coefs_mean_CI$species<-all_sp[s] %>% gsub('.RData','',.)
  
  rownames(coefs_mean_CI)<-NULL
  return(coefs_mean_CI)
}