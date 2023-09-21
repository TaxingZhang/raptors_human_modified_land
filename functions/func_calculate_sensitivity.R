calculate_sens<-function (s,input_coefs_path) {
  
  coef_conf_file<-paste0(input_coefs_path,gsub('.RData','.csv',all_sp[s]))
  
  coef_conf<-fread(coef_conf_file)
  
  col_coef_conf<-coef_conf[grepl('col',coef_conf$V1),] %>% setDF()
  ext_coef_conf<-coef_conf[grepl('ext',coef_conf$V1),] %>% setDF()
  
  load(all_sp[s])
  
  env<-fm@data@yearlySiteCovs
  
  env$pop<-scale(env$pop)
  env$npp<-scale(env$npp)
  env$shdi<-scale(env$shdi)
  env$pasture<-scale(env$pasture)
  env$LC_agri<-scale(env$LC_agri)
  env$LC_sett<-scale(env$LC_sett)
  

  env <-ddply(env,'cell',summarise,
              pop=mean(pop),
              npp=mean(npp),
              shdi=mean(shdi),
              pasture=mean(pasture),
              LC_agri=mean(LC_agri),
              LC_sett=mean(LC_sett))
  
  
  beta_gamma_npp<-col_coef_conf[col_coef_conf$V1 %>% grepl('npp',.),'Est']
  beta_gamma_pop<-col_coef_conf[col_coef_conf$V1 %>% grepl('pop',.),'Est']
  beta_gamma_shdi<-col_coef_conf[col_coef_conf$V1 %>% grepl('shdi',.),'Est']
  beta_gamma_pasture<-col_coef_conf[col_coef_conf$V1 %>% grepl('pasture',.),'Est']
  beta_gamma_LC_agri<-col_coef_conf[col_coef_conf$V1 %>% grepl('LC_agri',.),'Est']
  beta_gamma_LC_sett<-col_coef_conf[col_coef_conf$V1 %>% grepl('LC_sett',.),'Est']
  int_gamma<-col_coef_conf[col_coef_conf$V1 %>% grepl('Int',.),'Est']
  
  beta_epsilon_npp<-ext_coef_conf[ext_coef_conf$V1 %>% grepl('npp',.),'Est']
  beta_epsilon_pop<-ext_coef_conf[ext_coef_conf$V1 %>% grepl('pop',.),'Est']
  beta_epsilon_shdi<-ext_coef_conf[ext_coef_conf$V1 %>% grepl('shdi',.),'Est']
  beta_epsilon_pasture<-ext_coef_conf[ext_coef_conf$V1 %>% grepl('pasture',.),'Est']
  beta_epsilon_LC_agri<-ext_coef_conf[ext_coef_conf$V1 %>% grepl('LC_agri',.),'Est']
  beta_epsilon_LC_sett<-ext_coef_conf[ext_coef_conf$V1 %>% grepl('LC_sett',.),'Est']
  int_epsilon<-ext_coef_conf[ext_coef_conf$V1 %>% grepl('Int',.),'Est']
  
  
  attach(env)
  
  X_gamma<-beta_gamma_npp*npp+
    beta_gamma_pop*pop+
    beta_gamma_shdi*shdi+
    beta_gamma_pasture*pasture+
    beta_gamma_LC_agri*LC_agri+
    beta_gamma_LC_sett*LC_sett+
    int_gamma
  
  X_epsilon<-beta_epsilon_npp*npp+
    beta_epsilon_pop*pop+
    beta_epsilon_shdi*shdi+
    beta_epsilon_pasture*pasture+
    beta_epsilon_LC_agri*LC_agri+
    beta_epsilon_LC_sett*LC_sett+
    int_epsilon
  
  
  gamma<-exp(X_gamma)/(1+exp(X_gamma))
  epsilon<-exp(X_epsilon)/(1+exp(X_epsilon))
  
  psi_a<-gamma/(epsilon+gamma)
  
  dpsi_a.dgamma<-epsilon/(gamma+epsilon)^2
  dpsi_a.depsilon<-(gamma/(gamma+epsilon)^2)*-1

  Xs<-c('pop','npp','shdi','LC_sett','LC_agri','pasture')
  
  list_dpsi.dXs<-list()
  for (v in 1:6) {
    beta_gamma_v<-col_coef_conf[col_coef_conf$V1 %>% grepl(Xs[v],.),'Est']
    beta_epsilon_v<-ext_coef_conf[col_coef_conf$V1 %>% grepl(Xs[v],.),'Est']
    
    dgamma.dXs<- beta_gamma_v*exp(X_gamma)/(exp(2*(X_gamma))+2*exp(X_gamma)+1)
    
    depsilon.dXs<- beta_epsilon_v*exp(X_epsilon)/(exp(2*(X_epsilon))+2*exp(X_epsilon)+1)
    
    list_dpsi.dXs[[v]]<-(dpsi_a.dgamma*dgamma.dXs+dpsi_a.depsilon*depsilon.dXs)/psi_a
    
  }
  dpsi.dXs<-bind_cols(list_dpsi.dXs)
  
  colnames(dpsi.dXs)<-paste0('dpsi.d',Xs[1:6])
  dpsi.dXs <- dpsi.dXs[!(is.infinite(rowSums(dpsi.dXs))|is.na(rowSums(dpsi.dXs))),]
  
  
  dpsi.dXs<-data.frame(mean=apply(dpsi.dXs, 2, function(x){mean(x,na.rm=T)}),
                       L=apply(dpsi.dXs, 2, function(x){quantile(x,0.025,na.rm=T)}),
                       U=apply(dpsi.dXs, 2, function(x){quantile(x,0.975,na.rm=T)}))
  
  dpsi.dXs$sens_name<-rownames(dpsi.dXs)
  detach()
  
  colnames(dpsi.dXs)[1:3]<-c('mean_sens','L_sens','U_sens')
  
  dpsi.dXs$species<-all_sp[s] %>% gsub('.RData','',.)

  return(dpsi.dXs)

}