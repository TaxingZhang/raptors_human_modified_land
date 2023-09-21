batch_pgls<-function(land_type,parameter) {
  
  niche<-c("Vertivore","Scavenger","Aquatic predator","Omnivore","Invertivore")
  mod_list<-list()
  all_sp_one_coef_plot<-all_sp_coef_plot[c(grepl(land_type,all_sp_coef_plot$coefs_name)),]
  all_sp_one_sens_plot<-senss[c(grepl(land_type,senss$sens_name)),]
  all_sp_one_coef_sens_plot<-left_join(all_sp_one_coef_plot, all_sp_one_sens_plot,by='species') %>%
    ddply('TipLabel',summarise,
          mean_coef=mean(mean_coef),mean_sens=mean(mean_sens),
          L_coef=mean(L_coef),U_coef=mean(U_coef),L_sens=mean(L_sens),U_sens=mean(U_sens),
          Mass=mean(Mass),Range.Size=mean(Range.Size),Migration=unique(Migration),
          Primary.Lifestyle=unique(Primary.Lifestyle),Habitat.Density=unique(Habitat.Density),
          Trophic.Niche=unique(Trophic.Niche),Trophic.Level=unique(Trophic.Level),Order=unique(Order1))
  all_sp_one_coef_sens_plot$signf_coef<-ifelse((all_sp_one_coef_sens_plot$L_coef * all_sp_one_coef_sens_plot$U_coef)>0,'S','NS')
  all_sp_one_coef_sens_plot$signf_sens<-ifelse((all_sp_one_coef_sens_plot$L_sens * all_sp_one_coef_sens_plot$U_sens)>0,'S','NS')
  
  # mean(all_sp_one_coef_sens_plot$mean_sens)
  # sd(all_sp_one_coef_sens_plot$mean_sens)/sqrt(425)
  
  # sum(all_sp_one_coef_sens_plot$signf_sens=='S')
  # sum(all_sp_one_coef_sens_plot$signf_sens=='S')/425
  
  # sum(all_sp_one_coef_sens_plot$signf_sens=='S' & all_sp_one_coef_sens_plot$mean_sens>0)
  # sum(all_sp_one_coef_sens_plot$signf_sens=='S' & all_sp_one_coef_sens_plot$mean_sens>0)/425
  # 
  # sum(all_sp_one_coef_sens_plot$signf_sens=='S' & all_sp_one_coef_sens_plot$mean_sens<0)
  # sum(all_sp_one_coef_sens_plot$signf_sens=='S' & all_sp_one_coef_sens_plot$mean_sens<0)/425
  
  #
  # sum(all_sp_one_coef_sens_plot$signf_coef=='S')
  # sum(all_sp_one_coef_sens_plot$signf_coef=='S')/425
  #
  # sum(all_sp_one_coef_sens_plot$signf_coef=='S' & all_sp_one_coef_sens_plot$mean_coef>0)
  # sum(all_sp_one_coef_sens_plot$signf_coef=='S' & all_sp_one_coef_sens_plot$mean_coef>0)/425
  #
  # sum(all_sp_one_coef_sens_plot$signf_coef=='S' & all_sp_one_coef_sens_plot$mean_coef<0)
  # sum(all_sp_one_coef_sens_plot$signf_coef=='S' & all_sp_one_coef_sens_plot$mean_coef<0)/425
  #
  # mean(all_sp_one_coef_sens_plot$mean_coef)
  # sd(all_sp_one_coef_sens_plot$mean_coef)/sqrt(425)
  
  
  if(parameter=='coef'){
    all_sp_one_coef_sens_plot$mean_plot<-all_sp_one_coef_sens_plot[,'mean_coef']
  } else {
    all_sp_one_coef_sens_plot$mean_plot<-all_sp_one_coef_sens_plot[,'mean_sens']
  }
  
  all_sp_one_coef_sens_plot$Mass<-(all_sp_one_coef_sens_plot$Mass)/1000
  all_sp_one_coef_sens_plot$Range.Size<-(all_sp_one_coef_sens_plot$Range.Size)/1000000
  
  
  pdata <- comparative.data(phy = tree_dropped_tip, data = all_sp_one_coef_sens_plot,
                            names.col = TipLabel, vcv = TRUE,
                            na.omit = FALSE, warn.dropped = TRUE)
  
  
  mod_list[[1]]<- pgls(mean_plot ~ Mass,data = pdata, lambda = "ML")
  mod_list[[2]] <- pgls(mean_plot ~ Range.Size,data = pdata, lambda = "ML")
  mod_list[[3]] <- pgls(mean_plot ~ as.factor(Habitat.Density),data = pdata, lambda = "ML")
  mod_list[[4]] <- pgls(mean_plot ~ as.factor(Migration),data = pdata, lambda = "ML")
  
  all_sp_one_coef_sens_plot<-all_sp_one_coef_sens_plot[all_sp_one_coef_sens_plot$Trophic.Niche!='Frugivore',]
  pdata <- comparative.data(phy =  tree_dropped_tip, data = all_sp_one_coef_sens_plot,
                            names.col = TipLabel, vcv = TRUE,
                            na.omit = FALSE, warn.dropped = TRUE)
  
  mod_list[[5]] <- pgls(mean_plot~ as.factor(Trophic.Niche),data = pdata, lambda = "ML")
  
  names(mod_list)<-c('Mass','Range.size','Habitat','Migration','Trophic Niche')
  return(mod_list)
}