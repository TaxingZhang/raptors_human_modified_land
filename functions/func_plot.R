batch_plot<-function (land_type,parameter) {
  
  
  modN<-grep(land_type,c('LC_sett','LC_agri','pasture'))
  
  if( parameter =='sens') {modN=3+modN}
  
  ps<-list()
  
  # print(summary(mod_list[[1]]))
  ps[[1]]<- ggplot(all_mod_list[[modN]][[1]]$data$data, aes(x = Mass,y = mean_plot))+
    geom_point(size = 1.3)+geom_abline(slope = coefficients(all_mod_list[[modN]][[1]])[2],
                                       intercept = coefficients(all_mod_list[[modN]][[1]])[1],size=1)+
    theme_bw()
  # print(summary(all_mod_list[[modN]][[2]]))
  ps[[2]] <- ggplot(data=all_mod_list[[modN]][[2]]$data$data, aes(x = Range.Size,y = mean_plot))+
    geom_point(size = 1.3)+geom_abline(slope = coefficients(all_mod_list[[modN]][[2]])[2],
                                       intercept = coefficients(all_mod_list[[modN]][[2]])[1],size=1)+
    theme_bw()
  # 
  ps[[3]] <- ggplot(all_mod_list[[modN]][[5]]$data$data, 
                    aes(x = Trophic.Niche, 
                        y = mean_plot)) +
    stat_boxplot(geom = 'errorbar',width=0.2,cex=1)+
    geom_boxplot(notch = F, outlier.shape = NA)+
    geom_jitter(aes(color = Trophic.Niche,), size = 1.3, alpha=0.8,position = position_jitter(0.2)) +
    # stat_summary(fun.y = mean_plot, colour="black", geom="point", 
    #              shape=18, size=3, show_guide = FALSE) +
    scale_color_manual(values = col_list)+
    theme_bw()
  
 
  ps[[4]] <- ggplot(all_mod_list[[modN]][[3]]$data$data, 
                    aes(x = Migration,group = Migration, 
                        y = mean_plot)) +
    stat_boxplot(geom = 'errorbar',width=0.2,cex=1)+
    geom_boxplot(notch = F, outlier.shape = NA)+
    geom_jitter(aes(color=as.factor(Migration)),size = 1.3, alpha=0.8,position = position_jitter(0.2)) +
    # stat_summary(fun.y = mean, colour="black", geom="point", 
    #              shape=18, size=3, show_guide = FALSE) +
    scale_color_uchicago()+
    theme_bw()
  
  ps[[5]] <- ggplot(all_mod_list[[modN]][[4]]$data$data, 
                    aes(x = Habitat.Density, 
                        y = mean_plot,group = Habitat.Density)) +
    stat_boxplot(geom = 'errorbar',width=0.2,cex=1)+
    geom_boxplot(notch = F, outlier.shape = NA)+
    geom_jitter(aes(color=as.factor(Habitat.Density)),size = 1.3, alpha=0.8,position = position_jitter(0.2)) +
    # stat_summary(fun.y = mean, colour="black", geom="point", 
    #              shape=18, size=3, show_guide = FALSE) +
    scale_color_npg()+
    theme_bw()
  return(ps)
  # print(X)
}