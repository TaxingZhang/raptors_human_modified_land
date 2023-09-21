make_data<-function (s,output_path) {
  
  sp=all_sp[s] %>% gsub('.csv','',.)
  target_file<-paste0(output_path,sp,'.RData')
  if(file.exists(target_file)){
    return()
  }

  if(length(unique(ebird_zf_cov$cell))<10){
    return()
  }
  
  obs_cov<-c("duration_minutes", 
             'min_obs_started',
             'julian_date',
             'number_observers',
             'expertise',
             'effort_distance_km')
  
  ebird_zf_cov$site<-paste0(ebird_zf_cov$cell,"_",ebird_zf_cov$year)
  
  ebird_zf_cov<-ebird_zf_cov %>% na.omit()
  
  occ_cov<-c("pop","npp", "shdi","pasture","LC_agri","LC_sett")
  
  occ_cov<-occ_cov[occ_cov %in% colnames(ebird_zf_cov)]
  
  ebird_zf_cov<-ebird_zf_cov[site!="",] %>% setDF()

  
  occ_wide <- format_unmarked_occu(ebird_zf_cov,
                                   site_id = "site", 
                                   response = "species_observed",
                                   site_covs = c('cell','year',occ_cov),
                                   obs_covs = obs_cov)#
  
  
  occ_wide<-as.data.table(occ_wide)
  
  
  all_cells<-unique(occ_wide[,'cell'])
  all_cells$cell_lon<-apply(all_cells[,'cell'], 1, function(x){dgSEQNUM_to_GEO(dggs,x)[[1]]})
  all_cells$cell_lat<-apply(all_cells[,'cell'], 1, function(x){dgSEQNUM_to_GEO(dggs,x)[[2]]})
  
  occ_wide<-left_join(occ_wide,all_cells,by=c('cell'))
  
  list_year_y<-list()
  list_year_obs.1<-list()
  list_year_obs.2<-list()
  list_year_obs.3<-list()
  list_year_obs.4<-list()
  list_year_obs.5<-list()
  list_year_obs.6<-list()
  nrep<-sum(grepl('y.\\d',colnames(occ_wide)))

  for (r in 1:nrep) {
    # r=1
    print(r)
    year_all_cell<-all_cells
    list_year_occ_wide<-list()
    for (y in unique(occ_wide$year)) {
      # y=2001
      need_colname<-c('cell',
                      paste0('y.',r),
                      c(occ_cov),
                      paste0(obs_cov[1],'.',r),
                      paste0(obs_cov[2],'.',r),
                      paste0(obs_cov[3],'.',r),
                      paste0(obs_cov[4],'.',r),
                      paste0(obs_cov[5],'.',r),
                      paste0(obs_cov[6],'.',r))
      
      
      list_year_occ_wide[[y-2000]]<-occ_wide[year==y,..need_colname]
      
      colnames(list_year_occ_wide[[y-2000]])[2:ncol(list_year_occ_wide[[y-2000]])]<-
        paste0(colnames(list_year_occ_wide[[y-2000]])[2:ncol(list_year_occ_wide[[y-2000]])],'_',y)
      
      year_all_cell<-left_join(year_all_cell,list_year_occ_wide[[y-2000]],by='cell')
      print(y)
    }
    year_all_cell<-setDF(year_all_cell)
    list_year_y[[r]]<-year_all_cell[,grep(paste0('y.',r,'_'),colnames(year_all_cell))] %>% as.matrix()
    list_year_obs.1[[r]]<-year_all_cell[,grep(paste0(obs_cov[1],'.',r,'_'),colnames(year_all_cell))] %>% as.matrix()
    list_year_obs.2[[r]]<-year_all_cell[,grep(paste0(obs_cov[2],'.',r,'_'),colnames(year_all_cell))] %>% as.matrix()
    list_year_obs.3[[r]]<-year_all_cell[,grep(paste0(obs_cov[3],'.',r,'_'),colnames(year_all_cell))] %>% as.matrix()
    list_year_obs.4[[r]]<-year_all_cell[,grep(paste0(obs_cov[4],'.',r,'_'),colnames(year_all_cell))] %>% as.matrix()
    list_year_obs.5[[r]]<-year_all_cell[,grep(paste0(obs_cov[5],'.',r,'_'),colnames(year_all_cell))] %>% as.matrix()
    list_year_obs.6[[r]]<-year_all_cell[,grep(paste0(obs_cov[6],'.',r,'_'),colnames(year_all_cell))] %>% as.matrix()

    if(r==1){
      occ.covs<-list()
      for (c in 1:length(occ_cov)) {
        # c=1
        occ.covs[[c]]<-year_all_cell[,grep(occ_cov[c],colnames(year_all_cell))] %>% as.matrix()
        occ.covs[[c]][is.na(occ.covs[[c]])]<-0
        names(occ.covs)[c]<-occ_cov[c]
      }
    }
  }
  
  array_year_y<-array(unlist(list_year_y),dim=c(nrow(all_cells),length(unique(occ_wide$year)),sum(grepl('y.\\d',colnames(occ_wide)))))
  array_year_obs.1<-array(unlist(list_year_obs.1),dim=c(nrow(all_cells),length(unique(occ_wide$year)),sum(grepl('y.\\d',colnames(occ_wide)))))
  array_year_obs.2<-array(unlist(list_year_obs.2),dim=c(nrow(all_cells),length(unique(occ_wide$year)),sum(grepl('y.\\d',colnames(occ_wide)))))
  array_year_obs.3<-array(unlist(list_year_obs.3),dim=c(nrow(all_cells),length(unique(occ_wide$year)),sum(grepl('y.\\d',colnames(occ_wide)))))
  array_year_obs.4<-array(unlist(list_year_obs.4),dim=c(nrow(all_cells),length(unique(occ_wide$year)),sum(grepl('y.\\d',colnames(occ_wide)))))
  array_year_obs.5<-array(unlist(list_year_obs.5),dim=c(nrow(all_cells),length(unique(occ_wide$year)),sum(grepl('y.\\d',colnames(occ_wide)))))
  array_year_obs.6<-array(unlist(list_year_obs.6),dim=c(nrow(all_cells),length(unique(occ_wide$year)),sum(grepl('y.\\d',colnames(occ_wide)))))

  # dm= the time spent observing birds for that checklist
  # tos= the time of the day that observations were started
  # od = julian_date
  # no = number observers
  # ex = expertise,observer score
  # ed = effort_distance_km
  data.list<-list(y=array_year_y,
                    occ.covs=occ.covs,
                    det.covs=list(tos=array_year_obs.1,
                                  dm=array_year_obs.2,
                                  od=array_year_obs.3,
                                  no=array_year_obs.4,
                                  ex=array_year_obs.5,
                                  ed=array_year_obs.6),
                    coords=all_cells[,2:3] %>% as.matrix())

  save(data.list,file = target_file)
  rm(data.list)
  rm(ebird_zf_cov)
  gc()
}





