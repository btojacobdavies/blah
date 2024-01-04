require(rucrdtw)

#Make library of reference curves for each phenological category
sim_data<-data.frame("wintering"=rep(0,366),
                     "passage"=rep(0,366),
                     "breeding"=rep(0,366))
sim_data$wintering[1:100]<-1
sim_data$wintering[267:366]<-1
sim_data$passage[75:125]<-1
sim_data$passage[250:300]<-1
sim_data$breeding[125:250]<-1
sim_data$wintering_passage<-sim_data$wintering
sim_data$wintering_passage[1:50]<-0.5
sim_data$wintering_passage[317:366]<-0.5

sim_data$breeding_passage<-sim_data$breeding
sim_data$breeding_passage[180:200]<-0.5
sim_mat<-as.matrix(sim_data)

plot(sim_data$wintering,xlab="Julian day",ylab="Occurrence probability",main="Wintering")
plot(sim_data$passage,xlab="Julian day",ylab="Occurrence probability",main="Passage")
plot(sim_data$breeding,xlab="Julian day",ylab="Occurrence probability",main="Breeding")
plot(sim_data$wintering_passage,xlab="Julian day",ylab="Occurrence probability",main="Wintering passage")
plot(sim_data$breeding_passage,xlab="Julian day",ylab="Occurrence probability",main="Breeding passage")

table(occ_df$scheme)/nrow(occ_df)

#Plot all records
png("H:/JGD H drive backup 160320/Africa migration phenology/Plots/Data availability summary/all lists.png",7,7,res=100,units="in")
ggplot()+
  geom_sf(data=world_shp_s_sf, fill="#FFFFCC", colour="wheat3")+
  geom_sf(data=af_15deg_grid_sf, fill=NA, colour="gray50")+
  coord_sf(xlim=extent(af_15deg_r)[1:2],ylim=extent(af_15deg_r)[3:4])+
  geom_point(data=occ_df,aes(x=lon,y=lat),alpha=0.3)
dev.off()

table(is.na(occ_df$cell_15deg))

temp_sp_phen_stats_list<-list()
for(j in 1:length(target_spp)){
  print(j)
  print(target_spp[j])
  
  temp_occ_df<-occ_df[(is.na(occ_df[,target_spp_underscore[j]])==F) & 
                        occ_df$list_complete==1, #NB Complete lists only
                      c(target_spp_underscore[j],"jul","duration_hrs","start_hr","cell_15deg","lon","lat")]
  temp_occ_df$obs_bin<-temp_occ_df[,target_spp_underscore[j]]
  
  if(sum(temp_occ_df$obs_bin)==0) next #Subsetting to complete lists can leave zero records for some spp e.g. ring ouzel
  
  #Which 15-deg cells have more than 100 occurrences for this species?
  temp_cellfreq_df<-data.frame(table(temp_occ_df$cell_15deg[temp_occ_df[,target_spp_underscore[j]]==1]))
  temp_cellfreq_df[order(temp_cellfreq_df$Freq),]
  data_rich_cells<-as.numeric(as.character(temp_cellfreq_df[temp_cellfreq_df$Freq>100,]$Var1))
  
  if(length(data_rich_cells)==0) next
  
  temp_sp_phen_stats<-data.frame(matrix(ncol=10,nrow=length(data_rich_cells)))
  names(temp_sp_phen_stats)<-c("species","cell","phen_status","dist",
                               "spr_pass_peak","aut_pass_peak",
                               "spr_dep_0.1","spr_dep_0.9",
                               "aut_arr_0.1","aut_arr_0.9")
  temp_sp_phen_stats$species<-target_spp[j]
  temp_sp_phen_stats$cell<-data_rich_cells
  
  for(i in 1:length(data_rich_cells)){
    m1<-gam(obs_bin~s(jul,bs="cc")+s(duration_hrs)+s(start_hr), #bs=="cc" gives a *CYCLIC* cubic regression spline: ends match (don't do this for hour, it'll assume that 2000 runs immediately into 0500 if there are no data from between those times)
            family="binomial",data=temp_occ_df[temp_occ_df$cell_15deg==data_rich_cells[i],])
    
    newdata<-data.frame(jul=1:366,
                        duration_hrs=rep(median(temp_occ_df$duration_hrs,na.rm=T),366),
                        start_hr=rep(median(temp_occ_df$start_hr,na.rm=T),366))
    pred<-predict(m1,newdata,type="iterms",se=T) 
    pred_pn<-plogis(pred$fit[,1]+coef(m1)[1])
    pred_uci<-plogis(pred$fit[,1]+(1.96*pred$se.fit[,1])+coef(m1)[1])
    pred_lci<-plogis(pred$fit[,1]-(1.96*pred$se.fit[,1])+coef(m1)[1])
    plot_df<-data.frame(jul=newdata$jul,est=pred_pn,lci=pred_lci,uci=pred_uci)
    
    #Fill in phen status and dtw distance
    temp_sp_phen_stats$dist[i]<-ucrdtw_mv(sim_mat,plot_df$est,0.05)$distance
    temp_phen_ind<-ucrdtw_mv(sim_mat,plot_df$est,0.05)$location
    temp_sp_phen_stats$phen_status[i]<-names(sim_data)[temp_phen_ind]
    
    #Identifying spring and autumn peaks for pure passage squares
    #If phen_status==passage and if dtw_dist<10
    if(temp_sp_phen_stats$phen_status[i]=="passage" & temp_sp_phen_stats$dist[i]<12){
      plot(est~jul,plot_df,
           main=paste0(data_rich_cells[i],": passage; dist = ",round(temp_sp_phen_stats$dist[i],2)))
      spring_peak_y<-max(plot_df$est[plot_df$jul<200])
      spring_peak_x<-plot_df$jul[which(plot_df$jul<200 & plot_df$est==spring_peak_y)]
      temp_sp_phen_stats$spr_pass_peak[i]<-spring_peak_x
      abline(v=spring_peak_x,col="green",lwd=2)
      autumn_peak_y<-max(plot_df$est[plot_df$jul>200])
      autumn_peak_x<-plot_df$jul[which(plot_df$jul>200 & plot_df$est==autumn_peak_y)]
      temp_sp_phen_stats$aut_pass_peak[i]<-autumn_peak_x
      abline(v=autumn_peak_x,col="brown",lwd=2)
    }
    
    #Identifying spring departure from and autumn arrival to wintering grounds
    #If phen_status==wintering and if dtw_dist<10
    if(temp_sp_phen_stats$phen_status[i]=="wintering" & temp_sp_phen_stats$dist[i]<12){
      plot(est~jul,plot_df,
           main=paste0(data_rich_cells[i],": wintering; dist = ",round(temp_sp_phen_stats$dist[i],2)))
      #Spring first
      min_occ_spring<-min(plot_df$est[plot_df$jul<200])
      max_occ_spring<-max(plot_df$est[plot_df$jul<200])
      range_occ_spring<-max_occ_spring-min_occ_spring
      occ_spring_0.1<-max_occ_spring-(0.1*range_occ_spring)
      spring_depart_0.1<-max(plot_df$jul[which(plot_df$est>occ_spring_0.1 & plot_df$jul<200)])
      temp_sp_phen_stats$spr_dep_0.1[i]<-spring_depart_0.1
      abline(v=spring_depart_0.1,col="light green",lwd=2)
      occ_spring_0.9<-max_occ_spring-(0.9*range_occ_spring)
      spring_depart_0.9<-min(plot_df$jul[which(plot_df$est<occ_spring_0.9 & plot_df$jul<200)])
      temp_sp_phen_stats$spr_dep_0.9[i]<-spring_depart_0.9
      abline(v=spring_depart_0.9,col="dark green",lwd=2)
      
      #Now autumn
      min_occ_autumn<-min(plot_df$est[plot_df$jul>200])
      max_occ_autumn<-max(plot_df$est[plot_df$jul>200])
      range_occ_autumn<-max_occ_autumn-min_occ_autumn
      occ_autumn_0.1<-min_occ_autumn+(0.1*range_occ_autumn)
      autumn_arrival_0.1<-max(plot_df$jul[which(plot_df$est<occ_autumn_0.1 & plot_df$jul>200)])
      temp_sp_phen_stats$aut_arr_0.1[i]<-autumn_arrival_0.1
      abline(v=autumn_arrival_0.1,col="orange1",lwd=2)
      occ_autumn_0.9<-min_occ_autumn+(0.9*range_occ_autumn)
      autumn_arrival_0.9<-min(plot_df$jul[which(plot_df$est>occ_autumn_0.9 & plot_df$jul>200)])
      temp_sp_phen_stats$aut_arr_0.9[i]<-autumn_arrival_0.9
      abline(v=autumn_arrival_0.9,col="orangered1",lwd=2)
    }
  }
  
  # median(temp_sp_phen_stats$spr_dep_0.9-temp_sp_phen_stats$spr_dep_0.1,na.rm=T)
  # median(temp_sp_phen_stats$aut_arr_0.9-temp_sp_phen_stats$aut_arr_0.1,na.rm=T)
  
  temp_sp_phen_stats_list[[j]]<-temp_sp_phen_stats
}


phen_stats_df<-do.call(rbind,temp_sp_phen_stats_list)
phen_stats_df$spr_dep_dur<-phen_stats_df$spr_dep_0.9-phen_stats_df$spr_dep_0.1
phen_stats_df$aut_arr_dur<-phen_stats_df$aut_arr_0.9-phen_stats_df$aut_arr_0.1
hist(phen_stats_df$spr_dep_dur,breaks=20,xlab="Duration of spring departure from wintering grounds",main="")
hist(phen_stats_df$aut_arr_dur,breaks=20,xlab="Duration of autumn arrival to wintering grounds",main="")

phen_stats_agg_df<-aggregate(phen_stats_df[,5:ncol(phen_stats_df)],list(phen_stats_df$species),median,na.rm=T)
names(phen_stats_agg_df)[1]<-"species"
hist(phen_stats_agg_df$spr_dep_dur,breaks=20,xlab="Duration of spring departure from wintering grounds",main="")
hist(phen_stats_agg_df$aut_arr_dur,breaks=20,xlab="Duration of autumn arrival to wintering grounds",main="")

#Add in classification from Jenni's paper
phen_stats_agg_df$euro_cluster<-character(nrow(phen_stats_agg_df))
phen_stats_agg_df$euro_cluster[phen_stats_agg_df$species %in% c("Phylloscopus trochilus",
                                                                "Phoenicurus phoenicurus",
                                                                "Anthus trivialis",
                                                                "Oenanthe oenanthe",
                                                                "Cuculus canorus",
                                                                "Curruca communis",
                                                                "Sylvia borin",
                                                                "Saxicola rubetra",
                                                                "Phylloscopus sibilatrix",
                                                                "Ficedula hypoleuca")]<-"C1"
phen_stats_agg_df$euro_cluster[phen_stats_agg_df$species %in% c("Delichon urbicum",
                                                                "Upupa epops",
                                                                "Phylloscopus collybita",
                                                                "Jynx torquilla",
                                                                "Curruca curruca",
                                                                "Turdus torquatus",
                                                                "Sylvia atricapilla",
                                                                "Apus apus",
                                                                "Hirundo rustica",
                                                                "Motacilla flava",
                                                                "Riparia riparia",
                                                                "Acrocephalus schoenobaenus")]<-"C2"
phen_stats_agg_df$euro_cluster[phen_stats_agg_df$species %in% c("Ficedula albicollis",
                                                                "Caprimulgus europaeus",
                                                                "Streptopelia turtur",
                                                                "Acrocephalus palustris",
                                                                "Emberiza hortulana",
                                                                "Oriolus oriolus",
                                                                "Lanius collurio")]<-"C3"

#Add in difference between mean passage peak and departure from/arrival to breeding ground
phen_stats_agg_df$spr_dur<-phen_stats_agg_df$spr_pass_peak-phen_stats_agg_df$spr_dep_0.1
phen_stats_agg_df$aut_dur<-phen_stats_agg_df$aut_arr_0.9-phen_stats_agg_df$aut_pass_peak
# save(phen_stats_agg_df,file="H:/JGD H drive backup 160320/Africa migration phenology/Output data/phen_stats_agg_df.RData")
# load("H:/JGD H drive backup 160320/Africa migration phenology/Output data/phen_stats_agg_df.RData")

phen_stats_agg_df[order(phen_stats_agg_df$spr_dur),]

#What are the relationships among phenology statistics?
pairs(phen_stats_agg_df[,c("spr_pass_peak","spr_dep_0.1","spr_dep_0.9","spr_dep_dur")])
pairs(phen_stats_agg_df[,c("aut_pass_peak","aut_arr_0.1","aut_arr_0.9","aut_arr_dur")])

head(phen_stats_agg_df)

plot(spr_dep_dur~spr_pass_peak,phen_stats_agg_df)
plot(aut_arr_dur~aut_pass_peak,phen_stats_agg_df)
plot(aut_pass_peak~spr_pass_peak,phen_stats_agg_df)
plot(spr_dep_dur~aut_arr_dur,phen_stats_agg_df) #Maybe not as impressive as it looks. These are not estimated independently. A curve that is steep on one side may well be steep on the other.
m1<-lm(spr_dep_dur~aut_arr_dur,phen_stats_agg_df);summary(m1)
plot(spr_dur~aut_dur,phen_stats_agg_df) 
plot(spr_dur~spr_dep_0.1,phen_stats_agg_df) #Not sure this means anything interesting
plot(aut_dur~aut_arr_0.9,phen_stats_agg_df) #Not sure this means anything interesting

sd(phen_stats_agg_df$spr_pass_peak,na.rm=T)
sd(phen_stats_agg_df$spr_dep_0.1,na.rm=T)


#How do phenology statistics differ among different European migrant types (early/slow, late/fast)
boxplot(spr_pass_peak~euro_cluster,phen_stats_agg_df,ylab="Julian day",main="Spring passage peak")
boxplot(aut_pass_peak~euro_cluster,phen_stats_agg_df,ylab="Julian day",main="Autumn passage peak")
boxplot(spr_dur~euro_cluster,phen_stats_agg_df,ylab="Days",main="Duration of spring passage")
boxplot(aut_dur~euro_cluster,phen_stats_agg_df,ylab="Days",main="Duration of autumn passage")
#Next plots are potential core of paper.
boxplot(spr_dep_0.1~euro_cluster,phen_stats_agg_df,ylab="Julian day",main="Spring departure from wintering grounds (10%)")
boxplot(spr_dep_0.9~euro_cluster,phen_stats_agg_df,ylab="Julian day",main="Spring departure from wintering grounds (90%)")
boxplot(aut_arr_0.1~euro_cluster,phen_stats_agg_df,ylab="Julian day",main="Autumn arrival to wintering grounds (10%)")
boxplot(aut_arr_0.9~euro_cluster,phen_stats_agg_df,ylab="Julian day",main="Autumn arrival to wintering grounds (90%)")
boxplot(spr_dep_dur~euro_cluster,phen_stats_agg_df,ylab="Days",main="Duration of spring departure from wintering grounds")
boxplot(aut_arr_dur~euro_cluster,phen_stats_agg_df,ylab="Days",main="Duration of autumn arrival to wintering grounds")

phen_stats_agg_df$euro_cluster_col<-"black"
phen_stats_agg_df$euro_cluster_col[phen_stats_agg_df$euro_cluster=="C1"]<-"red"
phen_stats_agg_df$euro_cluster_col[phen_stats_agg_df$euro_cluster=="C2"]<-"green"
phen_stats_agg_df$euro_cluster_col[phen_stats_agg_df$euro_cluster=="C3"]<-"blue"
plot(spr_dep_dur~spr_dep_0.1,phen_stats_agg_df,col=phen_stats_agg_df$euro_cluster_col,pch=20,cex=2,
     xlab="Spring departure from wintering grounds (10%)",ylab="Duration of spring departure from wintering grounds")
legend("topright",legend=c("N/A","C1","C2","C3"),pch=c(20,20,20,20),col=c("black","red","green","blue"),pt.cex=c(2,2,2,2))
plot(aut_arr_dur~aut_arr_0.1,phen_stats_agg_df,col=phen_stats_agg_df$euro_cluster_col,pch=20,cex=2,
     xlab="Autumn arrival to wintering grounds (10%)",ylab="Duration of autumn arrival to wintering grounds")
legend("topright",legend=c("N/A","C1","C2","C3"),pch=c(20,20,20,20),col=c("black","red","green","blue"),pt.cex=c(2,2,2,2))

phen_stats_agg_df[order(phen_stats_agg_df$spr_dep_0.9),]

#Add in stats from Jenni's paper
Europe_uk_spr_filter_by_0.1<-readRDS("H:/JGD H drive backup 160320/Africa migration phenology/Jenni European migration phenology/Europe_uk_spr_filter_by_0.1.rds")
Europe_uk_spr_filter_by_0.1$spr_dur<-Europe_uk_spr_filter_by_0.1$end_spr-Europe_uk_spr_filter_by_0.1$start_spr
euro_phen_summary_df<-aggregate(Europe_uk_spr_filter_by_0.1$spr_dur,list(Europe_uk_spr_filter_by_0.1$engname),median)
names(euro_phen_summary_df)<-c("species","median_spr_dur")
euro_phen_summary_df$median_start_spr<-aggregate(Europe_uk_spr_filter_by_0.1$start_spr,list(Europe_uk_spr_filter_by_0.1$engname),median)$x
euro_phen_summary_df$median_end_spr<-aggregate(Europe_uk_spr_filter_by_0.1$end_spr,list(Europe_uk_spr_filter_by_0.1$engname),median)$x

EURING_conv_tab<-read.csv("H:/JGD H drive backup 160320/EFSA AI/Cloned MMT code/ico02_efsa_migration_mapping_tool/data/euring_lookup_tables/species_list_allocated_euring_codes.csv")
head(EURING_conv_tab)

euro_phen_summary_df$EURING_species<-euro_phen_summary_df$species
euro_phen_summary_df$EURING_species[euro_phen_summary_df$species=="Blackcap"]<-"Eurasian Blackcap"
euro_phen_summary_df$EURING_species[euro_phen_summary_df$species=="Chiffchaff"]<-"Common Chiffchaff"
euro_phen_summary_df$EURING_species[euro_phen_summary_df$species=="Cuckoo"]<-"Common Cuckoo"
euro_phen_summary_df$EURING_species[euro_phen_summary_df$species=="Golden Oriole"]<-"Eurasian Golden Oriole"
euro_phen_summary_df$EURING_species[euro_phen_summary_df$species=="Hoopoe"]<-"Eurasian Hoopoe"
euro_phen_summary_df$EURING_species[euro_phen_summary_df$species=="House Martin"]<-"Common House Martin"
euro_phen_summary_df$EURING_species[euro_phen_summary_df$species=="Nightingale"]<-"Common Nightingale"
euro_phen_summary_df$EURING_species[euro_phen_summary_df$species=="Nightjar"]<-"European Nightjar"
euro_phen_summary_df$EURING_species[euro_phen_summary_df$species=="Pied Flycatcher"]<-"European Pied Flycatcher"
euro_phen_summary_df$EURING_species[euro_phen_summary_df$species=="Redstart"]<-"Common Redstart"
euro_phen_summary_df$EURING_species[euro_phen_summary_df$species=="Swallow"]<-"Barn Swallow"
euro_phen_summary_df$EURING_species[euro_phen_summary_df$species=="Swift"]<-"Common Swift"
euro_phen_summary_df$EURING_species[euro_phen_summary_df$species=="Turtle Dove"]<-"European Turtle Dove"
euro_phen_summary_df$EURING_species[euro_phen_summary_df$species=="Wheatear"]<-"Northern Wheatear"
euro_phen_summary_df$EURING_species[euro_phen_summary_df$species=="Whitethroat"]<-"Common Whitethroat"
euro_phen_summary_df$EURING_species[euro_phen_summary_df$species=="Wryneck"]<-"Eurasian Wryneck"
euro_phen_summary_df$EURING_species[euro_phen_summary_df$species=="Yellow Wagtail"]<-"Western Yellow Wagtail"

euro_phen_summary_df$sci_name<-""
for(i in 1:nrow(euro_phen_summary_df)){
  euro_phen_summary_df$sci_name[i]<-EURING_conv_tab[EURING_conv_tab$common.name==euro_phen_summary_df$EURING_species[i],]$species.sci.name.ioc
}
euro_phen_summary_df$sci_name[euro_phen_summary_df$sci_name=="Sylvia curruca"]<-"Curruca curruca"
euro_phen_summary_df$sci_name[euro_phen_summary_df$sci_name=="Sylvia communis"]<-"Curruca communis"

phen_stats_agg_df$euro_spr_end<-phen_stats_agg_df$euro_spr_start<-phen_stats_agg_df$euro_spr_dur<-NA

for(i in 1:nrow(euro_phen_summary_df)){
  phen_stats_agg_df[phen_stats_agg_df$species==euro_phen_summary_df$sci_name[i],]$euro_spr_dur<-euro_phen_summary_df$median_spr_dur[i]  
  phen_stats_agg_df[phen_stats_agg_df$species==euro_phen_summary_df$sci_name[i],]$euro_spr_start<-euro_phen_summary_df$median_start_spr[i]  
  phen_stats_agg_df[phen_stats_agg_df$species==euro_phen_summary_df$sci_name[i],]$euro_spr_end<-euro_phen_summary_df$median_end_spr[i]  
}
plot(euro_spr_dur~spr_dep_dur,phen_stats_agg_df,xlab="Median duration of departure from African wintering grounds",ylab="Median duration of arrival to European breeding grounds")
m1<-lm(euro_spr_dur~spr_dep_dur,phen_stats_agg_df);summary(m1);abline(m1)

plot(euro_spr_start~spr_dep_0.1,phen_stats_agg_df)
plot(euro_spr_start~spr_dep_0.9,phen_stats_agg_df)
plot(euro_spr_end~spr_dep_0.1,phen_stats_agg_df)
plot(euro_spr_end~spr_dep_0.9,phen_stats_agg_df)


write.csv(phen_stats_agg_df[-c(8,9)],"H:/JGD H drive backup 160320/Africa migration phenology/Output data/phen_stats_agg_df_191023.csv",row.names=F)
?write.csv
