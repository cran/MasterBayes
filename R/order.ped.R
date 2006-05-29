"order.ped"<-function(ped, time_born=NULL, ...){

  reorder<-ped[,1]
  
  if(length(time_born)!=0){
    reorder<-ped[,1][order(time_born)]
  }else{
    reorder<-ped[,1][order((is.na(ped[,2])==FALSE & is.na(ped[,3])==FALSE))]
  }

  for(i in 1:length(ped[,1])){
    off_i<-match(ped[,1][i], reorder)
    dam_i<-match(ped[,2][i], reorder)
    sire_i<-match(ped[,3][i], reorder)
    if((is.na(dam_i)==FALSE & dam_i>off_i) | (is.na(sire_i)==FALSE & sire_i>off_i)){
      max_par<-max(dam_i, sire_i, na.rm=T)
      reorder<-append(reorder, ped[,1][i], max_par)
      reorder<-reorder[-off_i]
    }   
  }
ped[match(reorder, ped[,1]),]
}
