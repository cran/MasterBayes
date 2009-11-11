insertPed<-function(ped, founders=NULL){
   mmothers<-na.omit(ped[,2][which(ped[,2]%in%ped[,1]==FALSE)])
   mfathers<-na.omit(ped[,3][which(ped[,3]%in%ped[,1]==FALSE)])
   if(is.null(founders)==FALSE){
     founders<-na.omit(founders[which(founders%in%ped[,1]==FALSE)])
   }
   mparents<-unique(c(mmothers, mfathers, founders))
   nped<-rbind(cbind(mparents, rep(NA, length(mparents)), rep(NA, length(mparents))),ped)
   colnames(nped)<-colnames(ped)
   nped
}

