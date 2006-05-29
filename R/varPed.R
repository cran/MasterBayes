"varPed" <-
function(x, gender=NULL, lag=c(0,0), relational=FALSE, lag_relational=c(0,0), restrict=NULL, keep=FALSE, USvar=NULL, merge=FALSE, ...){

      if(relational!=FALSE){
        if(relational!="OFFSPRING" & relational!="MATE"){stop("relational must either be 'OFFSPRING' or 'MATE'")}
      }
     
      sex<-with(parent.frame(), sex)               
      id<-with(parent.frame(), id)
      off_record<-with(parent.frame(), off_record)                
      data<-with(parent.frame(), data)                 # these are objects taken from the environment in which
      keepDam<-with(parent.frame(), keepDam)           # varPed is called - typically MCMCped
      keepSire<-with(parent.frame(), keepSire)  
      time_var<-with(parent.frame(), timevar)  
      namevar<-x
      x<-data[,x]                                      # gets variable(s)
      if(length(USvar)>0 & relational==FALSE){
        x[which(is.na(x)==TRUE)]<-USvar                # Fills missing values if specified in USvar
      }     
      hermaphrodite<-length(sex)==0                    # is this an hermaphroditic system
      sex_specific<-length(gender)==1                  # is the variable sex-specific

#####################################################################################################################
###################################  restricting variables ##########################################################
#####################################################################################################################

if(length(restrict)!=0){  

  PedDesMatrix<-list(Dam=list(id=NULL), Sire=list(id=NULL), Dam_restrict=list(id=NULL), Sire_restrict=list(id=NULL)) 

  not_after_off<-c(time_var<=time_var[off_record] | is.na(time_var))

  if(hermaphrodite==FALSE){
     PedDesMatrix$Dam$id<-unique(id[which(sex=="Female" & not_after_off)])
     PedDesMatrix$Sire$id<-unique(id[which(sex=="Male"  & not_after_off)])
     PedDesMatrix$Dam_restrict$id<-unique(id[which(sex=="Female" & not_after_off)])
     PedDesMatrix$Sire_restrict$id<-unique(id[which(sex=="Male" & not_after_off)])
  }else{
     PedDesMatrix$Dam$id<-unique(id[which(not_after_off)])
     PedDesMatrix$Sire$id<-unique(id[which(not_after_off)])
     PedDesMatrix$Dam_restrict$id<-unique(id[which(not_after_off)])
     PedDesMatrix$Sire_restrict$id<-unique(id[which(not_after_off)])
  }

  if(relational==FALSE){
    if("Female"%in%gender){ 
      PedDesMatrix$Dam_restrict$id<-unique(id[which(sex=="Female" & (x==restrict | is.na(x)==TRUE) & not_after_off)])
      if(keep==FALSE){
         PedDesMatrix$Dam$id<-PedDesMatrix$Dam_restrict$id
      }
    }
    if("Male"%in%gender){
      PedDesMatrix$Sire_restrict$id<-unique(id[which(sex=="Male" & (x==restrict | is.na(x)==TRUE) & not_after_off)])
      if(keep==FALSE){
         PedDesMatrix$Sire$id<-PedDesMatrix$Sire_restrict$id
      }
    }
    if(sex_specific==FALSE & hermaphrodite==FALSE){
      PedDesMatrix$Dam_restrict$id<-unique(id[which(sex=="Female" & (x==restrict | is.na(x)==TRUE) & not_after_off)])
      PedDesMatrix$Sire_restrict$id<-unique(id[which(sex=="Male" & (x==restrict | is.na(x)==TRUE) & not_after_off)])
      if(keep==FALSE){
         PedDesMatrix$Dam$id<-PedDesMatrix$Dam_restrict$id
         PedDesMatrix$Sire$id<-PedDesMatrix$Sire_restrict$id
      }
    }
    if(hermaphrodite==TRUE){
      PedDesMatrix$Dam_restrict$id<-unique(id[which(x==restrict & not_after_off)])
      PedDesMatrix$Sire_restrict$id<-unique(id[which(x==restrict & not_after_off)])
      if(keep==FALSE){
         PedDesMatrix$Dam$id<-PedDesMatrix$Dam_restrict$id
         PedDesMatrix$Sire$id<-PedDesMatrix$Sire_restrict$id
      }
    }
  }

  if(relational=="OFFSPRING"){
    if("Female"%in%gender){
      PedDesMatrix$Dam_restrict$id<-unique(id[which(sex=="Female" & ((x==x[off_record])==restrict | is.na(x)==TRUE) & not_after_off)])
      if(keep==FALSE){
         PedDesMatrix$Dam$id<-PedDesMatrix$Dam_restrict$id
      }
    }
    if("Male"%in%gender){
      PedDesMatrix$Sire_restrict$id<-unique(id[which(sex=="Male" & ((x==x[off_record])==restrict | is.na(x)==TRUE) & not_after_off)])
      if(keep==FALSE){
         PedDesMatrix$Sire$id<-PedDesMatrix$Sire_restrict$id
      }
    }
    if(sex_specific==FALSE & hermaphrodite==FALSE){
      PedDesMatrix$Dam_restrict$id<-unique(id[which(sex=="Female" & ((x==x[off_record])==restrict | is.na(x)==TRUE) & not_after_off)])
      PedDesMatrix$Sire_restrict$id<-unique(id[which(sex=="Male" & ((x==x[off_record])==restrict | is.na(x)==TRUE) & not_after_off)])
      if(keep==FALSE){
         PedDesMatrix$Dam$id<-PedDesMatrix$Dam_restrict$id
         PedDesMatrix$Sire$id<-PedDesMatrix$Sire_restrict$id
      }
    }
    if(hermaphrodite==TRUE){
      PedDesMatrix$Dam_restrict$id<-unique(id[which(((x==x[off_record])==restrict | is.na(x)==TRUE) & not_after_off)])
      PedDesMatrix$Sire_restrict$id<-unique(id[which(((x==x[off_record])==restrict | is.na(x)==TRUE) & not_after_off)])
      if(keep==FALSE){
         PedDesMatrix$Dam$id<-PedDesMatrix$Dam_restrict$id
         PedDesMatrix$Sire$id<-PedDesMatrix$Sire_restrict$id
      }
    }
  }
}
#####################################################################################################################
##########################################  true variables ##########################################################
#####################################################################################################################


if(length(restrict)==0){

    PedDesMatrix<-list(Dam=list(var_name=NULL, var_type=NULL, X=NULL, merge=FALSE), Sire=list(var_name=NULL, var_type=NULL, X=NULL, merge=FALSE), DamSire=list(var_name=NULL, var_type=NULL, X=NULL, merge=FALSE))    # output; design matrix+info

    if(length(dim(x)[1])==0){                          # gets class of variable(s)
      facnum<-class(x)
      if(facnum=="integer"){
        facnum<-"numeric"
      }
    }else{
      facnum<-apply(x, 2, class)
      facnum<-facnum[which(facnum=="integer")]<-"numeric"
    }        

    predict_ped=NULL

    off_time<-time_var[off_record]
    off_var<-as.matrix(x)[off_record,]

############### Covariates of fecundity, or distance from offspring ###############################

    if(relational!="MATE"){

      for(g in 1:c(2-length(gender))){

        if(sex_specific==FALSE & g==1){gender<-"Female"}
        if(sex_specific==FALSE & g==2){gender<-"Male"}
        if("Female"%in%gender){ 
          var_tmp<-subset(x, id%in%keepDam==TRUE)
          time_tmp<-subset(time_var, id%in%keepDam==TRUE)
          id_tmp<-subset(id, id%in%keepDam==TRUE)
        }
        if("Male"%in%gender){
          var_tmp<-subset(x, id%in%keepSire==TRUE)
          time_tmp<-subset(time_var, id%in%keepSire==TRUE)
          id_tmp<-subset(id, id%in%keepSire==TRUE)
        }

        time_for_P<-which((time_tmp>=(lag[1]+off_time) & time_tmp<=(lag[2]+off_time)) | is.na(time_tmp))

        # OFFSPRING RELATIONAL - NUMERIC

        if(relational=="OFFSPRING"  & "numeric"%in%facnum){
          var_tmp<-as.matrix(var_tmp)[time_for_P,]
          if(is.null(dim(var_tmp))){
            var_tmp<-t(as.matrix(var_tmp)) 
          }
          id_tmp<-id_tmp[time_for_P]
          dup_off_var<-rep(1, length(var_tmp[,1]))%*%t(off_var)
          predict_ped<-rowSums((var_tmp-dup_off_var)^2)^0.5
          if(length(USvar)>0){
            predict_ped[which(is.na(predict_ped)==T)]<-USvar          # Fills missing values if specified in USvar
          }
          namex<-namevar
          predict_ped<-tapply(predict_ped, id_tmp, mean, na.rm=T)
          id_tmp<-names(predict_ped)
          predict_ped<-matrix(predict_ped, length(predict_ped),1)
        }

# OFFSPRING RELATIONAL - FACTOR

        if(relational=="OFFSPRING" & "factor"%in%facnum){
          var_tmp<-var_tmp[time_for_P]
          id_tmp<-id_tmp[time_for_P]
          NAvec<-which(is.na(var_tmp)==TRUE)
          var_tmp<-var_tmp%in%off_var
          if(length(USvar)>0){
            var_tmp[NAvec]<-USvar      # Fills missing values if specified in USvar
          }else{
            var_tmp[NAvec]<-NA
          }    
          predict_ped<-tapply(var_tmp, id_tmp, mean)>0
          namex<-paste(namevar, c(TRUE, FALSE), sep=".")
          id_tmp<-names(predict_ped)
          predict_ped<-matrix(predict_ped, length(predict_ped),1)      
        }

# FECUNDITY - NUMERIC

        if(relational==FALSE & "numeric"%in%facnum){
          var_tmp<-var_tmp[time_for_P]
          id_tmp<-id_tmp[time_for_P]
          var_tmp<-tapply(var_tmp, id_tmp, mean, na.rm=T)
          id_tmp<-names(var_tmp)
          predict_ped<-var_tmp
          if(length(USvar)>0){
            predict_ped[which(is.na(predict_ped)==T)]<-USvar                # Fills missing values if specified in USvar
          }
          predict_ped<-matrix(predict_ped, length(predict_ped),1)
          namex<-namevar
        }

# FECUNDITY - FACTOR

        if(relational==FALSE & "factor"%in%facnum){
          var_tmp<-var_tmp[time_for_P]
          id_tmp<-id_tmp[time_for_P]
          NAvec<-which(is.na(var_tmp)==TRUE)
          if(length(USvar)>0){
            var_tmp[NAvec]<-USvar                # Fills missing values if specified in USvar
          }else{
            var_tmp[NAvec]<-levels(var_tmp)[1]
          }
          predict_ped<-as.matrix(model.matrix(~var_tmp)[,-1])
          namex<-paste(namevar, levels(var_tmp)[-1], sep=".")
          if(length(USvar)==0){
            predict_ped[NAvec,]<-NA
          }
          colnames(predict_ped)<-levels(var_tmp)[-1]
        }

        if(gender=="Female"){
          PedDesMatrix$Dam$X<-predict_ped
          if(sex_specific==TRUE){
            PedDesMatrix$Dam$var_name<-namex
          }else{
            PedDesMatrix$Dam$var_name<-paste(namex, "linked", sep=".")
          }
          PedDesMatrix$Dam$var_type<-facnum[1]
          if(merge==TRUE){
            PedDesMatrix$Dam$merge<-TRUE
          }
        }

        if(gender=="Male"){
          PedDesMatrix$Sire$X<-predict_ped
          if(sex_specific==TRUE){
            PedDesMatrix$Sire$var_name<-namex
          }else{
            PedDesMatrix$Sire$var_name<-paste(namex, "linked", sep=".")
          }
          PedDesMatrix$Sire$var_type<-facnum[1]
          if(merge==TRUE){
            PedDesMatrix$Sire$merge<-TRUE
          }
        }
      }
    }
############### Covariates of distance from mate ###############################

    if(relational=="MATE"){
  
      var_tmpF<-subset(x, id%in%keepDam==TRUE)
      time_tmpF<-subset(time_var, id%in%keepDam==TRUE)
      id_tmpF<-subset(id, id%in%keepDam==TRUE)
 
      var_tmpM<-subset(x, id%in%keepSire==TRUE)
      time_tmpM<-subset(time_var, id%in%keepSire==TRUE)
      id_tmpM<-subset(id, id%in%keepSire==TRUE)
 
      if(sex_specific==FALSE){gender="Female"}

      if(gender=="Female"){
        lagF<-lag
        lagM<-lag_relational
      }else{
        lagM<-lag
        lagF<-lag_relational
      }

      timePM<-c((time_tmpM>=(lagM[1]+off_time) & time_tmpM<=(lagM[2]+off_time)) | is.na(time_tmpM))
      timePF<-c((time_tmpF>=(lagF[1]+off_time) & time_tmpF<=(lagF[2]+off_time)) | is.na(time_tmpM))

      var_tmpM<-as.matrix(subset(var_tmpM, timePM))
      var_tmpF<-as.matrix(subset(var_tmpF, timePF))

      id_tmpM<-subset(id_tmpM, timePM)
      id_tmpF<-subset(id_tmpF, timePF)
 
      distmat<-matrix(0, nrow=length(id_tmpM), ncol=length(id_tmpF))
      id<-paste(rep(id_tmpF, each=length(id_tmpM)), rep(id_tmpM, length(id_tmpF)))

      if("numeric"%in%facnum){

        for(d in 1:length(var_tmpM[1,])){
          distmat<-distmat+(outer(c(var_tmpM[,d]), c(var_tmpF[,d]), "-")^2)
        }
       
        predict_ped<-tapply(c(distmat^0.5), id, mean, na.rm=T)
        predict_ped<-predict_ped[match(unique(id),names(predict_ped))]

        if(length(USvar)>0){
          predict_ped[which(is.na(predict_ped)==T)]<-USvar                # Fills missing values if specified in USvar
        }
        namex<-namevar
        predict_ped[1]<-0
      }
         

      if("factor"%in%facnum){
        distmat<-outer(c(var_tmpM), c(var_tmpF), "==")
        predict_ped<-tapply(distmat, id, mean, na.rm=T)>0
        predict_ped<-predict_ped[match(unique(id),names(predict_ped))]
        NAvec<-which(is.na(predict_ped)==TRUE)
        if(length(USvar)>0){
          predict_ped[NAvec]<-USvar                # Fills missing values if specified in USvar
        }else{
          predict_ped[NAvec]<-NA
        }
        namex<-paste(namevar, unique(predict_ped), sep=".")
      }

      predict_ped<-matrix(predict_ped, length(predict_ped),1)
      PedDesMatrix$DamSire$X<-predict_ped
      PedDesMatrix$DamSire$var_name<-namex
      PedDesMatrix$DamSire$var_type<-facnum[1]
      if(merge==TRUE){
        PedDesMatrix$DamSire$merge<-TRUE
      }
    }
  }
PedDesMatrix
}
