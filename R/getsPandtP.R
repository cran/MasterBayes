"getsPandtP"<-function(sP, tP, PdP, GdP, X.list, nbeta, unique_id, checkP, ...){

  if(sP$estP==TRUE & sP$estG==FALSE & length(sP$G)==0){CERVUS<-TRUE}else{CERVUS<-FALSE} 
  # use the cervus approximation for genotyping error

  if(is.null(sP$id)){
    sP$id<-X.list$id
  }

  if(sP$estA == TRUE | sP$estG==TRUE | CERVUS==TRUE){ 
    if(length(sP$A)==0){           
      if(length(sP$G)!=0){        
        sP$A<-extractA(sP$G)
      }else{                      
        sP$A<-extractA(GdP$G)
      }
    }
  }

  No.E<-length(unique(GdP$categories))*length(GdP$G)*GdP$perlocus+((1-GdP$perlocus)*length(unique(GdP$categories)))

  if(is.null(sP$E1)){       
    if(is.null(GdP$categories)){    
      sP$E1<-0.005
    }else{
      sP$E1<-rep(0.005,No.E)
      if(GdP$perlocus==FALSE){
        names(sP$E1)<-unique(GdP$categories)
      }else{
        names(sP$E1)<-paste(unique(GdP$categories), rep(names(GdP$G), each=length(unique(GdP$categories))), sep=".")
      }
    }
  }else{
    if(sP$estG==TRUE | any(sP$E1==0)){
      sP$E1[which(sP$E1==0)]<-1e-5
    }
  }
  
  if(is.null(sP$E2)){       
    if(is.null(GdP$categories)){    
      sP$E2<-0.005
    }else{
      sP$E2<-rep(0.005,No.E)
      if(GdP$perlocus==FALSE){
        names(sP$E2)<-unique(GdP$categories)
      }else{
        names(sP$E2)<-paste(unique(GdP$categories), rep(names(GdP$G), each=length(unique(GdP$categories))), sep=".")
      }
    }
  }else{
    if(sP$estG==TRUE | any(sP$E2==0)){
      sP$E2[which(sP$E2==0)]<-1e-5
    }
  }



  # if sP$G does not exist obtain from the first genotype of each individual in GdP$G

  if(is.null(sP$G)){
    sPGgiven<-FALSE
  }else{
    sPGgiven<-TRUE
  }

    if(sPGgiven == FALSE){ 
      sP$G<-lapply(GdP$G, function(x){x[-duplicated(GdP$id)==FALSE]}) 
    }

    ped<-matrix(NA, length(sP$id), 3)

    ped[,1]<-sP$id

    if(is.null(sP$dam)==FALSE){
      ped[,2]<-match(sP$dam, unique_id)
    }
    if(is.null(sP$sire)==FALSE){
      ped[,3]<-match(sP$sire, unique_id)
    }

    if(sP$estUSdam==TRUE | sP$estUSsire==TRUE){
      if(sP$estUSsire=="USdam"){
        MLENus<-MLE.popsize(X.list, USdam=PdP$USdam, USsire="USdam", ped=ped)
      }else{
        MLENus<-MLE.popsize(X.list, USdam=PdP$USdam, USsire=PdP$USsire, ped=ped)
      }
    }

    if(length(PdP$USdam)==1 & PdP$USdam[1]==FALSE){
      nusd<-0
    }else{
      nusd<-length(unique(PdP$USdam))
    }

    if(length(PdP$USsire)==1 & PdP$USsire[1]==FALSE){
      sP$estUSsire<-FALSE
      nuss<-0
    }else{
      nuss<-length(unique(PdP$USsire))
    }

  if(sP$estUSdam==TRUE & is.null(sP$USdam)){
    sP$USdam<-MLENus$nUS[1:nusd]
  }

  if((sP$estUSsire==TRUE | sP$estUSsire=="USdam") & is.null(sP$USsire)){
    sP$USsire<-MLENus$nUS[(nusd+1):(nusd+nuss)]
  }

  if(sP$estP==TRUE){  
   ped<-MLE.ped(X.list, ped=ped, USdam=PdP$USdam, nUSdam=sP$USdam, USsire=PdP$USsire, nUSsire=sP$USsire, checkP=checkP)
  }

  sP$dam<-ped[,2]
  sP$sire<-ped[,3]
  
  if(sP$estbeta==TRUE){
    MLEestimates<-MLE.beta(X.list, ped=ped, beta=sP$beta, nUSdam=sP$USdam, nUSsire=sP$USsire)
  }

  if(sP$estbeta==TRUE){
    if(is.null(sP$beta)){
      sP$beta<-MLEestimates$beta
    }else{
      if(length(sP$beta)>1){
        sP$beta<-sP$beta[X.list$beta_map]
      }
    }
  }


  if(is.null(GdP$G)==FALSE & CERVUS==FALSE){
    if(sPGgiven==TRUE){
      if(is.null(PdP$timevar)){
        time_born=NULL
      }else{
        time_born = PdP$timevar[which(PdP$offspring==1)][match(ped[,1], PdP$id)] 
      }
      if(legalG(sP$G, sP$A, ped, time_born=time_born, marker.type=GdP$marker.type)$valid=="FALSE"){
        warning("sP$G does not have postive probability given possible starting pedigree")
        stop()
      }
    }else{   
      if(is.null(GdP$G)==FALSE){   
        sP$G<-legalG(sP$G, sP$A, ped, marker.type=GdP$marker.type)$G
      }  
    }
  }

  # get a legal configuration for sP$G 

############################################ get tuning parameters ###############################################

      if(sP$estE1){
        if(is.null(tP$E1)){  
          tP$E1<-chol(diag(No.E)*0.00003)
        }else{
          tP$E1<-chol(diag(No.E)*0.00003*tP$E1)
        }
      }

      if(sP$estE2==TRUE){
        if(is.null(tP$E2)){
          tP$E2<-chol(diag(No.E)*0.00003)
        }else{
          tP$E2<-chol(diag(No.E)*0.00003*tP$E2)
        }
      }


    if(sP$estbeta==TRUE){
      if(is.null(tP$beta)){  
        tP$beta<-rep(sqrt(10),length(sP$beta))
      }else{
        if(length(tP$beta)>1){
          tP$beta<-tP$beta[X.list$beta_map]*rep(sqrt(10),length(sP$beta))
        }else{
          tP$beta<-tP$beta*rep(sqrt(10),length(sP$beta))
        }
      }
      tP$beta<-sqrt(tP$beta%*%t(tP$beta))*MLEestimates$C
      tP$beta<-chol(tP$beta)
    }

    if(sP$estUSdam==TRUE | sP$estUSsire==TRUE){
       if(sum(diag(MLENus$C)<0)>0){
         warning("Hessian not positive-definite for MLE.popsize")
       }
    }

    if(sP$estUSdam){
      if(is.null(tP$USdam)){
        tP$USdam<-rep(10,nusd)
      }else{
        tP$USdam<-tP$USdam*rep(10, nusd)
      } 
      tP$USdam<-abs(tP$USdam*diag(MLENus$C)[1:nusd])
    }

     if(sP$estUSsire==TRUE | sP$estUSsire=="USdam"){
       if(is.null(tP$USsire)){
         tP$USsire<-rep(10, nuss)
       }else{
         tP$USsire<-tP$USsire*rep(10, nuss)
       }
       tP$USsire<-abs(tP$USsire*diag(MLENus$C)[nusd+(1:nuss)])
    }


list(sP=sP, tP=tP)
}
