MLE.ped<-function(X.list, ped=NULL, USdam=FALSE, nUSdam=NULL, USsire=FALSE, nUSsire=NULL, threshold=0, ...){

         if(is.null(ped)){
           ped<-matrix(NA, length(X.list$id), 3)
           ped[,1]<-as.character(X.list$id)
         }

         if(FALSE%in%is.na(ped[,2])){
           getd<-FALSE
         }else{
           getd<-TRUE
         }
         if(FALSE%in%is.na(ped[,3])){
           gets<-FALSE
         }else{
           gets<-TRUE
         }

        if(length(USdam)==1){
          if(USdam==TRUE){
            USdam<-rep(1, length(X.list$X))
            betaDcat<-1
          }else{
            USdam<-NULL
            betaDcat<-NULL
          }
        }else{
          betaDcat<-unique(USdam)
        }

        if(length(USsire)==1){
          if(USsire==TRUE){
            USsire<-rep(1, length(X.list$X))
            betaScat<-1
          }else{
            USsire<-NULL
            betaScat<-NULL
          }
        }else{
          betaScat<-unique(USsire)
        }

        nbetaD<-length(betaDcat)
        nbetaS<-length(betaScat)


         ndam<-unlist(lapply(X.list$X, function(x){length(x$restdam.id)}))
         nsire<-unlist(lapply(X.list$X, function(x){length(x$restsire.id)}))

         if(length(X.list$X[[1]]$G)==0){
           warning("X.list$X is missing genetic likelihoods")
           stop()
         }


         for(i in 1:length(X.list$X)){

           d_cat<-match(USdam[i], betaDcat)
           s_cat<-match(USsire[i], betaScat)

           pos_in_id<-as.numeric(names(X.list$X)[i])
 
           if(gets==TRUE | getd==TRUE){

             X<-X.list$X[[i]]$G[1:(ndam[i]*nsire[i])]
  
             nsampD<-ndam[i]-(nbetaD>0)
             nsampS<-nsire[i]-(nbetaS>0)

             sampD<-1:(nsampD*nsire[i])
             sampS<-1:(ndam[i]*nsire[i])
             if(length(USsire)>0){
               sampS<-sampS[-c(1:ndam[i])*nsire[i]]
             }
 
             DandS<-intersect(sampS, sampD)
             DnotS<-setdiff(sampD, sampS)
             SnotD<-setdiff(sampS, sampD)
             notDS<-if(nbetaD>0 & nbetaS>0){ndam[i]*nsire[i]} 
 
             X[DnotS]<-X[DnotS]*if(nbetaS>0){nUSsire[s_cat]}else{0}
             X[SnotD]<-X[SnotD]*if(nbetaD>0){nUSdam[d_cat]}else{0}
             X[notDS]<-X[notDS]*if(nbetaD>0 & nbetaS>0){nUSdam[d_cat]*nUSsire[s_cat]}else{0}
 
             X<-t(matrix(X, nsire[i], ndam[i]))
           }
     
           if(getd==TRUE & gets==TRUE){      # if neither parents have starting parameteristion
             MLpar<-which.max(t(X))[1]
             MLdam<-ceiling(MLpar/nsire[i])
             MLsire<-MLpar-((ceiling(MLpar/nsire[i])-1)*nsire[i])
             if((t(X)[MLpar]/sum(X))<threshold){
               MLdam<-NA
               MLsire<-NA
             }
           }else{
             if(gets==FALSE){
               MLsire<-match(ped[,3][pos_in_id], X.list$id[X.list$X[[i]]$restsire.id])
               if(is.na(MLsire)==TRUE & nbetaS>0){MLsire<-nsire[i]}  
               if(getd==TRUE){   
                 MLdam<-which.max(X[,MLsire])
                 if(MLdam==length(X.list$X[[i]]$restdam.id) & nbetaD>0){MLdam<-NA}  
                 if((X[,MLsire][MLdam]/sum(X[,MLsire]))<threshold){MLdam<-NA}  
               }  
             }
             if(getd==FALSE){
               MLdam<-match(ped[,2][pos_in_id], X.list$id[X.list$X[[i]]$restdam.id])
               if(is.na(MLdam)==TRUE & nbetaD>0){MLdam<-ndam[i]}
               if(gets==TRUE){        
                 MLsire<-which.max(X[,MLdam])
                 if(MLsire==length(X.list$X[[i]]$restsire.id) & nbetaS>0){MLsire<-NA}
                 if((X[,MLdam][MLsire]/sum(X[,MLdam]))<threshold){MLsire<-NA}  

               }
             }
           }
           if((MLdam==ndam[i] & nbetaD>0) | is.na(MLdam)==TRUE){
             ped[,2][pos_in_id]<-NA
           }else{
             ped[,2][pos_in_id]<-as.character(X.list$id[X.list$X[[i]]$restdam.id[MLdam]])
           }
           if((MLsire==nsire[i] & nbetaS>0) | is.na(MLsire)==TRUE){
             ped[,3][pos_in_id]<-NA
           }else{ 
             ped[,3][pos_in_id]<-as.character(X.list$id[X.list$X[[i]]$restsire.id[MLsire]])
           }
         }  
ped
}   
