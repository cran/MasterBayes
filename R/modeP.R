"modeP"<-function(postP, threshold=0, ...){
  if(is.list(postP)){
    ped<-matrix(NA, length(postP), 3)
    ped[,1]<-names(postP)
    ped[,2]<-unlist(lapply(postP, function(x){colnames(x)[which.max(colSums(x))]})) 
    post_prob<-unlist(lapply(postP, function(x){colSums(x)[which.max(colSums(x))]/sum(x)})) 
    ped[,2][which(post_prob<threshold)]<-NA
    ped[,2][which(ped[,2]=="USdam")]<-"USdam"
    ped[,3]<-unlist(lapply(postP, function(x){rownames(x)[which.max(rowSums(x))]})) 
    post_prob<-unlist(lapply(postP, function(x){rowSums(x)[which.max(rowSums(x))]/sum(x)})) 
    ped[,3][which(post_prob<threshold)]<-NA
    ped[,3][which(ped[,3]=="USsire")]<-"USsire"
  }else{
    ped<-matrix(NA, dim(postP)[1], 3)
    lpost<-dim(postP)[2]
    ped[,1]<-rownames(postP)
    postP<-apply(postP, 1, function(x){table(paste(x[seq(1,lpost,2)], x[seq(2,lpost,2)]))})
    ped[,2]<-unlist(lapply(postP, function(x){strsplit(names(x)[which.max(x)], " ")[[1]][1]})) 
    ped[,3]<-unlist(lapply(postP, function(x){strsplit(names(x)[which.max(x)], " ")[[1]][2]}))
    post_prob<-unlist(lapply(postP, function(x){x[which.max(x)]/sum(x)})) 
    ped[,2:3][which(post_prob<threshold),]<-NA
    ped[,2][which(ped[,2]=="USdam")]<-"USdam"
    ped[,3][which(ped[,3]=="USsire")]<-"USsire"
  }
ped
}
 
