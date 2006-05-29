"modeP"<-function(postP, threshold=0, ...){

  ped<-matrix(NA, length(postP), 3)
  ped[,1]<-names(postP)
  ped[,2]<-unlist(lapply(postP, function(x){colnames(x)[which.max(colSums(x))]})) 
  post_prob<-unlist(lapply(postP, function(x){colSums(x)[which.max(colSums(x))]/sum(x)})) 
  ped[,2][which(post_prob<threshold | ped[,2]=="USdam")]<-NA
  ped[,3]<-unlist(lapply(postP, function(x){rownames(x)[which.max(rowSums(x))]})) 
  post_prob<-unlist(lapply(postP, function(x){rowSums(x)[which.max(rowSums(x))]/sum(x)})) 
  ped[,3][which(post_prob<threshold | ped[,3]=="USsire")]<-NA

ped
}
 
