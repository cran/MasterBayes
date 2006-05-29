"legalG"<-function(G, A, ped, time_born=NULL, ...){

  nind<-length(ped[,1])
  nloci<-length(A)
  nall<-unlist(lapply(A, length))
  maxall<-max(nall)
  namesG<-names(G)

  for(i in 1:length(A)){  
    G[[i]]<-matrix(match(allele(G[[i]]),names(A[[i]]), nomatch=-998), length(G[[i]]),2)
  } 

  #################################### order pedigree and data #############################################

  oped<-order.ped(ped, time_born=NULL)

  rearrange_data<-match(oped[,1], ped[,1])

  dam<-match(oped[,2], oped[,1])-1
  dam[which(is.na(dam)==T)]<-nind
  sire<-match(oped[,3], oped[,1])-1
  sire[which(is.na(sire)==T)]<-nind

  G<-lapply(G, function(x){x[rearrange_data,]})

  ###########################################################################################################

  G<-c(t(matrix(unlist(G), nind,2*nloci)))-1

  legal<-TRUE

  output<-.C("legalG",
	as.integer(nind),		 
	as.integer(dam),		
	as.integer(sire),		
	as.integer(nloci),		
	as.integer(nall),		
	as.integer(maxall),		
        as.double(unlist(A)),                   
        as.integer(G),                  
        as.logical(legal))

  tmp<-array(output[[8]], c(2, length(ped[,1]), length(A)))+1

  G<-as.data.frame(matrix(NA, length(ped[,1]), 2*length(A)))

  for(i in 1:length(A)){
    G[,c(((i*2)-1):(i*2))]<-t(tmp[,,i])
    G[,(i*2)-1]<-names(A[[i]])[G[,(i*2)-1]]
    G[,(i*2)]<-names(A[[i]])[G[,(i*2)]]
  }

  G<-G[match(ped[,1], oped[,1]),]
  G<-genotype.list(G)
  names(G)<-namesG
  list(G=G,valid=output[[9]]) 

}
       
