fillX.G<-function(X.list, A, G, E2=0.005){

       noff<-length(X.list$X)
       ndam<-c(unlist(lapply(X.list$X,function(x){length(x$dam.id)})))	
       nsire<-c(unlist(lapply(X.list$X,function(x){length(x$sire.id)})))	
       offid<-as.numeric(names(X.list$X))-1	                                
       damid<-c(unlist(lapply(X.list$X,function(x){x$dam.id})))-1
       sireid<-c(unlist(lapply(X.list$X,function(x){x$sire.id})))-1

       nall<-sapply(A, length)
       maxall<-max(nall)
       nloci<-length(nall)
       nind<-length(X.list$id)
       nbeta<-rep(0,5)

       if(is.null(E2)){
         E2<-0.005
       }
       E2<-E2*(2-E2)

       for(i in 1:length(A)){  
         G[[i]]<-matrix(match(allele(G[[i]]),names(A[[i]]), nomatch=-998), length(G[[i]]),2)
       } 

       G<-c(t(matrix(unlist(G), nind,2*nloci)))-1

       X_design_G<-rep(0, sum(ndam*nsire))

output<-.C("fillXG",
        as.integer(nind),       # number of individuals sampled
        as.integer(noff),       # number of non-base (offspring) individuals
        as.integer(ndam),       # number of candidate dams per offspring
        as.integer(nsire),      # number of candidate sires per offspring
        as.integer(nloci),	# number of loci
        as.integer(nall),       # number of alleles per locus
        as.integer(maxall),     # number of alleles at most polymorhic locus
        as.integer(nbeta),
        as.integer(offid),      # offspring id
        as.integer(damid),      # candidate dam id's for each offspring
        as.integer(sireid),	# candidate sire id's for each offspring	
        as.double(X_design_G),  # Mendelian transition probabilities dam and sire sampled			
        as.double(unlist(A)),	# starting allele frequencies
        as.double(E2),	        # starting values of E1 and E2
        as.integer(G)           # starting true genotypes    
)


startG<-1
startD<-1
startS<-1

for(i in 1:noff){
X.list$X[[i]]$G<-as.matrix(output[[12]][(startG-1)+1:(ndam[i]*nsire[i])])
startG<-startG+ndam[i]*nsire[i]
}
X.list
}
