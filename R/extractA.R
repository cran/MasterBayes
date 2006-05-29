"extractA"<-function(G,...){

   if(is.genotype(G[[1]])==FALSE){
     if("id"%in%colnames(G)){
       G<-G[,-which(colnames(G)=="id")]
     }
     if("categories"%in%colnames(G)){
       G<-G[,-which(colnames(G)=="categories")]
     }
     G<-genotype.list(G)
   }

  A<-lapply(G, function(x){summary(x)$allele.freq[,2][1:nallele(x)]})

  A
}
