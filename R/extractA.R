"extractA"<-function(G,marker.type="MS", ...){

   if(is.genotype(G[[1]])==FALSE & is.genotypeD(G[[1]])==FALSE){
     if("id"%in%colnames(G)){
       G<-G[,-which(colnames(G)=="id")]
     }
     if("categories"%in%colnames(G)){
       G<-G[,-which(colnames(G)=="categories")]
     }
     G<-genotype.list(G,marker.type=marker.type)
   }

  A<-lapply(G, function(x){summary(x)$allele.freq[,"Proportion"][1:nallele(x)]})

  A
}
