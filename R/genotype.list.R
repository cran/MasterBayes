"genotype.list"<-function(G, ...){
gens<-list()
for(i in 1:(length(G[1,])/2)){
gens[[i]]<-genotype(as.matrix(G[,((i*2)-1):(i*2)]))
names(gens)[i]<-names(G[i*2])} 
gens
}
