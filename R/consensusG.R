consensusG<-function(GdP, cat.levels=NULL){

  Gid<-GdP$id[-duplicated(GdP$id)==FALSE]
  G<-lapply(GdP$G, function(x){x[-duplicated(GdP$id)==FALSE]})

  if(is.null(cat.levels)){
    GdP$categories<-NULL
  }else{
    if(any(GdP$categories%in%cat.levels==FALSE)){stop("some GdP$categories not appearing in cat.levels")}
  }
  for(i in 1:length(G)){
    nG<-GdP$G[[i]][order(is.na(GdP$G[[i]])*10000+match(GdP$categories,cat.levels))]
    nid<-GdP$id[order(is.na(GdP$G[[i]])*10000+match(GdP$categories,cat.levels))]
    nG<-nG[-duplicated(nid)==FALSE]
    nid<-nid[-duplicated(nid)==FALSE]
    G[[i]][match(nid, Gid)]<-nG
  }
  GdataPed(id=Gid, G=G, perlocus=GdP$perlocus, marker.type=GdP$marker.type)
}


