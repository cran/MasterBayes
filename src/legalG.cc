#include "MCMCped.h" 

extern "C"{  

using namespace std;

void legalG(
	int *nindP,		 // number of individuals sampled
        int *damP,             // number of candidate dams per offspring
        int *sireP,            // number of candidate sires per offspring
        int *nlociP,		 // number of loci
        int *nallP,              // number of alleles per locus
        int *maxallP,
	double *st_AP,	         // starting allele frequencies
        int *st_GP,              // starting true genotypes   
        bool *legalP 
){        

// pointers to single variables are redefined

int 	nind = nindP[0], 	
	nloci = nlociP[0],      
        *nall = nallP,
        maxall = maxallP[0],
        i,
        l,
        o1,
        o2,
        d1,
        d2,
        s1,
        s2, 
        index=0;

// declare some variable sized arrays

int     *pG,
        **G;
double  *pA,
        **A;
 
pG = new(nothrow) int [2*nind*nloci];
if(pG==NULL)
{
Rprintf("NO MEMORY for G\n");
exit(1);
}
G = new(nothrow) int* [nind];
if(G==NULL)
{
Rprintf("NO MEMORY for G\n");
exit(1);
}
pA = new(nothrow) double [nloci*maxall];
if(pA==NULL)
{
Rprintf("NO MEMORY for A\n");
exit(1);
}
A = new(nothrow) double* [nloci];
if(A==NULL)
{
Rprintf("NO MEMORY for A\n");
exit(1);
}

for (i=0; i<nind; ++i){
G[i] = &pG[index];
index  += (2*nloci);
}

index = 0;

for (i=0; i<nloci; ++i){
A[i] = &pA[index];
index  += maxall;
}
        
    bool TacP;
    bool o1ind;
    bool o1ins;
    bool o2ind;
    bool o2ins;

    read_stG(st_GP, st_AP, nind, nloci, G, A,nall);

    for(i=0; i<nind; ++i){
      for(l=0; l<nloci; ++l){
        if(damP[i]!=nind || sireP[i]!=nind){

          o1 = G[i][l*2];
          o2 = G[i][(l*2)+1];     

          if(damP[i]!=nind && sireP[i]!=nind){

            d1 = G[damP[i]][l*2];
            d2 = G[damP[i]][(l*2)+1];     

            s1 = G[sireP[i]][l*2];
            s2 = G[sireP[i]][(l*2)+1];     

            TacP=false;

            if((d1==o1 && s1==o2) || (d1==o2 && s1==o1)){TacP=true;}
            if((d1==o1 && s2==o2) || (d1==o2 && s2==o1)){TacP=true;}
            if((d2==o1 && s1==o2) || (d2==o2 && s1==o1)){TacP=true;}
            if((d2==o1 && s2==o2) || (d2==o2 && s2==o1)){TacP=true;}

            if(TacP==false){

              legalP[0]=false;

              if(o1==d1 || o1==d2){
                o1ind = true;
              }else{
                o1ind = false;
              }

              if(o2==d1 || o2==d2){
                o2ind = true;
              }else{
                o2ind = false;
              }

              if(o1==s1 || o1==s2){
                o1ins = true;
              }else{
                o1ins = false;
              }

              if(o2==s1 || o2==s2){
                o2ins = true;
              }else{
                o2ins = false;
              }

              if(o1ind==true || o2ind==true || o1ins==true || o2ins==true){ 
                if(o1ind==true || o1ins==true){
                  if(o1ind==true){
                    if((int)rbinom(1.0, 0.5)==1){
                      G[i][(l*2)+1] = s1;
                    }else{
                      G[i][(l*2)+1] = s2;
                    }
                  }else{
                    if((int)rbinom(1.0, 0.5)==1){
                      G[i][(l*2)+1] = d1;
                    }else{
                      G[i][(l*2)+1] = d2;
                    }
                  }
                }else{
                  if(o2ind==true){
                    if((int)rbinom(1.0, 0.5)==1){
                      G[i][(l*2)] = s1;
                    }else{
                      G[i][(l*2)] = s2;
                    }
                  }else{
                    if((int)rbinom(1.0, 0.5)==1){
                      G[i][(l*2)] = d1;
                    }else{
                      G[i][(l*2)] = d2;
                    }
                  }
                }
              }else{
                if((int)rbinom(1.0, 0.5)==1){
                  G[i][(l*2)] = d1;
                }else{
                  G[i][(l*2)] = d2;
                }
                if((int)rbinom(1.0, 0.5)==1){
                  G[i][(l*2)+1] = s1;
                }else{
                  G[i][(l*2)+1] = s2;
                }
              }
            }
          }else{            
            if(damP[i]!=nind){
              d1 = G[damP[i]][l*2];
              d2 = G[damP[i]][(l*2)+1];     
              if(o1!=d1 &&  o1!=d2 && o2!=d1 && o2!=d2){
                legalP[0]=false;
                if((int)rbinom(1.0, 0.5)==1){
                  G[i][(l*2)+1] = d1;
                }else{
                  G[i][(l*2)+1] = d2;
                }
              }
            }else{
              s1 = G[sireP[i]][l*2];
              s2 = G[sireP[i]][(l*2)+1];     
              if(o1!=s1 &&  o1!=s2 && o2!=s1 && o2!=s2){
                legalP[0]=false;
                if((int)rbinom(1.0, 0.5)==1){
                  G[i][(l*2)+1] = s1;
                }else{
                  G[i][(l*2)+1] = s2;
                }
              }
            } 
           if(G[i][(l*2)]==-999){G[i][(l*2)] = rmultinom_size1(A[l], nall[l]);}           
          }
        }else{   // if parents have missing genotypes sample from allele frequencies
          if(G[i][l*2]==-999){G[i][l*2] = rmultinom_size1(A[l], nall[l]);}
          if(G[i][(l*2)+1]==-999){G[i][(l*2)+1] = rmultinom_size1(A[l], nall[l]);}
        }
      }
    }     
    
   int records = 0;           

   for(l = 0; l < nloci; l++){  
      for(i = 0; i < nind; i++){	                        
          st_GP[records] = G[i][(l*2)];  
          records ++;
          st_GP[records] = G[i][(l*2)+1];  
          records ++;
       }
    }
  }
}





















