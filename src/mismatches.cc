#include "MCMCped.h" 

extern "C"{  

using namespace std;

void mismatches(
	int *nindP,		 // number of individuals sampled
        int *noffP,              // number of non-base (offspring) individuals
        int *ndamP,              // number of candidate dams per offspring
        int *nsireP,             // number of candidate sires per offspring
        int *nlociP,		 // number of loci
        int *offidP,             // offspring id
        int *damidP,             // candidate dam id's for each offspring
        int *sireidP,		 // candidate sire id's for each offspring	
        int *mmDP,             // number of misamtches per dam
        int *mmSP,             // number of misamtches per sire     
        int *st_GP              // starting true genotypes    
){         
// pointers to single variables are redefined

int 	nind = nindP[0],
        noff = noffP[0],  	
	nloci = nlociP[0],      
        i,
        l=0;

// declare some variable sized arrays

int     *pG,
        **G;
 
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

for (i=0; i<nind; ++i){
G[i] = &pG[l];
l  += (2*nloci);
}

     
        Matrix<int> mmD [noff];
        Matrix<int> mmS [noff];

	map<int, int> Dams [noff];     // two way indexing vectors
	map<int, int> Sires [noff];    // map[i][dam_id] = n           dam_id is the n^th mother of the i^th individual
        Matrix<int> Dams_vec [noff];   // Matrix[i][n] = dam_id        the n^th mother of the i^th individul is dam.id 
	Matrix<int> Sires_vec [noff];
      
	int records = 0;	// itterates through genotypes
		
        for(i = 0; i < nind; i++){	
          for(l = 0; l < nloci; l++){                                
            G[i][(l*2)] = st_GP[records];  
            records ++;
            G[i][(l*2)+1] = st_GP[records];  
            records ++;                                           
	  }
        }

        read_stP(noff, ndamP, damidP, nsireP, sireidP,Dams,Sires,Dams_vec,Sires_vec);
         
        calcX_Gmm(mmD, mmS, offidP, noff , ndamP, nsireP, nind, Dams_vec, Sires_vec, G, nloci);
 
         int cnt=0;
  
	 for(i = 0; i <noff ; i++){
           for(int d = 0; d <ndamP[i]; d++){
            mmDP[cnt] = mmD[i][d];
            cnt++;
           }
         }
         
         cnt = 0;

 	 for(i = 0; i <noff ; i++){
           for(int s = 0; s <nsireP[i]; s++){
            mmSP[cnt] = mmS[i][s];
            cnt++;
           }
         }
}
}
