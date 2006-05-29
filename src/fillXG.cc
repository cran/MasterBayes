#include "MCMCped.h" 

extern "C"{  

using namespace std;

void fillXG(
	int *nindP,		 // number of individuals sampled
        int *noffP,              // number of non-base (offspring) individuals
        int *ndamP,              // number of candidate dams per offspring
        int *nsireP,             // number of candidate sires per offspring
        int *nlociP,		 // number of loci
        int *nallP,              // number of alleles per locus
        int *maxallP,            // number of alleles at most polymorhic locus
        int *nparP,             // 
        int *offidP,             // offspring id
        int *damidP,             // candidate dam id's for each offspring
        int *sireidP,		 // candidate sire id's for each offspring	
        double *X_design_GP,     // Mendelian transition probabilities dam and sire sampled
	double *AP,	         // starting allele frequencies
	double *EP,	         // starting values of E1 and E2
        int *GP              // starting true genotypes    
){         // if TRUE joint posterior distribution of P is written (default = marginal)

// pointers to single variables are redefined

int 	nind = nindP[0],
        noff = noffP[0],  	
	nloci = nlociP[0],      
        maxall = maxallP[0],
        *nall = nallP,
        *npar = nparP,
        i,
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
        
        Matrix<double> X_design_G [noff];

	map<int, int> Dams [noff];     // two way indexing vectors
	map<int, int> Sires [noff];    // map[i][dam_id] = n           dam_id is the n^th mother of the i^th individual
        Matrix<int> Dams_vec [noff];   // Matrix[i][n] = dam_id        the n^th mother of the i^th individul is dam.id 
	Matrix<int> Sires_vec [noff];
      
      
         read_stG(GP, AP, nind, nloci, G, A,nall);
         read_stP(noff, ndamP, damidP, nsireP, sireidP,Dams,Sires,Dams_vec,Sires_vec); 
         calcX_Gcervus(X_design_G, offidP, noff , ndamP, nsireP, nind, Dams_vec, Sires_vec, G, nloci, A, EP[0]);

         int cnt_ds=0;
         int p;

	 for(i = 0; i <noff ; i++){
           for(p = 0; p <(ndamP[i]*nsireP[i]); p++){
            X_design_GP[cnt_ds] = X_design_G[i][p];
            cnt_ds++;
           }
         }
  }
}
