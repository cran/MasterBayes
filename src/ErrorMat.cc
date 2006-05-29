#include "ErrorMat.h"

void Error_Mat(Matrix<double> E1, Matrix<double> E2, double **E_mat, int ncat, int *nall, int nloci, bool LogL, bool perlocus){

   int l;
   int cat;
   double e1;
   double e2;
   double ve2;
   int loop;

   for(l= 0; l < nloci; l++){                  // iterates through loci
      for(cat= 0; cat < ncat; cat++){          // iterates through loci
        if(perlocus==FALSE){
          ve2 = E2[cat];
          e1 =  E1[cat]/(1.0+E1[cat]);
          e2 = ve2/(double(nall[l])-1.0);
          loop = cat*7;
        }else{
          ve2 = E2[(l*ncat)+cat];
          e1 =  E1[(l*ncat)+cat]/(1.0+E1[(l*ncat)+cat]);
          e2 = ve2/(double(nall[l])-1.0);
          loop = cat*7;
        }
        if(LogL==true){
         E_mat[l][loop] = log(pow(1-ve2,2));
         E_mat[l][1+loop] = log(e2*(1-ve2)+(e1*pow(1-ve2-e2,2)));
         E_mat[l][2+loop] = log(pow(e2,2));
         E_mat[l][3+loop] = log(pow(1-ve2,2)+pow(e2,2)-(2*e1*pow(1-ve2-e2,2)));
         E_mat[l][4+loop] = log(2*e2*(1-ve2));
         E_mat[l][5+loop] = log(e2*(1-ve2+e2));
         E_mat[l][6+loop] = log(2*pow(e2,2));
       }else{
         E_mat[l][loop] = pow(1-ve2,2);
         E_mat[l][1+loop] = (e2*(1-ve2))+(e1*pow(1-ve2-e2,2));
         E_mat[l][2+loop] = pow(e2,2);
         E_mat[l][3+loop] = pow(1-ve2,2)+pow(e2,2)-(2*e1*pow(1-ve2-e2,2));
         E_mat[l][4+loop] = 2*e2*(1-ve2);
         E_mat[l][5+loop] = e2*(1-ve2+e2);
         E_mat[l][6+loop] = 2*pow(e2,2);
       }
    }
  }
}
 
