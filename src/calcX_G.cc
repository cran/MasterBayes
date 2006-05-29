#include "calcX_G.h"

void calcX_G(Matrix<double> X_design_G [], int *offid, int noff , int *ndam, int *nsire, int nind, Matrix<int> Dams_vec [], Matrix<int> Sires_vec [], int **G, int nloci, double **A){

        int i;
        int l;
        int d;
        int s;
        int records_off;
        int off_poss;
        int sire_poss;
        int dam_poss;
        int o1;
        int o2;
        int d1;
        int d2;
        int s1;
        int s2;
        double o1inp;
        double o2inp;
        double sl;
        double fl;
        double TacP;

    for(i = 0; i < noff ; i++){         // iterate through offspring
      X_design_G[i] = ones<double>(ndam[i]*nsire[i], 1); 
    }

       for(l = 0; l<nloci ; l++){             // iterates through loci  

         for(i = 0; i < noff ; i++){         // iterate through offspring
           off_poss = offid[i];   
           o1 = G[off_poss][l*2];
           o2 = G[off_poss][(l*2)+1];

            for(d = 0; d < ndam[i]; d++){    // gets position of offspring genotype in G

              dam_poss = Dams_vec[i][d];       // gets position of dam genotype in G

              if(dam_poss<nind){                 // if the dam is sampled (and not a zero_prob) 
               d1 = G[dam_poss][l*2];
               d2 = G[dam_poss][(l*2)+1]; 
              }

              for(s = 0; s < nsire[i]; s++){      // iterate through sires

               records_off = (nsire[i]*d)+s;
               sire_poss = Sires_vec[i][s];           // gets position of sire genotype in G   

               if(sire_poss<nind){     

                 s1 = G[sire_poss][l*2];
                 s2 = G[sire_poss][(l*2)+1]; 
               }
         
// both sexes sampled 
                 
               if(sire_poss<nind && dam_poss<nind){        
   
                TacP = 0.0;
       
                if((d1==o1 && s1==o2) || (d1==o2 && s1==o1)){TacP +=0.25;}
                if((d1==o1 && s2==o2) || (d1==o2 && s2==o1)){TacP +=0.25;}
                if((d2==o1 && s1==o2) || (d2==o2 && s1==o1)){TacP +=0.25;}
                if((d2==o1 && s2==o2) || (d2==o2 && s2==o1)){TacP +=0.25;}

                X_design_G[i][records_off] *=  TacP;
                }                              

// female sampled, male unsampled 

             if(sire_poss>=nind && dam_poss<nind){                                  
              sl = 0.0;
              fl = 0.0;
              if(d1==o1 || d1==o2){sl += 0.5;}
              if(d2==o1 || d2==o2){sl += 0.5;}
              if(sl>0.25){
                if(o1!=o2){
                  if(o1==d1 || o1==d2){o1inp = 1.0;}else{o1inp = 0.0;}
                  if(o2==d1 || o2==d2){o2inp = 1.0;}else{o2inp = 0.0;}
                  if(o1inp>0.5){
                   fl += A[l][o2];
                  }
                  if(o2inp>0.5){
                   fl += A[l][o1];
                  }
                  fl /= o1inp+o2inp;
                }else{
                  fl = A[l][o1];
                }
              }  
              X_design_G[i][records_off] *=  sl*fl;
              } 

// male sampled, female unsampled 

             if(sire_poss<nind && dam_poss>=nind){                                   
              sl = 0.0;
              fl = 0.0;
              if(s1==o1 || s1==o2){sl += 0.5;}
              if(s2==o1 || s2==o2){sl += 0.5;}
              if(sl>0.25){
                if(o1!=o2){
                  if(o1==s1 || o1==s2){o1inp = 1.0;}else{o1inp = 0.0;}
                  if(o2==s1 || o2==s2){o2inp = 1.0;}else{o2inp = 0.0;}
                  if(o1inp>0.5){
                   fl += A[l][o2];
                  }
                  if(o2inp>0.5){
                   fl += A[l][o1];
                  }
                  fl /= o1inp+o2inp;
                }else{
                  fl = A[l][o1];
                }
              }
              X_design_G[i][records_off] *=  sl*fl;
              }
              // both parents unsampled 

             if(sire_poss>=nind && dam_poss>=nind){
              X_design_G[i][records_off] *=  A[l][o1]*A[l][o2];
              if(o1!=o2){
              X_design_G[i][records_off] *=  2.0;
              }
            }
           // Pr(P|G)  
           }
        }
     }
  }  
}

