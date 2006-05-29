#include "calcX_G.h"

void calcX_Gcervus(Matrix<double> X_design_G [], int *offid, int noff , int *ndam, int *nsire, int nind, Matrix<int> Dams_vec [], Matrix<int> Sires_vec [], int **G, int nloci, double **A, double E_cervus){

        int i;
        int l;
        int d;
        int s;
        int records_off;
        int off_poss;
        int sire_poss;
        int dam_poss;
        int o1 = -999;
        int o2 = -999;
        int d1 = -999;
        int d2 = -999;
        int s1 = -999;
        int s2 = -999;
        double o1inp;
        double o2inp;
        double sl;
        double fl;
        double TacP;
        double pG;
        double sslfl;
        double dslfl;

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
              if(dam_poss<nind){
                d1 = G[dam_poss][l*2];
                d2 = G[dam_poss][(l*2)+1]; 
              }else{
                d1 = -999;
              }

              for(s = 0; s < nsire[i]; s++){      // iterate through sires

               records_off = (nsire[i]*d)+s;
               sire_poss = Sires_vec[i][s];           // gets position of sire genotype in G   

               if(sire_poss<nind){ 
                s1 = G[sire_poss][l*2];
                s2 = G[sire_poss][(l*2)+1]; 
               }else{
                s1 = -999;
               }
// all sampled and all genotypes known

               if(o1!=-999){
                  pG = A[l][o1]*A[l][o2];
                  if(o1!=o2){
                    pG *= 2.0;
                  }
                  if(d1!=-999){ // if dam has genotype
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
                    dslfl = sl*fl;
                  }
                  if(s1!=-999){ // if sire has genotype
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
                    sslfl = sl*fl;
                  }
                  if(s1==-999 && d1==-999){       // if nobody sampled
                    X_design_G[i][records_off] *= pG; 
                  }
                  if(d1!=-999 && s1==-999){      // if dam sampled only
                    X_design_G[i][records_off] *=  (pow(1-E_cervus,2.0)*dslfl)+pG*E_cervus*(2.0-E_cervus);
                  }
                  if(s1!=-999 &&  d1==-999){      // if sire sampled only
                    X_design_G[i][records_off] *=  (pow(1-E_cervus,2.0)*sslfl)+pG*E_cervus*(2.0-E_cervus);
                  }
                  if(s1!=-999 && d1!=-999){      // if both sampled only
                    TacP = 0.0;
       
                    if((d1==o1 && s1==o2) || (d1==o2 && s1==o1)){TacP +=0.25;}
                    if((d1==o1 && s2==o2) || (d1==o2 && s2==o1)){TacP +=0.25;}
                    if((d2==o1 && s1==o2) || (d2==o2 && s1==o1)){TacP +=0.25;}
                    if((d2==o1 && s2==o2) || (d2==o2 && s2==o1)){TacP +=0.25;}

                    TacP *= pow(1-E_cervus,3.0);
                    X_design_G[i][records_off] *=  ((E_cervus*pow(1-E_cervus,2.0)*(sslfl+dslfl+pG))+TacP+(pG*pow(E_cervus,2.0)*(4.0-(3.0*E_cervus))));
                  }
              }
           }
        }
     }
  }  

}

