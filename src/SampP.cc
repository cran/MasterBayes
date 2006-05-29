#include "SampP.h"

void sampP(int *offid, int noff, Matrix<double> X_design_G [], int *npar, int *DSuu, Matrix<double> X_design_betaDus [], Matrix<double> X_design_betaSus [], Matrix<double> X_design_betaDSus [], Matrix<double> X_design_betaDs [], Matrix<double> X_design_betaSs [], Matrix<double> X_design_betaDSs [], int *dam, int *sire, Matrix<double> beta, Matrix<double> us, int *usdamcat, int *ussirecat, int nusd, int nuss, int *ndam, int *nsire, int *ntdam, int *ntsire, Matrix<int> Dams_vec [], Matrix<int> Sires_vec [], int nind, int nmerge, int *mergeV, int *mergeUS, Matrix<double> mergeN []){

	int i;
	int d;
        int s;
        int b;
        int mvar;
        double n1;
        double n2;
        double inv_log_beta;
        int cnt=0;
        int newpar;
        int new_dam_pos;
        int new_sire_pos;
        double S_vec;
        double mean_vec;
        double n;
        double N;
        double tot_mean;
        bool damV;
        bool sireV;
        bool damsireV;
         
         Matrix<double> betaDus (npar[0],1); 

           for(i = 0; i < npar[0]; i++){   
              betaDus [i] = beta[cnt];
              cnt++;
           }
           
         Matrix<double> betaDs (npar[1],1); 

           for(i = 0; i < npar[1]; i++){   
              betaDs [i] = beta[cnt];
              cnt++;
           }

         if((npar[0]+npar[1])>0){damV=true;}else{damV=false;}

         Matrix<double> betaSus (npar[2],1); 
           for(i = 0; i < npar[2]; i++){   
              betaSus [i] = beta[cnt];
              cnt++;
           }

         Matrix<double> betaSs (npar[3],1); 
           for(i = 0; i < npar[3]; i++){   
              betaSs [i] = beta[cnt];
              cnt++;
           }

         if((npar[2]+npar[3])>0){sireV=true;}else{sireV=false;}

         Matrix<double> betaDSus (npar[4],1); 
           for(i = 0; i < npar[4]; i++){   
              betaDSus [i] = beta[cnt];
              cnt++;
           }

         Matrix<double> betaDSs (npar[5],1); 
           for(i = 0; i < npar[5]; i++){   
              betaDSs [i] = beta[cnt];
              cnt++;
           }

         if((npar[4]+npar[5])>0){damsireV=true;}else{damsireV=false;}

	Matrix<double> Dpred;
	Matrix<double> Spred;	
        Matrix<double> DSpred;
	Matrix<double> Dpreds;
	Matrix<double> Spreds;	
        Matrix<double> DSpreds;

	for(i = 0; i < noff; i++){   

/***********************/
/* reparameterise beta */
/***********************/

           if(nmerge>0){
             for(d = 0; d < nmerge; d++){  
               mvar = mergeV[d];
               n1 = mergeN[d](0,i);
               n2 = mergeN[d](1,i);
               inv_log_beta = exp(beta[d])/(1.0+exp(beta[d]));
               if(mvar<(npar[0]+npar[1])){
                 if(mergeUS[d]==0){
                   n1 += us[usdamcat[i]];
                 }
                 if(mergeUS[d]==1){
                   n2 += us[usdamcat[i]];
                 }
                 if(mvar<npar[0]){
                   betaDus[mvar] = inv_log_beta/n1;
                   betaDus[mvar] /= ((inv_log_beta/n1)+((1.0-inv_log_beta)/n2));
                   betaDus[mvar] = log(betaDus[mvar]/(1.0-betaDus[mvar]));
                 }else{
                   betaDs[mvar] = inv_log_beta/n1;
                   betaDs[mvar] /= ((inv_log_beta/n1)+((1.0-inv_log_beta)/n2));
                   betaDs[mvar] = log(betaDs[mvar]/(1.0-betaDs[mvar]));
                 }
               }else{
                 mvar -= npar[0]+npar[1];
                 if(mergeUS[d]==0){
                   n1 += us[nusd+ussirecat[i]];
                 }
                 if(mergeUS[d]==1){
                   n2 += us[nusd+ussirecat[i]];
                 }
                 if(mvar<npar[2]){
                   betaSus[mvar] = inv_log_beta/n1;
                   betaSus[mvar] /= ((inv_log_beta/n1)+((1.0-inv_log_beta)/n2));
                   betaSus[mvar] = log(betaSus[mvar]/(1.0-betaSus[mvar]));
                 }else{
                   betaSs[mvar] = inv_log_beta/n1;
                   betaSs[mvar] /= ((inv_log_beta/n1)+((1.0-inv_log_beta)/n2));
                   betaSs[mvar] = log(betaSs[mvar]/(1.0-betaSs[mvar]));
                 }
               }
             }    
           }

/***********************/
/* sample missing data */
/***********************/

            if(npar[0]!=0){   
              Dpred = exp(X_design_betaDus[i]*betaDus);
              mean_vec = 0.0;
              S_vec = 0.0;
              for(d=0; d<ntdam[i]; d++){
                mean_vec += Dpred[d];
                S_vec +=  pow(Dpred[d],2.0); 
              }
              mean_vec -= Dpred[ndam[i]-1];           
              S_vec -= pow(Dpred[ndam[i]-1], 2.0);
              mean_vec /= (ntdam[i]-1);
              S_vec  /= (ntdam[i]-1);
              S_vec -=  mean_vec*mean_vec; 
              n = double(ntdam[i]-1);
              N = n+us[usdamcat[i]];
              S_vec *= N/(n*(N-n));
              Dpred[ndam[i]-1] = rnorm(mean_vec, sqrt(S_vec)); 
            }
            if(npar[2]!=0){          
              Spred = exp(X_design_betaSus[i]*betaSus);
              mean_vec = 0.0;
              S_vec = 0.0;
              for(s=0; s<ntsire[i]; s++){
                mean_vec += Spred[s];
                S_vec +=  pow(Spred[s],2.0); 
              }
              mean_vec -= Spred[nsire[i]-1];        
              S_vec -= pow(Spred[nsire[i]-1],2.0);
              mean_vec /= (ntsire[i]-1);
              S_vec  /= (ntsire[i]-1);
              S_vec -=  mean_vec*mean_vec; 
              n = double(ntsire[i]-1);
              N = n+us[ussirecat[i]+nusd];
              S_vec *= N/(n*(N-n));
              Spred[nsire[i]-1] = rnorm(mean_vec, sqrt(S_vec));
            }
            if(npar[4]!=0){       
              DSpred = exp(X_design_betaDSus[i]*betaDSus);
              tot_mean=0.0;
              if(DSuu[0]==1){ 
                tot_mean=0.0;
                for(s=0; s<ntsire[i]; s++){
                  mean_vec = 0.0;
                  S_vec = 0.0;   
                  for(d=0; d<ntdam[i]; d++){                   
                    mean_vec += DSpred[(d*ntsire[i])+s];
                    S_vec += pow(DSpred[(d*ntsire[i])+s],2);
                  }
                  mean_vec -= DSpred[((ndam[i]-1)*ntsire[i])+s];
                  S_vec -= pow(DSpred[((ndam[i]-1)*ntsire[i])+s],2);                        
                  mean_vec /= (ntdam[i]-1);
                  S_vec  /= (ntdam[i]-1);
                  S_vec -=  mean_vec*mean_vec; 
                  n = double(ntdam[i]-1);
                  N = n+us[usdamcat[i]];
                  S_vec *= N/(n*(N-n));
                  DSpred[((ndam[i]-1)*ntsire[i])+s] = rnorm(mean_vec, sqrt(S_vec));
                  tot_mean += DSpred[((ndam[i]-1)*ntsire[i])+s];
                }
              }
              if(DSuu[1]==1){
                for(d=0; d<ntdam[i]; d++){
                  mean_vec = 0.0;
                  S_vec = 0.0;   
                  for(s=0; s<ntsire[i]; s++){                   
                    mean_vec += DSpred[(d*ntsire[i])+s];
                    S_vec += pow(DSpred[(d*ntsire[i])+s],2);
                  }
                  mean_vec -= DSpred[(d*ntsire[i])+(nsire[i]-1)];
                  S_vec -= pow(DSpred[(d*ntsire[i])+(nsire[i]-1)],2);          
                  mean_vec /= (ntsire[i]-1);
                  S_vec  /= (ntsire[i]-1);
                  S_vec -=  mean_vec*mean_vec; 
                  n = double(ntsire[i]-1);
                  N = n+us[ussirecat[i]+nusd];
                  S_vec *= N/(n*(N-n));
                  DSpred[(d*ntsire[i])+(nsire[i]-1)] = rnorm(mean_vec, sqrt(S_vec));
                  tot_mean += DSpred[(d*ntsire[i])+(nsire[i]-1)];
                }
              } 
              if(DSuu[0]==1 && DSuu[1]==1){
                DSpred[((ndam[i]-1)*ntsire[i])+(nsire[i]-1)] = tot_mean;
              }
            }

/* combine complete and incomplete design matrices */
 
            if(npar[1]>0){
              Dpreds = Matrix<double>(ndam[i],1);
              for(d=0; d<ndam[i]; d++){
                 for(b=0; b<npar[1]; b++){
                   Dpreds[d] += X_design_betaDs[i](d,b)*betaDs[b]; 
                 }
              }      
              Dpreds = exp(Dpreds);   
              if(npar[0]>0){  
                  for(d=0; d<ndam[i]; d++){
                    Dpreds[d] *= Dpred[d];
                }
              }
            }else{       
              if(npar[0]>0){ 
                Dpreds = Matrix<double>(ndam[i],1);     
                for(d=0; d<ndam[i]; d++){
                  Dpreds[d] = Dpred[d];
                }
              }
            }

            if(npar[3]>0){
              Spreds = Matrix<double>(nsire[i],1);
              for(s=0; s<nsire[i]; s++){
                for(b=0; b<npar[3]; b++){
                  Spreds[s] += X_design_betaSs[i](s,b)*betaSs[b]; 
                }
              }      
              Spreds = exp(Spreds);
              if(npar[2]>0){  
                for(s=0; s<nsire[i]; s++){
                  Spreds[s] *= Spred[s];
                }
              }
            }else{
              if(npar[2]>0){ 
                Spreds = Matrix<double>(nsire[i],1);
                for(s=0; s<nsire[i]; s++){
                  Spreds[s] = Spred[s];
                }
              }
            }

            if(npar[5]>0){
              DSpreds = Matrix<double>(nsire[i]*ndam[i],1);
              for(d=0; d<ndam[i]; d++){
                for(s=0; s<nsire[i]; s++){            
                  for(b=0; b<npar[5]; b++){
                    DSpreds[(d*nsire[i])+s] += X_design_betaDSs[i]((d*ntsire[i])+s, b)*betaDSs[b]; 
                  }
                }
              }   
              DSpreds = exp(DSpreds);
              if(npar[4]>0){  
                for(d=0; d<ndam[i]; d++){
                  for(s=0; s<nsire[i]; s++){            
                    DSpreds[(d*nsire[i])+s] *= DSpred[(d*ntsire[i])+s]; 
                  }
                }
              }      
            }else{
              if(npar[4]>0){  
                DSpreds = Matrix<double>(nsire[i]*ndam[i],1);
                for(d=0; d<ndam[i]; d++){
                  for(s=0; s<nsire[i]; s++){            
                    DSpreds[(d*nsire[i])+s] = DSpred[(d*ntsire[i])+s]; 
                  }
                }
              }
            }

            if(damsireV==false){ 
              if(sireV==false && damV==false){ 
                DSpreds = X_design_G[i];
              }else{
                if(sireV && damV){ 
                  DSpreds = vecc(Spreds*t(Dpreds));
                  for(d=0; d<(ndam[i]*nsire[i]); d++){
                    DSpreds[d] *= X_design_G[i][d];
                  }
                }else{
                  DSpreds = ones<double>(ndam[i]*nsire[i],1);
                  if(sireV){ 
                    cnt = 0;
                    for(d=0; d<ndam[i]; d++){
                      for(s=0; s<nsire[i]; s++){
                        DSpreds[cnt] = Spreds[s]*X_design_G[i][cnt];
                        cnt++;
                      }
                    }
                  }else{
                    cnt = 0;
                    for(d=0; d<ndam[i]; d++){
                      for(s=0; s<nsire[i]; s++){
                        DSpreds[cnt] = Dpreds[d]*X_design_G[i][cnt];
                        cnt++;
                      }
                    }
                  }
                }  
              }
            }else{
              if(sireV==false && damV ==false){
                for(d=0; d<(ndam[i]*nsire[i]); d++){
                  DSpreds[d] *= DSpreds[d]*X_design_G[i][d];
                }
              }else{
                if(sireV && damV){ 
                  cnt = 0;
                  for(d=0; d<ndam[i]; d++){
                    for(s=0; s<nsire[i]; s++){
                      DSpreds[cnt] *= Dpreds[d]*Spreds[s]*X_design_G[i][cnt];
                      cnt++;
                    }
                  }
                }else{
                  if(sireV){ 
                    cnt = 0;
                    for(d=0; d<ndam[i]; d++){
                      for(s=0; s<nsire[i]; s++){
                        DSpreds[cnt] *= Spreds[s]*X_design_G[i][cnt];
                        cnt++;
                      }
                    }
                  }else{
                    cnt = 0;
                    for(d=0; d<ndam[i]; d++){
                      for(s=0; s<nsire[i]; s++){
                        DSpreds[cnt] *= Dpreds[d]*X_design_G[i][cnt];
                        cnt++;
                      }
                    }
                  }
                }
              }
            }

            if(nuss>0){
              for(d=0; d<ndam[i]; d++){
                DSpreds[(nsire[i]*(d+1))-1] *= us[ussirecat[i]+nusd];
              }
            }


            if(nusd>0){
              for(s=0; s<nsire[i]; s++){
                DSpreds[(nsire[i]*(ndam[i]-1))+s] *= us[usdamcat[i]];
              }
            }
            newpar = rmultinom_size1M(DSpreds, nsire[i]*ndam[i]);      // samples a parental combination
 
            new_dam_pos = newpar/nsire[i];
            new_sire_pos = newpar-(new_dam_pos*nsire[i]);
		
            dam[offid[i]] = Dams_vec[i][new_dam_pos];     // puts new dam.id into pedigree
            sire[offid[i]] = Sires_vec[i][new_sire_pos];  // puts new sire.id into pedigre        
 
        }

}



