/* 
* array of Pr(Gobs|G, E1, E2 , n(allele)) for the seven possible combinations:
*
* No Mismatch   Gobs    G
*     0         Hom     Hom
*     1         Hom     Het
*     2         Hom      -
*     0         Het     Het
*     1         Het     Hom
*     1         Het     Het
*     2         Het      -
*
* Defined for each locus and each category of samples over which error rates vary.
* See Wang J.L. 2004 Genetics 166 4 1963-1979    
*/


#define _ERRORMAT_H
#ifndef _GENERAL_H
#include "General.h"
#endif
 
void Error_Mat(Matrix<double> E1, Matrix<double> E2, double **E_mat, int ncat, int *nall, int nloci, bool LogL, bool perlocus, int mtype);

/* 
* E1 [ncat] = vector of error rates for E1 
* E2 [ncat] = vector of error rates for E2 
* E_mat [nloci][ncat*7]= array of likelihoods
* ncat = number of categories over which error rates vary
* nall [nloci] = vector of number of alleles at ecah locus
* nloci = number of loci
* LogL = logical value indicating whether the log-likelihood or the likelihood is stored
*/
