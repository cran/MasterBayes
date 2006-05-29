#define _CALCX_G_H
#ifndef _GENERAL_H
#include "General.h"
#endif


void calcX_G(Matrix<double> X_design_G [], int *offid, int noff , int *ndam, int *nsire, int nind, Matrix<int> Dams_vec [], Matrix<int> Sires_vec [], int **G, int nloci, double **A);

void calcX_Gcervus(Matrix<double> X_design_G [], int *offid, int noff , int *ndam, int *nsire, int nind, Matrix<int> Dams_vec [], Matrix<int> Sires_vec [], int **G, int nloci, double **A, double E_cervus);

void calcX_Gmm(Matrix<int> mmD [], Matrix<int> mmS [], int *offid, int noff , int *ndam, int *nsire, int nind, Matrix<int> Dams_vec [], Matrix<int> Sires_vec[],int **G, int nloci);
