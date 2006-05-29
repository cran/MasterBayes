#define _READ_IN_DATA_H
#ifndef _GENERAL_H
#include "General.h"
#endif

void read_Gobs(int *GobsP, int nsamp, int nloci, int **Gobs);

void read_stG(int *st_GP, double *st_AP, int nind, int nloci, int **G, double **A, int *nall);

void read_X_beta(int noff, int *ndam, int *nsire, int *npar, double *X_design_betaDusP, double *X_design_betaSusP, double *X_design_betaDSusP,  double *X_design_betaDsP, double *X_design_betaSsP, double *X_design_betaDSsP, Matrix<double> X_design_betaDus [], Matrix<double> X_design_betaSus [], Matrix<double> X_design_betaDSus [], Matrix<double> X_design_betaDs [], Matrix<double> X_design_betaSs [], Matrix<double> X_design_betaDSs []);

void read_stP(int noff, int *ndam, int *damid, int *nsire, int *sireid, map<int,int> Dams [], map<int,int> Sires [], Matrix<int> Dams_vec [], Matrix<int> Sires_vec []);
