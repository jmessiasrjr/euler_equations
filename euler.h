#ifndef EULER_H_
#define EULER_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/*
void euler(double rhol, double ul, double pl, 
	   double rhor, double ur, double pr, 
           double L, double PD, double gamma,
           double TOUT, double MAXT, double CFL,
           double tol, int n, int BC1, int BC2,
	   int method);
*/
void euler(double dbl[14], int p_int[5], double *et,
           double *t, int *j, int *err);
void gammas(double *g, double gamma);
void vacuum(double *dbl, int *err);
void exact(double *dbl, int *p_int);
void godunov(double *dbl, int *p_int, double *t, int *j);
void CI(double W[][3], double CS[][3], double *dbl,
        int *p_int);
void boundary(double W[][3], int *p_int);
double cfl(double W[][3], double *dbl, double time,
           int *p_int, int m);
void godRP(double W[][3], double F[][3], double *dbl,
           int *p_int);
void roe(double W[][3], double F[][3], double CS[][3], double *dbl,
        double dtdx, int *p_int);
void starvals(double *star, double gamma, double rhok, double uk,
              double ek, double ak, double sig);
void HLL(double W[][3], double F[][3], double CS[][3], double *dbl,
         int *p_int);
void HLLC(double W[][3], double F[][3], double CS[][3], double *dbl,
         int *p_int);
void Rusanov(double W[][3], double F[][3], double CS[][3], double *dbl,
         int *p_int);
void update(double W[][3], double CS[][3], double F[][3],
            double dtdx, double *dbl, int *p_int);
void exactRP(double *var, double *star, double *g,
             double *dbl);
void star_pu(double *star, double *var, double *dbl, double *g);
void rhostar(double *star, double *var, double *g);
void guessp(double *var, double *dbl, double *g, double *pm,
            double *um);
double f(double *var, double *dbl, double *g, double p_old, int i);
double df(double *var, double *dbl, double *g, double p_old, int i);
void flux(double F[][3], double *var, double *dbl, int m);
void sample(double *var, double *star, double *g,
            double *dbl, double S);
void gnuplot(double *dbl, int ex);

#endif /* EULER_H_ */
