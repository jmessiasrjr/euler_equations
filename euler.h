#ifndef EULER_H_
#define EULER_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>

/*
void euler(double rhol, double ul, double pl,
       double rhor, double ur, double pr,
           double L, double PD, double gamma,
           double TOUT, double MAXT, double CFL,
           double tol, int n, int BC1, int BC2,
       int method);
*/
void euler(const std::vector<double> &values, const std::vector<int> &parameters, double &et,
        double &t, int &j, int &err);
void gammas(std::vector<double> &g, const double gamma);
void vacuum(const std::vector<double> &values, int &err);
void exact(const std::vector<double> &values, const std::vector<int> &parameters);
void godunov(const std::vector<double> &values, const std::vector<int> &parameters, double &t, int &j);
void CI(std::vector<std::vector<double> > &W, std::vector<std::vector<double> > &CS,
        const std::vector<double> &values, const std::vector<int> &parameters);
void boundary(std::vector<std::vector<double>> &W, const std::vector<int> &parameters);
double cfl(const std::vector<std::vector<double>> &W, const std::vector<double> &values,
           const std::vector<int> &parameters, const double time, const int m);
void godRP(std::vector<std::vector<double>> &W, std::vector<std::vector<double>> &F,
           const std::vector<double> &values, const std::vector<int> &parameters);
void roe(std::vector<std::vector<double>> &W, std::vector<std::vector<double>> &F,
         std::vector<std::vector<double>> &CS, const std::vector<double> &values,
         const std::vector<int> &parameters, double dtdx);
void starvals(std::vector<double> &star, const double gamma, const double rhok,
              const double uk, const double ek, const double ak, const double sig);
void HLL(std::vector<std::vector<double>> &W, std::vector<std::vector<double>> &F,
         std::vector<std::vector<double>> &CS, const std::vector<double> &values,
         const std::vector<int> &parameters);
void HLLC(std::vector<std::vector<double>> &W, std::vector<std::vector<double>> &F,
          std::vector<std::vector<double>> &CS, const std::vector<double> &values,
          const std::vector<int> &parameters);
void Rusanov(std::vector<std::vector<double>> &W, std::vector<std::vector<double>> &F,
             std::vector<std::vector<double>> &CS, const std::vector<double> &values,
             const std::vector<int> &parameters);
void update(std::vector<std::vector<double>> &W, std::vector<std::vector<double>> &F,
            std::vector<std::vector<double>> &CS, const std::vector<double> &values,
            const std::vector<int> &parameters, double dtdx);
void exactRP(const std::vector<double> &var, std::vector<double> &star,
             const std::vector<double> &g, const std::vector<double> &values);
void star_pu(const std::vector<double> &var, std::vector<double> &star,
              const std::vector<double> &g, const std::vector<double> &values);
void rhostar(const std::vector<double> &var, std::vector<double> &star,
             const std::vector<double> &g);
void guessp(const std::vector<double> &var, const std::vector<double> &values,
            const std::vector<double> &g, double &pm, double &um);
double f(const std::vector<double> &var, const std::vector<double> &values,
         const std::vector<double> &g, const double p_old, const int i);
double df(const std::vector<double> &var, const std::vector<double> &values,
          const std::vector<double> &g, double p_old, int i);
void flux(std::vector<std::vector<double>> &F, const std::vector<double> &var,
          const std::vector<double> &values, const int m);
void sample(std::vector<double> &var, std::vector<double> &star,
            const std::vector<double> &g, const std::vector<double> &values, const double S);
void gnuplot(const std::vector<double> &values, const int ex);

#endif /* EULER_H_ */
