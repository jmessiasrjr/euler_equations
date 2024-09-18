/*
    Solve one-dimensional Euler equations for an
    ideal gas by Exact solution and Godunov method
    with exact RP, Roe, HLL, HLLC or Rusanov
    Riemann solvers.
*/

#include "euler.h"

void euler(double rhol, double ul, double pl, 
	   double rhor, double ur, double pr, 
           double L, double PD, double gamma,
           double TOUT, double MAXT, double CFL,
           double tol, int n, int BC1, int BC2,
	   int method)
{
    int p_int[4];
    double dbl[13];
    struct timespec ti, tf;

    p_int[0] = n;
    p_int[1] = BC1;
    p_int[2] = BC2;
    p_int[3] = method;

    dbl[0] = rhol;
    dbl[1] = ul;
    dbl[2] = pl;
    dbl[3] = rhor;
    dbl[4] = ur;
    dbl[5] = pr;
    dbl[6] = L;
    dbl[7] = PD;
    dbl[8] = gamma;
    dbl[9] = TOUT;
    dbl[10] = MAXT;
    dbl[11] = CFL;
    dbl[12] = tol;

    clock_gettime(CLOCK_MONOTONIC_RAW, &ti);

    exact(dbl, p_int);

    godunov(dbl, p_int);

    clock_gettime(CLOCK_MONOTONIC_RAW, &tf);

    time_msec(ti, tf);
}

void gammas(double *g, double gamma)
{
    g[0] = (gamma - 1.)/(2.*gamma);
    g[1] = (gamma + 1.)/(2.*gamma);
    g[2] = (2.*gamma)/(gamma - 1.);
    g[3] = 2./(gamma - 1.);
    g[4] = 2./(gamma + 1.);
    g[5] = (gamma - 1.)/(gamma + 1.);
    g[6] = (gamma - 1.)/2.;
    g[7] = gamma - 1.;
};

void exact(double *dbl, int *p_int)
{
	int i, n;
	double rho, u, p, dx, L, PD, gamma, TOUT, S, e,
	       x[10*p_int[0]], g[8], var[6], star[4];
	FILE *exato;

    exato = fopen("exact.out","w");

	n = p_int[0];

    L = dbl[6];
    PD = dbl[7];
    gamma = dbl[8];
    TOUT = dbl[9];

    gammas(g, gamma);

    var[0] = dbl[0];
    var[1] = dbl[1];
    var[2] = dbl[2];
    var[3] = dbl[3];
    var[4] = dbl[4];
    var[5] = dbl[5];

    exactRP(var, star, g, dbl, p_int);

    dx = L/(10*n);

    fprintf(exato,"       x [m]           rho             u  "
            "           p             e\n");

    for(i = 0;i < (10*n);i++)
    {
        x[i] = (i + 0.5)*dx;
        S = ( x[i] - PD )/TOUT;
        fprintf(exato,"%12.4lf  ",x[i]);

        sample(var, star, g, dbl, S, i);

        rho = var[0];
        u = var[1];
        p = var[2];
        e = p/( (gamma - 1)*rho );

        fprintf(exato,"%12.6lf  %12.6lf  %12.6lf  %12.6lf\n",
                       rho, u, p, e);
    }

    fclose(exato);

};

void godunov(double *dbl, int *p_int)
{
    int i, n, method;
    double rho, u, p, x, dx, L, gamma, TOUT, MAXT, e, dt, time = 0.,
           CS[p_int[0]+2][3], F[p_int[0]+2][3], W[p_int[0]+2][3];
    FILE *saida;

    saida = fopen("result.out","w");

    n = p_int[0] + 2;
    method = p_int[3];

    L = dbl[6];
    gamma = dbl[8];
    TOUT = dbl[9];
    MAXT = dbl[10];

    dx = L/(n-2);

    CI(W, CS, dbl, p_int);

    if(method == 1) printf("Roe Solver\n");
    else if(method == 2) printf("HLL Solver\n");
    else if(method == 3) printf("HLLC Solver\n");
    else if(method == 4) printf("Rusanov Solver\n");
    else printf("Godunov Exact RP\n");

    for(i=0;i<MAXT;i++)
    {
        boundary(W, dbl, p_int);
        dt = cfl(W, dbl, time, p_int, i);

        if(method == 1) roe(W, F, CS, dbl, dt/dx, p_int);
        else if(method == 2) HLL(W, F, CS, dbl, p_int);
        else if(method == 3) HLLC(W, F, CS, dbl, p_int);
        else if(method == 4) Rusanov(W, F, CS, dbl, p_int);
        else godRP(W, F, dbl, p_int);

        update(W, CS, F, dt/dx, dbl, p_int);

        time += dt;
        if(fabs(time -TOUT) < 1e-6) break;
    }

    printf("Iterations number: %8d\n",i);
    printf("Time: %19.10g s\n",time);
    fprintf(saida,"       x [m]           rho             u  "
            "           p             e\n");

    for(i=1;i<(n-1);i++)
    {
        rho = W[i][0];
        u = W[i][1];
        p = W[i][2];

        x = (i - 0.5)*dx;
        e = p/( (gamma - 1)*rho );

        fprintf(saida,"%12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf\n",
                x, rho, u, p, e);
    }

    fclose(saida);

};

void CI(double W[][3], double CS[][3], double *dbl,
        int *p_int)
{
    int i, n;
    double x, dx, rhol, ul, pl, rhor, ur, pr, L, PD, gamma;

    n = p_int[0] + 2;

    rhol = dbl[0];
    ul = dbl[1];
    pl = dbl[2];
    rhor = dbl[3];
    ur = dbl[4];
    pr = dbl[5];
    L = dbl[6];
    PD = dbl[7];
    gamma = dbl[8];

    dx = L/(n-2);

    for(i=1;i<(n-1);i++)
    {
        x = (i-.5)*dx;
        if(x < PD)
        {
            W[i][0] = rhol;
            W[i][1] = ul;
            W[i][2] = pl;
        }
        else
        {
            W[i][0] = rhor;
            W[i][1] = ur;
            W[i][2] = pr;
        }

        CS[i][0] = W[i][0];
        CS[i][1] = W[i][0]*W[i][1];
        CS[i][2] = 0.5*CS[i][1]*W[i][1] + W[i][2]/
                          (gamma - 1.);
    }

};

void boundary(double W[][3], double *dbl, int *p_int)
{
    int n, BCL, BCR;

    n = p_int[0] + 2;
    BCL = p_int[1];
    BCR = p_int[2];

    if(BCL == 0) //transmissive
    {
        W[0][0] = W[1][0];
        W[0][1] = W[1][1];
        W[0][2] = W[1][2];
    }
    else //reflective
    {
        W[0][0] = W[1][0];
        W[0][1] = -W[1][1];
        W[0][2] = W[1][2];
    }

    if(BCR == 0) //transmissive
    {
        W[n-1][0] = W[n-2][0];
        W[n-1][1] = W[n-2][1];
        W[n-1][2] = W[n-2][2];
    }
    else //reflective
    {
        W[n-1][0] = W[n-2][0];
        W[n-1][1] = -W[n-2][1];
        W[n-1][2] = W[n-2][2];
    }

};

double cfl(double W[][3], double *dbl, double time,
           int *p_int, int m)
{
    int i, n;
    double a, rho, u, p, L, gamma, TOUT, CFL, dx, dt, S, Smax = 0.;

    n = p_int[0] + 2;

    L = dbl[6];
    gamma = dbl[8];
    TOUT = dbl[9];
    CFL = dbl[11];

    for(i=0;i<n;i++)
    {
        rho = W[i][0];
        u = W[i][1];
        p = W[i][2];
        a = sqrt(gamma*p/rho);
        S = fabs(u+a);
        if(S > Smax) Smax = S;
    }

    dx = L/(n-2);
    dt = CFL*dx/Smax;

    if(m<5) dt = 0.2*dt;
    if( (time+dt) > TOUT ) dt = TOUT - time;

    return(dt);

};

void godRP(double W[][3], double F[][3], double *dbl,
           int *p_int)
{
    int i, n;
    double gamma, var[6], g[8], star[4];

    n = p_int[0] + 2;

    gamma = dbl[8];

    gammas(g, gamma);

    for(i=0;i<(n-1);i++)
    {
        var[0] = W[i][0];
        var[1] = W[i][1];
        var[2] = W[i][2];
        var[3] = W[i+1][0];
        var[4] = W[i+1][1];
        var[5] = W[i+1][2];

        exactRP(var, star, g, dbl, p_int);

        sample(var, star, g, dbl, 0., i);

        flux(F, var, dbl, i);

    }

};

void roe(double W[][3], double F[][3], double CS[][3], double *dbl,
        double dtdx, int *p_int)
{
    int i, j, n;
    double gamma, rhol, ul, pl, el, hl, rhor, ur, pr, er, hr,
           al, ar, rm, rhom, um, hm, am, u_star, a_star, ak,
           sig, cflm, Sml, Smr, Snew, Tolson = 0.1, K[3],
           star[5], FD[p_int[0]+2][3];

    n = p_int[0] + 2;

    gamma = dbl[8];
    um = star[0];       // \u^~
    hm = star[1];       // \H^~
    am = star[2];       // \a^~
    u_star = star[3];
    a_star = star[4];

    for(i=0;i<n;i++)
    {
        if( i==0 || i==(n-1) )
        {
            CS[i][0] = W[i][0];
            CS[i][1] = W[i][0]*W[i][1];
            CS[i][2] = 0.5*( W[i][0]*W[i][1]*W[i][1]) +
                       W[i][2]/(gamma - 1.);
        }
        FD[i][0] = CS[i][1];
        FD[i][1] = CS[i][1]*W[i][1] + W[i][2];
        FD[i][2] = W[i][1]*( CS[i][2] + W[i][2] );
    }

    for(i=0;i<(n-1);i++)
    {
        rhol = W[i][0];
        ul = W[i][1];
        pl = W[i][2];
        el = CS[i][2];
        hl = (el + pl)/rhol;

        rhor = W[i+1][0];
        ur = W[i+1][1];
        pr = W[i+1][2];
        er = CS[i+1][2];
        hr = (er + pr)/rhor;

        al = sqrt(gamma*pl/rhol);
        ar = sqrt(gamma*pr/rhor);

//  Roe averages

        rm = sqrt(rhor/rhol);
        rhom = rm*rhol;     // \rho^~

        um = (ul + rm*ur)/(1. + rm);
        hm = (hl + rm*hr)/(1. + rm);
        am = sqrt((gamma - 1.)*(hm - 0.5*um*um));

//  identify wave pattern
        if(um > 0.)
        {
//  Contact wave goes to the right
            Snew = um - am;
            ak = ( (pr - pl) - rhom*am*(ur - ul) )/(2*am*am);
            cflm = (um - am)*dtdx;
//      Test for left sonic rarefaction
            if(fabs(cflm) < Tolson)
            {
                sig = 1.;
                starvals(star, gamma, rhol, ul, el, ak, sig);
                Sml = ul - al;          // \lambda_1^L
                Smr = u_star - a_star;  // \lambda_1^R

//       Left wave is a sonic rarefaction, speed is modified
                if(Sml < 0. && Smr > 0.)
                {
                    Snew = Sml*(Smr - (um - am))/(Smr - Sml);
                }
            }
//       Compute one-sided intercell flux from left side
            if(Snew < 0.)
            {
//              Right eigenvectors
                K[0] = 1.;
                K[1] = um - am;
                K[2] = hm - um*am;

                for(j=0;j<3;j++)
                {
                    F[i][j] = FD[i][j] + Snew*ak*K[j];
                }
            }
            else
            {
                for(j=0;j<3;j++)
                {
                    F[i][j] = FD[i][j];
                }
            }
        }
        else
        {
//       Contact wave goes to the left
            Snew = um + am;
            ak = ( (pr - pl) + rhom*am*(ur - ul) )/(2*am*am);
            cflm = (um + am)*dtdx;

//      Test for right sonic rarefaction
            if(fabs(cflm) < Tolson)
            {
                sig = -1.;
                starvals(star, gamma, rhor, ur, er, ak, sig);
                Sml = u_star + a_star;  // \lambda_3^L
                Smr = ur + ar;          // \lambda_3^R

//       Right wave is a sonic rarefaction, speed is modified
                if(Sml < 0. && Smr > 0.)
                {
                    Snew = Smr*((um + am) - Sml)/(Smr - Sml);
                }
            }
//       Compute one-sided intercell flux from right side
            if(Snew < 0.)
            {
//              Right eigenvectors
                K[0] = 1.;
                K[1] = um + am;
                K[2] = hm + um*am;

                for(j=0;j<3;j++)
                {
                    F[i][j] = FD[i+1][j] -Snew*ak*K[j];
                }
            }
            else
            {
                for(j=0;j<3;j++)
                {
                    F[i][j] = FD[i][j];
                }
            }
        }

    }
};

void starvals(double *star, double gamma, double rhok, double uk,
              double ek, double ak, double sig)
{
    double rho_star, p_star;

//  Roe-Averaged States
    rho_star = rhok + sig*ak;
    star[3] = (rhok*uk + sig*ak*(star[0] - sig*star[2]))/rho_star;   //u_star
    p_star = (gamma - 1.)*(ek + sig*ak*(star[1] - sig*star[0]*star[2]) -
         0.5*rho_star*star[3]*star[3]);
    star[4] = sqrt(gamma*p_star/rho_star);  //a_star_k

};

void HLL(double W[][3], double F[][3], double CS[][3], double *dbl,
         int *p_int)
{
    int i, j, n;
    double gamma, rhol, ul, pl, rhor, ur, pr, al, ar, pm, um,
           Sl, Sr, HLLF, g[8], FD[p_int[0]+2][3], var[6];

    n = p_int[0] + 2;

    gamma = dbl[8];

    gammas(g, gamma);

    for(i=0;i<n;i++)
    {
        if( i==0 || i==(n-1) )
        {
            CS[i][0] = W[i][0];
            CS[i][1] = W[i][0]*W[i][1];
            CS[i][2] = 0.5*( W[i][0]*W[i][1]*W[i][1]) +
                       W[i][2]/(gamma - 1.);
        }
        FD[i][0] = CS[i][1];
        FD[i][1] = CS[i][1]*W[i][1] + W[i][2];
        FD[i][2] = W[i][1]*( CS[i][2] + W[i][2] );
    }
    for(i=0;i<(n-1);i++)
    {
        rhol = W[i][0];
        ul = W[i][1];
        pl = W[i][2];
        rhor = W[i+1][0];
        ur = W[i+1][1];
        pr = W[i+1][2];

        var[0] = rhol;
        var[1] = ul;
        var[2] = pl;
        var[3] = rhor;
        var[4] = ur;
        var[5] = pr;

        al = sqrt(gamma*pl/rhol);
        ar = sqrt(gamma*pr/rhor);

//      Obtain p_* and u_*
        guessp(var, dbl, g, &pm, &um);

//      Estimate speeds S_L and S_R  Toro et al.
        if(pm <= pl)
            Sl = ul - al;
        else
            Sl = ul - al*sqrt(1. + g[1]*(pm/pl - 1.));

        if(pm <= pr)
            Sr = ur + ar;
        else
            Sr = ur + ar*sqrt(1. + g[1]*(pm/pr - 1.));

//      Calculate HLL intercell flux
        if(Sl >= 0.)
        {
            for(j=0;j<3;j++)
            {
                F[i][j] = FD[i][j]; //F_L
            }
        }

        if( (Sl <= 0.) && (Sr >= 0.) )
        {
            for(j=0;j<3;j++)
            {
                HLLF = Sr*FD[i][j] - Sl*FD[i+1][j] +
                       Sl*Sr*(CS[i+1][j] - CS[i][j]);
                F[i][j] = HLLF/(Sr - Sl); // F^{hll}
            }
        }

        if(Sr <= 0.)
        {
            for(j=0;j<3;j++)
            {
                F[i][j] = FD[i+1][j]; //F_R
            }
        }

    }

};

void HLLC(double W[][3], double F[][3], double CS[][3], double *dbl,
         int *p_int)
{
    int i, j, n;
    double gamma, rhol, ul, pl, rhor, ur, pr, al, ar, pm, um,
           Sl, Sr, Sm, e, CSK[3], g[8], FD[p_int[0]+2][3], var[6];

    n = p_int[0] + 2;

    gamma = dbl[8];

    gammas(g, gamma);

    for(i=0;i<n;i++)
    {
        if( i==0 || i==(n-1) )
        {
            CS[i][0] = W[i][0];
            CS[i][1] = W[i][0]*W[i][1];
            CS[i][2] = 0.5*( W[i][0]*W[i][1]*W[i][1]) +
                       W[i][2]/(gamma - 1.);
        }
        FD[i][0] = CS[i][1];
        FD[i][1] = CS[i][1]*W[i][1] + W[i][2];
        FD[i][2] = W[i][1]*( CS[i][2] + W[i][2] );
    }
    for(i=0;i<(n-1);i++)
    {
        rhol = W[i][0];
        ul = W[i][1];
        pl = W[i][2];
        rhor = W[i+1][0];
        ur = W[i+1][1];
        pr = W[i+1][2];

        var[0] = rhol;
        var[1] = ul;
        var[2] = pl;
        var[3] = rhor;
        var[4] = ur;
        var[5] = pr;

        al = sqrt(gamma*pl/rhol);
        ar = sqrt(gamma*pr/rhor);

//      Obtain p_* and u_*
        guessp(var, dbl, g, &pm, &um);

//      Estimate speeds S_L and S_R  Toro et al.
        if(pm <= pl)
            Sl = ul - al;
        else
            Sl = ul - al*sqrt(1. + g[1]*(pm/pl - 1.));

        Sm = um;

        if(pm <= pr)
            Sr = ur + ar;
        else
            Sr = ur + ar*sqrt(1. + g[1]*(pm/pr - 1.));

//      Calculate HLLC intercell flux
        if(Sl >= 0.)
        {
            for(j=0;j<3;j++)
            {
                F[i][j] = FD[i][j]; //F_L
            }
        }

        if( (Sl <= 0.) && (Sr >= 0.) )
        {
//      Subsonic flow
            if(Sm >= 0.)
            {
//              Subsonic flow to the right
                e = CS[i][2]/rhol + (Sm - ul)*
                    (Sm + pl/(rhol*(Sl - ul)));
                CSK[0] = rhol*(Sl - ul)/(Sl - Sm);
                CSK[1] = CSK[0]*Sm;
                CSK[2] = CSK[0]*e;

                for(j=0;j<3;j++)
                {
                    F[i][j] = FD[i][j] + Sl*(CSK[j] - CS[i][j]);  // F_{*L}
                }
            }
            else
            {
//              Subsonic flow to the left
                e = CS[i+1][2]/rhor + (Sm - ur)*
                    (Sm + pr/(rhor*(Sr - ur)));
                CSK[0] = rhor*(Sr - ur)/(Sr - Sm);
                CSK[1] = CSK[0]*Sm;
                CSK[2] = CSK[0]*e;

                for(j=0;j<3;j++)
                {
                    F[i][j] = FD[i+1][j] + Sr*(CSK[j] - CS[i+1][j]);  // F_{*R}
                }
            }
        }

        if(Sr <= 0.)
        {
            for(j=0;j<3;j++)
            {
                F[i][j] = FD[i+1][j]; //F_R
            }
        }

    }

};

void Rusanov(double W[][3], double F[][3], double CS[][3], double *dbl,
         int *p_int)
{
    int i, j, n;
    double gamma, rhol, ul, pl, rhor, ur, pr, al, ar, pm, um,
           Sl, Sr, Splus, RusF, g[8], FD[p_int[0]+2][3], var[6];

    n = p_int[0] + 2;

    gamma = dbl[8];

    gammas(g, gamma);

    for(i=0;i<n;i++)
    {
        if( i==0 || i==(n-1) )
        {
            CS[i][0] = W[i][0];
            CS[i][1] = W[i][0]*W[i][1];
            CS[i][2] = 0.5*( W[i][0]*W[i][1]*W[i][1]) +
                       W[i][2]/(gamma - 1.);
        }
        FD[i][0] = CS[i][1];
        FD[i][1] = CS[i][1]*W[i][1] + W[i][2];
        FD[i][2] = W[i][1]*( CS[i][2] + W[i][2] );
    }
    for(i=0;i<(n-1);i++)
    {
        rhol = W[i][0];
        ul = W[i][1];
        pl = W[i][2];
        rhor = W[i+1][0];
        ur = W[i+1][1];
        pr = W[i+1][2];

        var[0] = rhol;
        var[1] = ul;
        var[2] = pl;
        var[3] = rhor;
        var[4] = ur;
        var[5] = pr;

        al = sqrt(gamma*pl/rhol);
        ar = sqrt(gamma*pr/rhor);

//      Obtain p_* and u_*
        guessp(var, dbl, g, &pm, &um);

//      Estimate speeds S_L and S_R  Toro et al.
        if(pm <= pl)
            Sl = ul - al;
        else
            Sl = ul - al*sqrt(1. + g[1]*(pm/pl - 1.));

        if(pm <= pr)
            Sr = ur + ar;
        else
            Sr = ur + ar*sqrt(1. + g[1]*(pm/pr - 1.));

//      Calculate Rusanov intercell flux
        if(Sl >= 0.)
        {
            for(j=0;j<3;j++)
            {
                F[i][j] = FD[i][j]; //F_L
            }
        }

        Splus = fmax(fabs(Sl), fabs(Sr));

        if( (Sl <= 0.) && (Sr >= 0.) )
        {
            for(j=0;j<3;j++)
            {
                RusF = 0.5*(FD[i][j] + FD[i+1][j]);
                F[i][j] = RusF + 0.5*Splus*(CS[i][j] - CS[i+1][j]); // F_{RUS}
            }
        }

        if(Sr <= 0.)
        {
            for(j=0;j<3;j++)
            {
                F[i][j] = FD[i+1][j]; //F_R
            }
        }

    }

};

void update(double W[][3], double CS[][3], double F[][3],
            double dtdx, double *dbl, int *p_int)
{
    int i, j, n;
    double gamma;

    n = p_int[0] + 2;

    gamma = dbl[8];

    for(i=1;i<(n-1);i++)
    {
        for(j=0;j<3;j++)
        {
            CS[i][j] = CS[i][j] + dtdx*( F[i-1][j] - F[i][j] );
        }
        W[i][0] = CS[i][0];
        W[i][1] = CS[i][1]/W[i][0];
        W[i][2] = (gamma - 1.)*( CS[i][2] - 0.5*CS[i][1]*W[i][1] );

    }

};

void exactRP(double *var, double *star, double *g,
             double *dbl, int *p_int)
{

    star_pu(star, var, dbl, g);
    rhostar(star, var, dbl, g);

};

void star_pu(double *star, double *var, double *dbl, double *g)
{
    int i;
    double change, ul, ur, fL, dfL, fR, dfR, tol, p_start, p_old, um;

    ul = var[1];
    ur = var[4];
    tol = dbl[12];

    guessp(var, dbl, g, &p_start, &um);

    p_old = p_start;

    for(i = 0;i < 20;i++)
    {
        fL = f(var, dbl, g, p_old, 0);
        dfL = df(var, dbl, g, p_old, 0);
        fR = f(var, dbl, g, p_old, 1);
        dfR = df(var, dbl, g, p_old, 1);
        p_start = p_old - (fL + fR + ur - ul)/(dfL + dfR);
        change = 2.*fabs( (p_start - p_old)/(p_start + p_old) );

        if(change < tol) break;
        p_old = p_start;
    }

    if(change > tol) printf("Iteracoes no Newton-Raphson "
            "insuficientes\n");
    star[0] = p_start;
    star[1] = 0.5*( ur + ul + fR - fL );

};

// Solve /rho_L^* e /rho_R^*
void rhostar(double *star, double *var, double *dbl, double *g)
{
    double rhol, pl, rhor, pr, pm, pmL, pmR;

    rhol = var[0];
    pl = var[2];
    rhor = var[3];
    pr = var[5];

    pm = star[0];

    //rholm -- /rho_L^*
    if(pm <= pl)
        star[2] = rhol*pow( (pm/pl), (1./(g[7]+1.)) );
    else
    {
        pmL = pm/pl;
        star[2] = rhol*(pmL + g[5])/(pmL*g[5] + 1.);
    }
    //rhorm -- /rho_R^*
    if(pm <= pr)
        star[3] = rhor*pow( (pm/pr), (1./(g[7]+1.)) );
    else
    {
        pmR = pm/pr;
        star[3] = rhor*(pmR + g[5])/(pmR*g[5] + 1.);
    }

};

void guessp(double *var, double *dbl, double *g, double *pm,
            double *um)
{
    double gamma, rhol, ul, pl, rhor, ur, pr, aL, aR,
           aup, ppv, p_min, p_max, qmax, pq, ptL, ptR,
           geL, geR, quser = 2.;

    rhol = var[0];
    ul = var[1];
    pl = var[2];
    rhor = var[3];
    ur = var[4];
    pr = var[5];
    gamma = dbl[8];

    aL = sqrt(gamma*pl/rhol);
    aR = sqrt(gamma*pr/rhor);

//  Compute guess pressure from PVRS Riemann solver
    aup = 0.25*(rhol + rhor)*(aL + aR);
    ppv = 0.5*( (pl + pr) + (ul - ur)*aup );
    ppv = fmax(0., ppv);
    p_min = fmin(pl, pr);
    p_max = fmax(pl, pr);
    qmax = p_max/p_min;

    if( (qmax <= quser) && ((p_min <= ppv) && (ppv <= p_max)) )
    {
//      Select PVRS Riemann solver
        (*pm) = ppv;
        (*um) = 0.5*(ul + ur) + 0.5*(pl - pr)/aup;
    }
    else
    {
        if(ppv < p_min)
        {
//          Select Two-Rarefaction Riemann solver
            pq = pow( (pl/pr), g[0] );
            (*um) = (pq*ul/aL + ur/aR + g[3]*(pq - 1.))/(pq/aL + 1./aR);
            ptL = 1. + g[6]*(ul - (*um))/aL;
            ptR = 1. + g[6]*((*um) - ur)/aR;
            (*pm) = 0.5*(pl*pow(ptL, g[2]) + pr*pow(ptR, g[2]));
        }
        else
        {
//          Select Two-Shock Riemann solver with PVRS as estimate
            geL = sqrt((g[4]/rhol)/(g[5]*pl + ppv));
            geR = sqrt((g[4]/rhor)/(g[5]*pr + ppv));
            (*pm) = (geL*pl + geR*pr - (ur - ul))/(geL + geR);
            (*um) = 0.5*(ul + ur) + 0.5*(geR*((*pm) - pr) -
                    geL*((*pm) - pl));
        }
    }

};

double f(double *var, double *dbl, double *g, double p_old, int i)
{
    double gamma, a, p, rho, A, B;

    gamma = dbl[8];
    p = var[2 + 3*i];
    rho = var[0 + 3*i];

    a = sqrt(gamma*p/rho);

    if(p_old < p)
    return( g[3]*a*( pow((p_old/p), g[0]) - 1. ) );
    else
    {
        A = g[4]/rho;
        B = g[5]*p;
        return( (p_old - p)*sqrt(A/(B + p_old)) );
    }
};

double df(double *var, double *dbl, double *g, double p_old, int i)
{
    double gamma, a, p, rho, A, B;

    gamma = dbl[8];
    p = var[2 + 3*i];
    rho = var[0 + 3*i];

    a = sqrt(gamma*p/rho);

    if(p_old < p)
    return( (1./(a*rho))*pow((p_old/p), -g[1]) );
    else
    {
        A = g[4]/rho;
        B = g[5]*p;
        return( ( 1. - 0.5*(p_old - p)/(B + p_old) )*
               sqrt(A/(B + p_old)) );
    }
};

void flux(double F[][3], double *var, double *dbl, int m)
{
    double gamma, rho, u, p, e;

    rho = var[0];
    u = var[1];
    p = var[2];

    gamma = dbl[8];

	F[m][0] = rho*u;
    F[m][1] = rho*u*u + p;
    e = 0.5*u*u*rho + p/(gamma - 1.);
    F[m][2] = u*(e + p);

};

void sample(double *var, double *star, double *g,
            double *dbl, double S, int m)
{
    double rho, u, p, rhol, ul, pl, rhor, ur, pr, pm, um,
           rholm, rhorm, a, aL, aR, gamma, SHL, STL, SHR, STR,
           amL, amR, pmL, pmR, SL, SR;

    rhol = var[0];
    ul = var[1];
    pl = var[2];
    rhor = var[3];
    ur = var[4];
    pr = var[5];
    gamma = dbl[8];

    pm = star[0];
    um = star[1];
    rholm = star[2];
    rhorm = star[3];

    aL = sqrt(gamma*pl/rhol);
    aR = sqrt(gamma*pr/rhor);

    if(S <= um)
    {
//Sampling point lies to the left of the contact discontinuity
        if(pm <= pl)
        {
//LEFT RAREFACTION
            SHL = ul - aL;
            if(S <= SHL)
            {
                //Sampled point is left data state
                rho = rhol;
                u = ul;
                p = pl;
            }
            else
            {
                amL = aL*pow( (pm/pl), g[0] );
                STL = um - amL;
                if(S > STL)
                {
                    //Sampled point is Star Left state
                    rho = rholm;
                    u = um;
                    p = pm;
                }
                else
                {
                    //Sampled point is inside left fan
                    a = g[4]*( aL + g[6]*(ul - S) );
                    rho = rhol*pow( (a/aL), g[3] );
                    u = g[4]*( aL + g[6]*ul + S );
                    p = pl*pow( (a/aL), g[2] );
                }
            }
        }
        else
        {
//LEFT SHOCK
            pmL = pm/pl;
            SL = ul - aL*sqrt(g[1]*pmL + g[0]);
            if(S <= SL)
            {
                //Sampled point is left data state
                rho = rhol;
                u = ul;
                p = pl;
            }
            else
            {
                //Sampled point is Star Left state
                rho = rholm;
                u = um;
                p = pm;
            }
        }
    }
    else
    {
//Sampling point lies to the right of the contact discontinuity
        if(pm > pr)
        {
//RIGHT SHOCK
            pmR = pm/pr;
            SR = ur + aR*sqrt(g[1]*pmR + g[0]);
            if(S >= SR)
            {
                //Sampled point is right data state
                rho = rhor;
                u = ur;
                p = pr;
            }
            else
            {
                //Sampled point is Star Right state
                rho = rhorm;
                u = um;
                p = pm;
            }
        }
        else
        {
//RIGHT RAREFACTION
            SHR = ur + aR;
            if(S >= SHR)
            {
                //Sampled point is right data state
                rho = rhor;
                u = ur;
                p = pr;
            }
            else
            {
                amR = aR*pow( (pm/pr), g[0] );
                STR = um + amR;
                if(S <= STR)
                {
                    //Sampled point is Star Right state
                    rho = rhorm;
                    u = um;
                    p = pm;
                }
                else
                {
                    //Sampled point is inside left fan
                    a = g[4]*( aR - g[6]*(ur - S) );
                    rho = rhor*pow( (a/aR), g[3] );
                    u = g[4]*( -aR + g[6]*ur + S );
                    p = pr*pow( (a/aR), g[2] );
                }
            }
        }
    }

    var[0] = rho;
    var[1] = u;
    var[2] = p;

};

void time_msec(struct timespec ti, struct timespec tf)
{
    int ts;
    double tempo, tms;
    tempo = (tf.tv_sec - ti.tv_sec) * 1e3 +
            (tf.tv_nsec - ti.tv_nsec) / 1e6;
//  display time in seconds and milisseconds
    if(tempo > 1e3)
    {
        ts = tempo / 1e3;
        tms = tempo - ts*1e3;
        fprintf(stdout, "Elapsed time: %d s %g ms\n", ts, tms);
    }
    else //time in milisseconds
        fprintf(stdout, "Elapsed time: %g ms\n", tempo);

};
