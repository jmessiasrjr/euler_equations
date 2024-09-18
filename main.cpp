#include "mainwindow.h"
#include <QApplication>
#include <vector>

#define sign(x) ((x >= 0.) ? " " : "-" )

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    return a.exec();
}

void euler(const std::vector<double> &values, const std::vector<int> &parameters, double& et,
           double& t, int& j, int& err)
{
    int i;
    double time, elpst;
    struct timespec ti, tf;

    clock_gettime(CLOCK_MONOTONIC_RAW, &ti);

    vacuum(values, err);

    if(err == -1) return;

    if(parameters.back() == 1) exact(values, parameters);

    godunov(values, parameters, time, i);

    gnuplot(values,parameters[4]);

    clock_gettime(CLOCK_MONOTONIC_RAW, &tf);

    elpst = (tf.tv_sec - ti.tv_sec) +
            (tf.tv_nsec - ti.tv_nsec) / 1e9;

    j = i;
    t = time;
    et = elpst;
}

void gammas(std::vector<double>& g, const double gamma)
{
    g[0] = (gamma - 1.)/(2.*gamma);
    g[1] = (gamma + 1.)/(2.*gamma);
    g[2] = (2.*gamma)/(gamma - 1.);
    g[3] = 2./(gamma - 1.);
    g[4] = 2./(gamma + 1.);
    g[5] = (gamma - 1.)/(gamma + 1.);
    g[6] = (gamma - 1.)/2.;
    g[7] = gamma - 1.;
}

void vacuum(const std::vector<double>& values, int& err)
{
    double gamma, rhol, ul, pl, rhor, ur, pr, al, ar;

    rhol = values[0];
    ul = values[1];
    pl = values[2];
    rhor = values[3];
    ur = values[4];
    pr = values[5];
    gamma = values[9];

    al = sqrt(gamma*pl/rhol);
    ar = sqrt(gamma*pr/rhor);

    if((2./(gamma - 1.))*(al + ar) <= (ur - ul))
    {
        err = -1;
        return;
    }

    err = 0;
}

void exact(const std::vector<double> &values, const std::vector<int> &parameters)
{
    double X, RHO, U, P, E, rho, u, p, dx, xi, xf, PD, gamma, TOUT, S, e;
    FILE *exato;

    exato = fopen("exact.out","w");

    const int n = parameters[0];

    std::vector<double> var(6), g(8), star(5), x(10*n);

    xi = values[6];
    xf = values[7];
    PD = values[8];
    gamma = values[9];
    TOUT = values[10];

    var[0] = values[0];
    var[1] = values[1];
    var[2] = values[2];
    var[3] = values[3];
    var[4] = values[4];
    var[5] = values[5];

    gammas(g, gamma);

    exactRP(var, star, g, values);

    dx = (xf - xi)/(10*n);

    for(int i = 0;i < (10*n);i++)
    {
        x[i] = xi + (i + 0.5)*dx;
        S = ( x[i] - PD )/TOUT;

        sample(var, star, g, values, S);

        X = x[i] - (int)x[i];
        fprintf(exato,"%s%d.%06d  ",sign(x[i]), abs((int)x[i]), abs((int)(1e6*X)));

        rho = var[0]; RHO = rho - (int)rho;
        u = var[1]; U = u - (int)(u);
        p = var[2]; P = p - (int)p;
        e = p/( (gamma - 1)*rho ); E = e - (int)e;

        fprintf( exato,"%s%d.%06d  %s%d.%06d  %s%d.%06d  %s%d.%06d  \n",
                       sign(rho), abs((int)rho), abs((int)(1e6*RHO)),
                       sign(u), abs((int)u), abs((int)(1e6*U)),
                       sign(p), abs((int)p), abs((int)(1e6*P)),
                       sign(e), abs((int)e), abs((int)(1e6*E)) );
    }

    fclose(exato);

}

void godunov(const std::vector<double>& values, const std::vector<int>& parameters, double &t, int &j)
{
    int i;
    double X, RHO, U, P, E, rho, u, p, x, dx, e, dt, time = 0.;
    FILE *saida;

    saida = fopen("result.out","w");

    const int &n        = parameters[0] + 2;
    const int &method   = parameters[3];

    const double &xi    = values[6];
    const double &xf    = values[7];
    const double &gamma = values[9];
    const double &TOUT  = values[10];
    const double &MAXT  = values[11];

    dx = (xf - xi)/(n-2);

    std::vector<std::vector<double>> W(n, std::vector<double>(3)),
            CS(n, std::vector<double>(3)), F(n, std::vector<double>(3));

    CI(W, CS, values, parameters);

    for(i=0; i<MAXT; i++)
    {
        boundary(W, parameters);
        dt = cfl(W, values, parameters, time, i);

        switch(method)
        {
            case(1): roe(W, F, CS, values, parameters, dt/dx);
            break;

            case(2): HLL(W, F, CS, values, parameters);
            break;

            case(3): HLLC(W, F, CS, values, parameters);
            break;

            case(4): Rusanov(W, F, CS, values, parameters);
            break;

            default: godRP(W, F, values, parameters);
        }

        update(W, F, CS, values, parameters, dt/dx);

        time += dt;
        if(fabs(time -TOUT) < 1e-6) break;
    }

    j = i;
    t = time;

    for(i=1; i<(n-1); i++)
    {
        rho = W[i][0]; RHO = rho - (int)rho;
        u = W[i][1]; U = u - (int)u;
        p = W[i][2]; P = p - (int)p;

        x = xi + (i - 0.5)*dx; X = x - (int)x;
        e = p/( (gamma - 1)*rho ); E = e - (int)e;

        fprintf(saida,"%s%d.%06d  %s%d.%06d  %s%d.%06d  %s%d.%06d  %s%d.%06d \n",
                      sign(x), abs((int)x), abs((int)(1e6*X)),
                      sign(rho), abs((int)rho), abs((int)(1e6*RHO)),
                      sign(u), abs((int)u), abs((int)(1e6*U)),
                      sign(p), abs((int)p), abs((int)(1e6*P)),
                      sign(e), abs((int)e), abs((int)(1e6*E)));
    }

    fclose(saida);
}

void CI(std::vector<std::vector<double>> &W, std::vector<std::vector<double>> &CS, const std::vector<double> &values,
        const std::vector<int> &parameters)
{
    double x, dx;

    const int n = parameters[0] + 2;

    const double &rhol  = values[0];
    const double &ul    = values[1];
    const double &pl    = values[2];
    const double &rhor  = values[3];
    const double &ur    = values[4];
    const double &pr    = values[5];
    const double &xi    = values[6];
    const double &xf    = values[7];
    const double &PD    = values[8];
    const double &gamma = values[9];

    dx = (xf - xi)/(n-2);

    for(int i=1; i<(n-1); i++)
    {
        x = xi + (i-.5)*dx;
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

}

void boundary(std::vector<std::vector<double>> &W, const std::vector<int> &parameters)
{
    const int &n   = parameters[0] + 2;
    const int &BCL = parameters[1];
    const int &BCR = parameters[2];

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
}

double cfl(const std::vector<std::vector<double>> &W, const std::vector<double> &values,
           const std::vector<int> &parameters, const double time, const int m)
{
    double a, rho, u, p, xi, xf, gamma, TOUT, CFL, dx, dt, S, Smax = 0.;

    const int n = parameters[0] + 2;

    xi = values[6];
    xf = values[7];
    gamma = values[9];
    TOUT = values[10];
    CFL = values[12];

    for(int i=0;i<n;i++)
    {
        rho = W[i][0];
        u = W[i][1];
        p = W[i][2];
        a = sqrt(gamma*p/rho);
        S = fabs(u)+a;
        if(S > Smax) Smax = S;
    }

    dx = (xf - xi)/(n-2);
    dt = CFL*dx/Smax;

    if(m<5) dt = 0.2*dt;
    if( (time+dt) > TOUT ) dt = TOUT - time;

    return(dt);

}

void godRP(std::vector<std::vector<double> > &W, std::vector<std::vector<double> > &F,
           const std::vector<double> &values, const std::vector<int> &parameters)
{
    const int n = parameters[0] + 2;
    const double gamma = values[9];

    std::vector<double> var(6), g(8), star(5);

    gammas(g, gamma);

    for(int i=0; i<(n-1); i++)
    {
        var[0] = W[i][0];
        var[1] = W[i][1];
        var[2] = W[i][2];
        var[3] = W[i+1][0];
        var[4] = W[i+1][1];
        var[5] = W[i+1][2];

        exactRP(var, star, g, values);

        sample(var, star, g, values, 0.);

        flux(F, var, values, i);
    }

}

void roe(std::vector<std::vector<double>> &W, std::vector<std::vector<double>> &F, std::vector<std::vector<double>> &CS,
         const std::vector<double> &values, const std::vector<int> &parameters, double dtdx)
{
    double rhol, ul, pl, el, hl, rhor, ur, pr, er, hr,
           al, ar, rm, rhom, um, hm, am, u_star, a_star, ak,
           sig, cflm, Sml, Smr, Snew;

    const int n      = parameters[0] + 2;
    const int Tolson = 0.1;

    std::vector<double> K(3), star(5);
    std::vector<std::vector<double>> FD(n, std::vector<double>(3));

    const double &gamma = values[9];

    for(int i=0; i<n; i++)
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

    for(int i=0; i<(n-1); i++)
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

        star[0] = um;       // \u^~
        star[1] = hm;       // \H^~
        star[2] = am;       // \a^~

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
                u_star = star[3];
                a_star = star[4];
                Sml = ul - al;          // \lambda_1^L
                Smr = u_star - a_star;  // \lambda_1^R

                if((Sml < 0.) && (Smr > 0.))
                {
//       Left wave is a sonic rarefaction, speed is modified
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

                for(int j=0; j<3; j++)
                {
                    F[i][j] = FD[i][j] + Snew*ak*K[j];
                }
            }
            else
            {
                for(int j=0; j<3; j++)
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
                u_star = star[3];
                a_star = star[4];
                Sml = u_star + a_star;  // \lambda_3^L
                Smr = ur + ar;          // \lambda_3^R

//       Right wave is a sonic rarefaction, speed is modified
                if((Sml < 0.) && (Smr > 0.))
                {
                    Snew = Smr*((um + am) - Sml)/(Smr - Sml);
                }
            }
//       Compute one-sided intercell flux from right side
            if(Snew > 0.)
            {
//              Right eigenvectors
                K[0] = 1.;
                K[1] = um + am;
                K[2] = hm + um*am;

                for(int j=0; j<3; j++)
                {
                    F[i][j] = FD[i+1][j] - Snew*ak*K[j];
                }
            }
            else
            {
                for(int j=0; j<3; j++)
                {
                    F[i][j] = FD[i+1][j];
                }
            }

        }

    }
}

void starvals(std::vector<double> &star, const double gamma, const double rhok, const double uk,
              const double ek, const double ak, const double sig)
{
    double um, hm, am, rho_star, p_star;

    um = star[0];
    hm = star[1];
    am = star[2];

//  Roe-Averaged States
    rho_star = rhok + sig*ak;
    star[3] = (rhok*uk + sig*ak*(um - sig*am))/rho_star;   //u_star
    p_star = (gamma - 1.)*(ek + sig*ak*(hm - sig*um*am) -
         0.5*rho_star*star[3]*star[3]);
    star[4] = sqrt(gamma*p_star/rho_star);  //a_star_k

}

void HLL(std::vector<std::vector<double> > &W, std::vector<std::vector<double> > &F, std::vector<std::vector<double> > &CS,
         const std::vector<double> &values, const std::vector<int> &parameters)
{
    double rhol, ul, pl, rhor, ur, pr, al, ar, pm, um, Sl, Sr, HLLF;

    const int n = parameters[0] + 2;
    const double gamma = values[9];

    std::vector<double> g(8), var(6);
    std::vector<std::vector<double>> FD(n, std::vector<double>(3));

    gammas(g, gamma);

    for(int i=0; i<n; i++)
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
    for(int i=0; i<(n-1); i++)
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
        guessp(var, values, g, pm, um);

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
            for(int j=0; j<3; j++)
            {
                F[i][j] = FD[i][j]; //F_L
            }
        }

        if( (Sl <= 0.) && (Sr >= 0.) )
        {
            for(int j=0; j<3; j++)
            {
                HLLF = Sr*FD[i][j] - Sl*FD[i+1][j] +
                       Sl*Sr*(CS[i+1][j] - CS[i][j]);
                F[i][j] = HLLF/(Sr - Sl); // F^{hll}
            }
        }

        if(Sr <= 0.)
        {
            for(int j=0; j<3; j++)
            {
                F[i][j] = FD[i+1][j]; //F_R
            }
        }

    }

}

void HLLC(std::vector<std::vector<double>> &W, std::vector<std::vector<double>> &F,
          std::vector<std::vector<double>> &CS, const std::vector<double> &values,
          const std::vector<int> &parameters)
{
    double rhol, ul, pl, rhor, ur, pr, al, ar, pm, um, Sl, Sr, Sm, e;

    const int n = parameters[0] + 2;
    const double gamma = values[9];

    std::vector<double> CSK(3), g(8), var(6);
    std::vector<std::vector<double>> FD(n, std::vector<double>(3));

    gammas(g, gamma);

    for(int i=0; i<n; i++)
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
    for(int i=0; i<(n-1); i++)
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
        guessp(var, values, g, pm, um);

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
            for(int j=0; j<3; j++)
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

                for(int j=0;j<3;j++)
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

                for(int j=0; j<3; j++)
                {
                    F[i][j] = FD[i+1][j] + Sr*(CSK[j] - CS[i+1][j]);  // F_{*R}
                }
            }
        }

        if(Sr <= 0.)
        {
            for(int j=0; j<3; j++)
            {
                F[i][j] = FD[i+1][j]; //F_R
            }
        }

    }

}

void Rusanov(std::vector<std::vector<double>> &W, std::vector<std::vector<double>> &F,
             std::vector<std::vector<double>> &CS, const std::vector<double> &values,
             const std::vector<int> &parameters)
{
    double rhol, ul, pl, rhor, ur, pr, al, ar, pm, um, Sl, Sr, Splus, RusF;

    const int n = parameters[0] + 2;
    const double gamma = values[9];

    std::vector<double> g(8), var(6);
    std::vector<std::vector<double>> FD(n, std::vector<double>(3));

    gammas(g, gamma);

    for(int i=0; i<n; i++)
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
    for(int i=0; i<(n-1); i++)
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
        guessp(var, values, g, pm, um);

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
            for(int j=0; j<3; j++)
            {
                F[i][j] = FD[i][j]; //F_L
            }
        }

        Splus = fmax(fabs(Sl), fabs(Sr));

        if( (Sl <= 0.) && (Sr >= 0.) )
        {
            for(int j=0; j<3; j++)
            {
                RusF = 0.5*(FD[i][j] + FD[i+1][j]);
                F[i][j] = RusF + 0.5*Splus*(CS[i][j] - CS[i+1][j]); // F_{RUS}
            }
        }

        if(Sr <= 0.)
        {
            for(int j=0; j<3; j++)
            {
                F[i][j] = FD[i+1][j]; //F_R
            }
        }

    }

}

void update(std::vector<std::vector<double>> &W, std::vector<std::vector<double>> &F,
            std::vector<std::vector<double>> &CS, const std::vector<double> &values,
            const std::vector<int> &parameters, double dtdx)
{
    const int n = parameters[0] + 2;
    const double gamma = values[9];

    for(int i=1; i<(n-1); i++)
    {
        for(int j=0; j<3; j++)
        {
            CS[i][j] = CS[i][j] + dtdx*( F[i-1][j] - F[i][j] );
        }
        W[i][0] = CS[i][0];
        W[i][1] = CS[i][1]/W[i][0];
        W[i][2] = (gamma - 1.)*( CS[i][2] - 0.5*CS[i][1]*W[i][1] );
    }
}

void exactRP(const std::vector<double> &var, std::vector<double> &star,
             const std::vector<double> &g, const std::vector<double> &values)
{
    star_pu(var, star, g, values);
    rhostar(var, star, g);
}

void star_pu(const std::vector<double> &var, std::vector<double> &star,
             const std::vector<double> &g, const std::vector<double> &values)
{
    double change, ul, ur, fL, dfL, fR, dfR, tol, p_start, p_old, um;

    ul = var[1];
    ur = var[4];
    tol = values[13];

    guessp(var, values, g, p_start, um);

    p_old = p_start;

    for(int i=0; i < 20; i++)
    {
        fL = f(var, values, g, p_old, 0);
        dfL = df(var, values, g, p_old, 0);
        fR = f(var, values, g, p_old, 1);
        dfR = df(var, values, g, p_old, 1);
        p_start = p_old - (fL + fR + ur - ul)/(dfL + dfR);
        change = 2.*fabs( (p_start - p_old)/(p_start + p_old) );

        if(change < tol) break;
        p_old = p_start;
    }

    if(change > tol) printf("Iteracoes no Newton-Raphson "
            "insuficientes\n");
    star[0] = p_start;
    star[1] = 0.5*( ur + ul + fR - fL );
}

// Solve /rho_L^* e /rho_R^*
void rhostar(const std::vector<double> &var, std::vector<double> &star, const std::vector<double> &g)
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

}

void guessp(const std::vector<double> &var, const std::vector<double> &values,
            const std::vector<double> &g, double &pm, double &um)
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
    gamma = values[9];

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
        pm = ppv;
        um = 0.5*(ul + ur) + 0.5*(pl - pr)/aup;
    }
    else
    {
        if(ppv < p_min)
        {
//          Select Two-Rarefaction Riemann solver
            pq = pow( (pl/pr), g[0] );
            um = (pq*ul/aL + ur/aR + g[3]*(pq - 1.))/(pq/aL + 1./aR);
            ptL = 1. + g[6]*(ul - um)/aL;
            ptR = 1. + g[6]*(um - ur)/aR;
            pm = 0.5*(pl*pow(ptL, g[2]) + pr*pow(ptR, g[2]));
        }
        else
        {
//          Select Two-Shock Riemann solver with PVRS as estimate
            geL = sqrt((g[4]/rhol)/(g[5]*pl + ppv));
            geR = sqrt((g[4]/rhor)/(g[5]*pr + ppv));
            pm = (geL*pl + geR*pr - (ur - ul))/(geL + geR);
            um = 0.5*(ul + ur) + 0.5*(geR*(pm - pr) -
                    geL*(pm - pl));
        }
    }

}

double f(const std::vector<double> &var, const std::vector<double> &values,
         const std::vector<double> &g, const double p_old, const int i)
{
    double A, B;

    const double gamma = values[9];
    const double p     = var[2 + 3*i];
    const double rho   = var[0 + 3*i];

    const double a = sqrt(gamma*p/rho);

    if(p_old < p)
        return( g[3]*a*( pow((p_old/p), g[0]) - 1. ) );
    else
    {
        A = g[4]/rho;
        B = g[5]*p;
        return( (p_old - p)*sqrt(A/(B + p_old)) );
    }
}

double df(const std::vector<double> &var, const std::vector<double> &values,
          const std::vector<double> &g, const double p_old, const int i)
{
    double A, B;

    const double gamma = values[9];
    const double p     = var[2 + 3*i];
    const double rho   = var[0 + 3*i];

    const double a = sqrt(gamma*p/rho);

    if(p_old < p)
        return( (1./(a*rho))*pow((p_old/p), -g[1]) );
    else
    {
        A = g[4]/rho;
        B = g[5]*p;
        return( ( 1. - 0.5*(p_old - p)/(B + p_old) )*
               sqrt(A/(B + p_old)) );
    }
}

void flux(std::vector<std::vector<double>> &F, const std::vector<double> &var,
          const std::vector<double> &values, const int m)
{
    double gamma, rho, u, p, e;

    rho = var[0];
    u = var[1];
    p = var[2];

    gamma = values[9];

    F[m][0] = rho*u;
    F[m][1] = rho*u*u + p;
    e = 0.5*u*u*rho + p/(gamma - 1.);
    F[m][2] = u*(e + p);

}

void sample(std::vector<double> &var, std::vector<double> &star,
            const std::vector<double> &g, const std::vector<double> &values, const double S)
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
    gamma = values[9];

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

}

void gnuplot(const std::vector<double> &values, const int ex)
{
    FILE *plot;

    const double xi = values[6];
    const double xf = values[7];

    plot = fopen("gnuplot.gp","w");

    fprintf(plot,"set terminal pdfcairo enhanced color dashed size 6 in, 4 in \n");
    fprintf(plot,"set encoding iso_8859_1 \n");
    if(ex==1) fprintf(plot,"datafile1 = 'exact.out' \n");
    fprintf(plot,"datafile2 = 'result.out' \n");
    fprintf(plot,"set format x '%%8.2g' \n");
    fprintf(plot,"set format y '%%12.4g' \n");
    fprintf(plot,"unset key \n");
    fprintf(plot,"set xtics %d.%04d \n",(int)((xf - xi)/4.), (int)(1e4*(xf - xi)/4.));
    fprintf(plot,"set xlabel ('x') \n");
    fprintf(plot,"set output 'data.pdf' \n");
    fprintf(plot,"set multiplot layout 2,2 rowsfirst \n");

    fprintf(plot,"set ylabel ('{/Symbol r}') \n");
    if(ex == 0) fprintf(plot,"plot datafile2 using 1:2 with lines ls 1 lw 2\n");
    else fprintf(plot,"plot datafile1 using 1:2 with lines ls 1 lw 2, \\\n"
            "datafile2 using 1:2 with lines lc 1 dt 2 lw 2 \n");

    fprintf(plot,"set ylabel ('u') \n");
    if(ex == 0) fprintf(plot,"plot datafile2 using 1:3 with lines ls 2 lw 2\n");
    else fprintf(plot,"plot datafile1 using 1:3 with lines ls 2 lw 2, \\\n"
            "datafile2 using 1:3 with lines lc 2 dt 2 lw 2 \n");

    fprintf(plot,"set ylabel ('p') \n");
    if(ex == 0) fprintf(plot,"plot datafile2 using 1:4 with lines ls 7 lw 2\n");
    else fprintf(plot,"plot datafile1 using 1:4 with lines ls 7 lw 2, \\\n"
            "datafile2 using 1:4 with lines lc 7 dt 2 lw 2 \n");

    fprintf(plot,"set ylabel ('e') \n");
    if(ex == 0) fprintf(plot,"plot datafile2 using 1:5 with lines ls 8 lw 2\n");
    else fprintf(plot,"plot datafile1 using 1:5 with lines ls 8 lw 2, \\\n"
            "datafile2 using 1:5 with lines lc 8 dt 2 lw 2 \n");

    fclose(plot);
}
