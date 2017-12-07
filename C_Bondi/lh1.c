/* ==================================================================================== */
/* ====================== Description du programme lh1.c ===================== */
/* ==================================================================================== */
/* Auteurs & modifications                                                              */
/*    - Eliot AYACHE    : 24/10/17                                                      */
/*                                                                                      */
/* Fonctionnement                                                                       */
/*    - This Lagrangian code follozs the adiabatic expansion of a hot spherical bubble  */
/* in a unifor ambient medium. The present setup is adjusted to spherical geometry      */
/* The code can also work in cartesian geometry. To switch go to subroutine inincond    */
/* and set cartesian to true                                                            */
/*                                                                                      */
/* Acknowledgments                                                                      */
/*    - This code is based on the Fortran version included in                           */
/*                  Numerical Methods in astrophysics: An Introduction                  */
/* */

/* ==================================================================================== */
/* =============================== En tetes necessaires =============================== */
/* ==================================================================================== */
#include <stdio.h>
#include <stdlib.h>
#include <math.h> /*utiliser -lm lors de la compilation*/
/*#include <err.h>*/

/* ==================================================================================== */
/* ============================== Constantes symboliques ============================== */
/* ==================================================================================== */

#define DUMPRES     250         // Data dump step
#define NX          1000        // Total number of cells
#define CFLFACTOR   0.1         // Courant factor, must be lower than 1 
#define DFACTOR     (pow(1.5,2))// Diffusion factor, must be greater than 1
#define GAMMA       (4./3.)     // Heat ratio
#define PI          M_PI        // Pi
#define G           1e0         // Gravitational constant
#define R           1.          // Bondi ball radius

/* ==================================================================================== */
/* ================================= Global variables ================================= */
/* ==================================================================================== */

FILE    *fdata  = NULL; // input-data file   
FILE    *fout1  = NULL; // Final model output file   
FILE    *fout2  = NULL; // Incremented output file 

// Parameters:
int     nsteps  = 0;    // Max number of time steps
double  tmax    = 0;    // Max integration time
double  efac    = 0.;   // fireball internal energy / ambient internal energy
int     il      = 0;    // Number of grid cells in the fireball
double  q       = 0.;   // artificial viscosity parameter
double  fmratio = 0.;   // mass ratio ejecta / ambient
int     bondi   = 1;    // switch to include bondi accretion
double  smass   = 1.;   // star mass


// Simulation (12 represents cell edges in staggered grid):
double  t = 0.;                     // time
double  dt, dt12;                   // time steps
double  *fm, *dm, *dm12,            // Mass
        *r,  *dr, *dr12, *rold,     // Radii
        *v,  *a,  *at12, *ak12,     // Volume, a = 4*PI(r**2)
        *u,  *p,  *rho,  *eps,      // Speed, Thermodynamics
        *w,  *wt12, 
        *aux ;


/* ==================================================================================== */
/* ============================= Definitions de fonctions ============================= */
/* ==================================================================================== */
void tmsc();
void initiate();
void multiallocate(int);
void emptybuffer(FILE*);








int main(void)
{
    /* ----------------------------------------------------- */
    /* Variables                                             */
    /* ----------------------------------------------------- */
    // Physical:
    double  ethe, ekin, epot, etot;   // Thermal, Kinetic, Total energies
    double  etot0;              // Initial total energy

    // Practical:
    char    buffer[100];        // buffer used to skip lines
    int     istep = 0;          // Time step
    int     ix;                 // iteration integer

    /* ----------------------------------------------------- */
    /* Program Starts Here                                   */
    /* ----------------------------------------------------- */
    fdata = fopen("lh1.dat", "r");
    fout1 = fopen("lh1.out", "w");
    fout2 = fopen("lh2.out", "w");

    fgets(buffer, 100, fdata);
    fscanf(fdata, "%d", &nsteps);
    emptybuffer(fdata);
    fscanf(fdata, "%lf", &tmax);
    emptybuffer(fdata);
    fscanf(fdata, "%lf", &efac);
    emptybuffer(fdata);
    fscanf(fdata, "%d", &il);
    emptybuffer(fdata);
    fscanf(fdata, "%lf", &q);
    emptybuffer(fdata);
    fscanf(fdata, "%lf", &fmratio);
    emptybuffer(fdata);
    fscanf(fdata, "%d", &bondi);
    emptybuffer(fdata);
    fscanf(fdata, "%lf", &smass);
    emptybuffer(fdata);


    // Initiating Grid:
    initiate();
    if (bondi) printf("Using Bondi accretion!\n");

    // Begin loop over time-steps:
    for (istep = 0; istep < nsteps; ++istep)
    {
        // determine time-step:
        tmsc();



        /* ----------------------------------------------------- */
        /* Uptating variables                                    */
        /* ----------------------------------------------------- */

        // update velocities:
        if (!bondi)
        {
            for (ix = 1; ix < NX; ++ix)
                u[ix] = u[ix] - a[ix]*(p[ix] - p[ix-1])*dt/dm[ix]
                            -0.5 * (w[ix]   * (3.*ak12[ix]   - a[ix])
                                  - w[ix-1] * (3.*ak12[ix-1] - a[ix]))
                            * dt/dm[ix];
            u[0] = 0.;  // boundary condition
        } else {
            for (ix = 1; ix < NX; ++ix)
                u[ix] = u[ix] - a[ix]*(p[ix] - p[ix-1])*dt/dm[ix]
                            -0.5 * (w[ix]   * (3.*ak12[ix]   - a[ix])
                                  - w[ix-1] * (3.*ak12[ix-1] - a[ix]))
                            * dt/dm[ix]
                            - G * smass * dt / (r[ix]*r[ix]);
            for (ix = NX-1; ix >= 0; --ix)
                if (r[ix] < R/10000.) u[ix] = u[ix+1];
        }




        // update radii, surfaces, and volumes:
        for (ix = 0; ix < NX;   ++ix) rold[ix] = r[ix];
        for (ix = 0; ix < NX;   ++ix) r[ix]    = rold[ix] + u[ix]*dt12;
        for (ix = 0; ix < NX-1; ++ix) dr12[ix] = r[ix+1] - r[ix];
        for (ix = 0; ix < NX  ; ++ix)
        {
            at12[ix] = 4. * PI * pow(0.5 * (r[ix]+rold[ix]), 2);
            a[ix]    = 4. * PI * pow(r[ix], 2);
            v[ix]    = 4./3. * PI * pow(r[ix], 3);
        }

        for (ix = 0; ix < NX-1  ; ++ix)
            ak12[ix] = 0.5 * (at12[ix+1] + at12[ix]);
        


        // update densities:
        for (ix = 0; ix < NX-1; ++ix) rho[ix] = dm12[ix]/(v[ix+1]-v[ix]);
        rho[NX-1] = rho[NX-2];


        // artificial viscosity:
        for (ix = 0; ix < NX-1; ++ix)
        {
            w[ix] = - pow(q, 2) * rho[ix] * fabs(u[ix+1] - u[ix])
                        *(u[ix+1] * (1.-at12[ix+1]/3./ak12[ix])
                        - u[ix]   * (1.-at12[ix]  /3./ak12[ix]));
        }

        for (ix = 0; ix < NX-1; ++ix)
        {
            if (u[ix+1] > u[ix]) w[ix]=0;
        }


        // update internal energies and pressures:
        for (ix = 0; ix < NX-1; ++ix)
            aux[ix] = eps[ix] - p[ix]
                        * (  at12[ix+1]*u[ix+1]
                           - at12[ix]  *u[ix]  )*dt12/dm12[ix];

        for (ix = 0; ix < NX-1; ++ix)
            p[ix]   = 0.5 * (p[ix] + (GAMMA-1.) * rho[ix] * aux[ix]);

        for (ix = 0; ix < NX-1; ++ix)
            eps[ix] = eps[ix] - p[ix]
                        * (  at12[ix+1]*u[ix+1]
                           - at12[ix]  *u[ix]  )*dt12/dm12[ix];


        // Contribution from artificial viscosity:
        for (ix = 0; ix < NX-1   ; ++ix)
        {
            eps[ix] = eps[ix] - 0.5*w[ix] * dt12/dm12[ix]
                        * (u[ix+1] * (3.*ak12[ix] - at12[ix+1])
                         - u[ix]   * (3.*ak12[ix] - at12[ix]  ));
        }

        for (ix = 0; ix < NX-1; ++ix)
            p[ix] = (GAMMA-1.) * rho[ix] * eps[ix];

        p[NX-1]   = p[NX-2];
        eps[NX-1] = eps[NX-2];


        // End time step:
        t = t + dt12;


        /* ----------------------------------------------------- */
        /* Check Energy conservation                             */
        /* ----------------------------------------------------- */

        ethe = 0.;
        ekin = 0.;
        epot = 0.;
        for (ix = 1; ix < NX-1; ++ix)
        {
            ethe = ethe + eps[ix]*dm[ix];
            ekin = ekin + 0.5 * pow(0.5 * (u[ix+1]+u[ix]), 2)*dm[ix];
            epot = epot - dm[ix] * G * smass / r[ix];
            // printf("%lf %lf %lf %lf\n", ethe, ekin, dm[ix], eps[ix]);
        }

        etot = ethe + ekin + epot;
        if (istep == 1) 
        {
            etot0 = etot;
        }
        etot = etot/etot0;
        ethe = ethe/etot0;
        ekin = ekin/etot0;
        epot = epot/etot0;

        if (istep%10 == 0.)
            fprintf(fout2, "%d %lf %lf %lf %lf %lf\n",
                        istep, t, etot, ethe, ekin, epot);
        if (istep%100 == 0.)
            printf("istep = %d, \
                            t = %lf, \
                            etot = %lf, \
                            ethe = %lf, \
                            ekin = %lf, \
                            epot = %lf\n",
                        istep, t, etot, ethe, ekin, epot);

        if (istep%DUMPRES == 0.)
            for (ix = 0; ix < NX-1; ++ix)
                fprintf(fout1, "%d %lf %d %lf %lf %lf %lf %lf %lf %lf\n", 
                        istep, t, ix, fm[ix], r[ix], u[ix],
                        rho[ix], p[ix], eps[ix], w[ix]);
    }


    /* ----------------------------------------------------- */
    /* Closing files                                         */
    /* ----------------------------------------------------- */
    fclose(fdata);
    fclose(fout1);
    fclose(fout2);
    fdata = NULL;
    fout1 = NULL;
    fout2 = NULL;

    return 0;
}





/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* This function computes the adequate time-scale based on cell and diffusion limits    */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void tmsc()
{
    /* ----------------------------------------------------- */
    /* Variables                                             */
    /* ----------------------------------------------------- */
    double  dtc = 1.e30;    // Include description
    double  dtd = 1.e-30;   // Include description 
    int     ix;             // iteration integer


    for (ix = 0; ix < NX-1; ++ix)
        dtc = fmin(dtc, 
                dr12[ix] / (fabs(u[ix]) + sqrt(GAMMA * eps[ix])));


    dtc = CFLFACTOR * dtc;
    if (t+dtc > tmax) dtc = tmax-t;

    // diffusion limit:

    for (ix = 0; ix < NX-1; ++ix)
        dtd = fmax(dtd, 
                fabs(at12[ix+1] * u[ix+1] - at12[ix] * u[ix])
                / (v[ix+1] - v[ix]));


    dtd  = 0.5/dtd/DFACTOR;
    dtc  = fmin(dtc,dtd);
    dt   = 0.5 * (dt12+dtc);
    dt12 = dtc;

    return;
}






/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* This function sets-up the first step of the simulation                               */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void initiate()
{
    /* ----------------------------------------------------- */
    /* Variables                                             */
    /* ----------------------------------------------------- */
    double  rhoej, dmej;    // ejecta mass density, shell mass
    double  rhoamb, dmamb;  // ambiant mass density
    int     ix;             // iteration integer



    /* ----------------------------------------------------- */
    /* Actions                                               */
    /* ----------------------------------------------------- */
    // Allocating the grid:
    multiallocate(NX);


    // Filling up ejecta medium:
    rhoej = ( pow((float)NX/(float)il, 3) - 1.) * fmratio ;
    dmej  = 4./3. * PI * pow((float)il/(float)NX, 3) * rhoej / (float) il ;

    for (ix = 0; ix < il; ++ix)
    {
        dm12[ix] = dmej;    // Original dm definition
        // dm12[ix] = 4./3. * PI * rhoej * (float)(pow(ix+1,3.) - pow(ix,3.)) / pow((float)il, 3);
        // Using radius different from 1.
        // dm12[ix] = pow(R,3.) * 4./3. * PI * rhoej * 
        //     (float)(pow(ix+1,3.) - pow(ix,3.)) / pow((float)il, 3);
        rho [ix] = rhoej;
        eps [ix] = efac;
        p   [ix] = (GAMMA - 1.) * rho[ix] * eps[ix] ;   // eq. (6.4)
        u   [ix] = 0.;
    }


    // Filling up the ambient medium:
    rhoamb = 1.;    // normalised mass density
    dmamb = 4./3. * PI 
        * (1. - pow((float)il/(float)NX, 3)) / (float)(NX - il)
        * rhoamb ;

    for (ix = il; ix < NX; ++ix)
    {
        dm12[ix] = dmamb;
        // dm12[ix] = 4./3. * PI * rhoamb * (float)(pow(ix+1,3.) - pow(ix,3.)) / pow((float)NX, 3);
        // Using radius different from 1.
        // dm12[ix] = pow(R,3.) * 4./3. * PI * rhoamb * 
        //     (float)(pow(ix+1,3.) - pow(ix,3.)) / pow((float)NX, 3);
        rho [ix] = rhoamb;
        eps [ix] = 1.;
        p   [ix] = (GAMMA - 1.) * rho[ix] * eps[ix] ;   // eq. (6.4)
        u   [ix] = 0.;
    }

    // If bondi accretion we fill up the whole medium with
    // the external properties:
    if (bondi)
    {
        dmamb = 4./3. * PI / (float)(NX) * rhoamb ;
        for (ix = 0; ix < NX; ++ix)
        {
            // dm12[ix] = dmamb;
            // dm12[ix] = 4./3. * PI * rhoamb * (float)(pow(ix+1,3.) - pow(ix,3.)) / pow((float)NX, 3);
            // Using radius different from 1.
            dm12[ix] = pow(R,3.) * 4./3. * PI * rhoamb * 
                (float)(pow(ix+1,3.) - pow(ix,3.)) / pow((float)NX, 3);
            rho [ix] = rhoamb;
            eps [ix] = 1.;
            p   [ix] = (GAMMA - 1.) * rho[ix] * eps[ix] ;   // eq. (6.4)
            u   [ix] = 0.;
        }        
    }

    fm[0] = dm12[0];
    for (ix = 1; ix < NX; ++ix)
    {
        fm[ix] = fm[ix-1] + dm12[ix-1];
    }

    for (ix = 1; ix < NX; ++ix)
    {
        dm[ix] = 0.5 * (dm12[ix] + dm[ix-1]);
    }

    r[0] = 0.;
    v[0] = 0.;


    for (ix = 1; ix < NX; ++ix)
    {
        v[ix] = v[ix-1] + dm12[ix-1] / rho[ix-1];
        // spherical case only (update to cartesian later):
        r[ix] = pow(v[ix] /(4./3.*PI), 1./3.);
        a[ix] = 4.*PI*pow(r[ix],2);
    }

    for (ix = 0; ix < NX-1; ++ix)
    {
        dr12[ix] = r[ix+1] - r[ix];
    }

    // Initializing simulation:
    t = 0.;

    return;
}





/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* This function allocates memory for all the global arrays in the program              */
/*  -> The argument is the size of the array                                            */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void multiallocate(int n)
{
    fm      = calloc(n, sizeof(double));
    dm      = calloc(n, sizeof(double));
    dm12    = calloc(n, sizeof(double));
    r       = calloc(n, sizeof(double));
    dr      = calloc(n, sizeof(double));
    dr12    = calloc(n, sizeof(double));
    rold    = calloc(n, sizeof(double));
    v       = calloc(n, sizeof(double));
    a       = calloc(n, sizeof(double));
    at12    = calloc(n, sizeof(double));
    ak12    = calloc(n, sizeof(double));
    u       = calloc(n, sizeof(double));
    p       = calloc(n, sizeof(double));
    rho     = calloc(n, sizeof(double));
    eps     = calloc(n, sizeof(double));
    w       = calloc(n, sizeof(double));
    wt12    = calloc(n, sizeof(double));
    aux     = calloc(n, sizeof(double));
    return;
}






/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++













++++++++++++++++++++++++++++++ */
/* This function empties the buffer associated with a stream                            */
/*    -> The stream adress is given by the user                                         */
/*    -> The function doesn't return anything, no errors are handled                    */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
void emptybuffer ( FILE *stream )
{
   char cc = 0 ;
   while (cc != '\n' && cc != EOF)
   {
      cc = fgetc(stream) ;
   }

   return ;
}









