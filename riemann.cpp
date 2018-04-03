#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace std;
double sgn(double );


/*
 
 Solve riemann shock tube problem for a general equation of state using
 the method of Colella and Glaz.  Use a two shock approximation, and
 linearly interpolation between the head and tail of a rarefaction to
 treat rarefactions.
 
 Take as input the effective left and right states, obtained by
 integrating the parabolic reconstructions of the data over the domain
 of dependence of each of the characteristics on the respective side
 of the interface.  This is accomplished by states().  Return the
 solution to Riemann's problem on axis -- this is used in computing
 the fluxes.
 
 The Riemann problem for the Euler's equation produces 4 states,
 separated by the three characteristics (u - cs, u, u + cs):
 
 
 l_1      t    l_2       l_3
 \       ^   .       /
 \  *L  |   . *R   /
 \     |  .     /
 \    |  .    /
 L    \   | .   /    R
 \  | .  /
 \ |. /
 \|./
 ----------+----------------> x
 
 l_1 = u - cs   eigenvalue
 l_2 = u        eigenvalue (contact)
 l_3 = u + cs   eigenvalue
 
 only density jumps across l_2
 
 References:
 
 CG:   Colella & Glaz 1985, JCP, 59, 264.
 
 CW:   Colella & Woodward 1984, JCP, 54, 174.
 
 Fry:  Fryxell et al. 2000, ApJS, 131, 273.
 
 Toro: Toro 1999, ``Riemann Solvers and Numerical Methods for Fluid
 Dynamcs: A Practical Introduction, 2nd Ed.'', Springer-Verlag
 
 */


int Ugly_riemann_(double *gamma_pass,
             double *rho_l_pass, double *u_l_pass, double *p_l_pass,
             double *rho_r_pass, double *u_r_pass, double *p_r_pass,
             double *F_dens_pass, double *F_xmom_pass, double *F_ener_pass) {
    
    
    double gamma, rho_l, u_l, p_l, rho_r, u_r, p_r;
    
    double riemann_tol;
    int nriem;
    
    double scratch, scratch2;
    double ustar_sgn;
    
    double tau_l, c_l;
    double tau_r, c_r;
    
    double pstar, pstar1, pstar2;
    double w_l, w_l1;
    double w_r, w_r1;
    
    double wes, westar;
    double rhos, rhostr, us, ustar, ps;
    
    double smlrho, smallp, smallu;
    
    double ustar_l1, ustar_r1, ustar_l2, ustar_r2, delu1, delu2;
    
    double rhoav, uav, pav;
    
    double pres_err;
    int has_converged;
    
    double vs, vstar, ces, cestar, ws;
    
    int n;
    
    
    /* some parameters */
    riemann_tol = 1.e-5;
    nriem = 15;
    
    gamma = *gamma_pass;
    
    rho_l = *rho_l_pass;
    u_l = *u_l_pass;
    p_l = *p_l_pass;
    
    rho_r = *rho_r_pass;
    u_r = *u_r_pass;
    p_r = *p_r_pass;
    
    smlrho = 1.e-10;
    smallp = 1.e-10;
    smallu = 1.e-10;
    
    // specific volume
    tau_l = 1.0/fmax(rho_l, smlrho);
    tau_r = 1.0/fmax(rho_r, smlrho);
    
    p_l = fmax(smallp, p_l);
    p_r = fmax(smallp, p_r);
    
    c_l = sqrt(gamma*p_l*rho_l);
    c_r = sqrt(gamma*p_r*rho_r);
    
    
    /*
     construct first guess for secant iteration by assuming that the
     nonlinear wave speed is equal to the sound speed -- the resulting
     expression is the same as Toro, Eq. 9.28 in the Primitive Variable
     Riemann Solver (PVRS).  See also Fry Eq. 72.
     */
    
    pstar1 = p_r - p_l - c_r*(u_r - u_l);
    pstar1 = p_l + pstar1*(c_l/(c_l + c_r));
    pstar1 = fmax(smallp, pstar1);
    
    
    
    /*
     calculate nonlinear wave speeds for the left and right moving
     waves based on the first guess for the pressure jump.  Again,
     there is a left and a right wave speed.  Compute this using CG
     Eq. 34.
     */
    
    /* note -- we simplify a lot here, assuming constant gamma */
    w_l1 = pstar1 + 0.5*(gamma - 1.0)*(pstar1 + p_l);
    w_l1 = sqrt(rho_l*fabs(w_l1));
    
    w_r1 = pstar1 + 0.5*(gamma - 1.0)*(pstar1 + p_r);
    w_r1 = sqrt(rho_r*fabs(w_r1));
    
    
    /*
     construct second guess for the pressure using the nonlinear wave
     speeds from the first guess.  This is basically the same thing we
     did to get pstar1, except now we are using the better wave speeds
     instead of the sound speed.
     */
    
    pstar2 = p_r - p_l - w_r1*(u_r - u_l);
    pstar2 = p_l + pstar2*w_l1/(w_l1 + w_r1);
    pstar2 = fmax(smallp, pstar2);
    
    
    /*
     begin the secant iteration -- see CG Eq. 17 for details.  We will
     continue to interate for convergence until the error falls below
     tol (in which case, things are good), or we hit nriem iterations
     (in which case we have a problem, and we spit out an error).
     */
    
    has_converged = 0;
    
    for (n = 0; n < nriem; n++) {
        
        /* new nonlinear wave speeds, using CG Eq. 34 */
        w_l = pstar2 + 0.5*(gamma  - 1.0)*(pstar2 + p_l);
        w_l = sqrt(rho_l*fabs(w_l));
        
        w_r = pstar2 + 0.5*(gamma - 1.0)*(pstar2 + p_r);
        w_r = sqrt(rho_r*fabs(w_r));
        
        
        /*
         compute the velocities in the "star" state -- using CG
         Eq. 18 -- ustar_l2 and ustar_r2 are the velocities they define
         there.  ustar_l1 and ustar_l2 seem to be the velocities at the
         last time, since pstar1 is the old 'star' pressure, and
         w_l1 is the old wave speed.
         */
        ustar_l1 = u_l - (pstar1 - p_l)/w_l1;
        ustar_r1 = u_r + (pstar1 - p_r)/w_r1;
        ustar_l2 = u_l - (pstar2 - p_l)/w_l;
        ustar_r2 = u_r + (pstar2 - p_r)/w_r;
        
        delu1 = ustar_l1 - ustar_r1;
        delu2 = ustar_l2 - ustar_r2;
        
        scratch = delu2  - delu1;
        
        if (fabs(pstar2 - pstar1) <= smallp) scratch = 0.0;
        
        if (fabs(scratch) < smallu) {
            delu2 = 0.0;
            scratch = 1.0;
        }
        
        /* pressure at the "star" state -- using CG Eq. 18 */
        
        pstar = pstar2 - delu2*(pstar2 - pstar1)/scratch;
        pstar = fmax(smallp, pstar);
        
        /* check for convergence of iteration, riemann_tol is a
         run-time parameter */
        
        pres_err = fabs(pstar - pstar2)/pstar;
        if (pres_err < riemann_tol) {
            has_converged = 1;
            break;
        }
        
        /* reset variables for next iteration */
        pstar1 = pstar2;
        pstar2 = pstar;
        
        w_l1 = w_l;
        w_r1 = w_r;
        
    }
    
    if (!has_converged) {
        
        printf("\n");
        printf("Nonconvergence in subroutine rieman!\n\n");
        printf("Pressure error   = %f\n\n", pres_err);
        printf("pL       = %f,  pR        = %f\n", p_l, p_r);
        printf("uL       = %f,  uR        = %f\n", u_l, u_r);
        printf("cL       = %f,  cR        = %f\n", c_l, c_r);
        printf("\n");
        printf("Terminating execution.\n");
        fflush(NULL);
        exit(-1);
        
    }
    
    
    /* end of secant iteration */
    
    /*
     calculate fluid velocity for the "star" state -- this comes from
     the shock jump equations, Fry Eq. 68 and 69.  The ustar velocity
     can be computed using either the jump eq. for a left moving or
     right moving shock -- we use the average of the two.
     */
    
    scratch = u_l - (pstar - p_l)/w_l;
    scratch2 = u_r + (pstar - p_r)/w_r;
    ustar = 0.5*(scratch + scratch2);
    
    ustar_sgn = sgn(ustar);
    
    /*
     decide which state is located at the zone iterface based on
     the values of the wave speeds.  This is just saying that if
     ustar > 0, then the state is U_L.  if ustar < 0, then the
     state on the axis is U_R.
     */
    
    scratch = 0.50*(1.0 + ustar_sgn);
    scratch2 = 0.5e0*(1.0 - ustar_sgn);
    
    ps = p_l*scratch + p_r*scratch2;
    us = u_l*scratch + u_r*scratch2;
    vs = tau_l*scratch + tau_r*scratch2;
    
    rhos = 1.0/vs;
    rhos = fmax(smlrho, rhos);
    
    vs = 1.0/rhos;
    ws = w_l*scratch + w_r*scratch2;
    ces = sqrt(gamma*ps*vs);
    
    /* compute rhostar, using the shock jump condition (Fry Eq. 80) */
    vstar = vs - (pstar - ps)/(ws*ws);
    rhostr = 1.0/ vstar;
    cestar = sqrt(gamma*pstar*vstar);
    
    /* compute some factors, Fry Eq. 81 and 82 */
    wes = ces - ustar_sgn*us;
    westar = cestar - ustar_sgn*ustar;
    
    scratch = ws*vs - ustar_sgn*us;
    
    if (pstar - ps >= 0.0) {
        wes = scratch;
        westar = scratch;
    }
    
    /* compute correct state for rarefaction fan by linear interpolation */
    scratch = fmax(wes - westar, wes + westar);
    scratch = fmax(scratch, smallu);
    
    scratch = (wes + westar)/scratch;
    
    scratch = 0.5*(1.0 + scratch);
    scratch2 = 1.0 - scratch;
    
    rhoav = scratch*rhostr + scratch2*rhos;
    uav = scratch*ustar + scratch2*us;
    pav = scratch*pstar + scratch2*ps;
    
    if (westar >= 0.0) {
        rhoav  = rhostr;
        uav = ustar;
        pav = pstar;
    }
    
    if (wes < 0.0) {
        rhoav = rhos;
        uav = us;
        pav = ps;
    }
    
    
    /* now compute the fluxes */
    *F_dens_pass = rhoav*uav;
    *F_xmom_pass = rhoav*uav*uav + pav;
    *F_ener_pass = uav*(pav/(gamma - 1.0) + 0.5*rhoav*uav*uav + pav);
    
    
    return 1;
    
}

inline double sgn(double x) {
    if (x < 0.0) {
        return -1.0;
    } else if (x == 0.0) {
        return 0.0;
    } else
        return 1.0;
}
