/*! \file exact.cpp
 *  \brief Function definitions for Toro exact Riemann solver. */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"global.h"
#include"exact.h"


//#define LOUD


/* \fn Calculate_Exact_Fluxes(Real cW[], Real fluxes[], Real gamma)
 * \brief Returns the density, momentum, and Energy fluxes at an interface.
   Inputs are an array containg left and right density, momentum, and Energy. */
void Calculate_Exact_Fluxes(Real cW[], Real fluxes[], Real gamma)
{
  Real dl, vxl, vyl, vzl, pl, cl; //density, velocity, pressure, sound speed (left)
  Real dr, vxr, vyr, vzr, pr, cr; //density, velocity, pressure, sound speed (right)
  Real ds, us, ps, Es; //sampled density, velocity, pressure, total energy
  Real um, pm; //velocity and pressure in the star region
  #ifdef DE
  Real gel, ger;
  #endif

  // calculate primative variables from input array
  dl = cW[0];
  dr = cW[1];
  vxl = cW[2] / dl;
  vxr = cW[3] / dr;
  vyl = cW[4] / dl;
  vyr = cW[5] / dr;
  vzl = cW[6] / dl;
  vzr = cW[7] / dr;
  pl = (cW[8] - 0.5*dl*(vxl*vxl + vyl*vyl + vzl*vzl)) * (gamma-1.0);
  pl = fmax(pl, TINY_NUMBER);
  pr = (cW[9] - 0.5*dr*(vxr*vxr + vyr*vyr + vzr*vzr)) * (gamma-1.0);
  pr = fmax(pr, TINY_NUMBER);
  #ifdef DE
  gel = cW[10] / dl;
  ger = cW[11] / dr;
  #endif


  //compute sound speeds in left (cell i-1) and right (cell i) regions
  cl = sqrt(gamma*pl/dl);
  cr = sqrt(gamma*pr/dr);


  //test for the pressure positivity condition
  /*
  if ((2.0 / (gamma - 1.0))*(cl+cr) <= (vxr-vxl))
  {
    //the initial data is such that vacuum is generated. Program stopped.
    printf("Vacuum is generated by initial data.\n");
    printf("%f %f %f %f %f %f\n", dl, vxl, pl, dr, vxr, pr);
    //exit(1);
  }
  */
 
  //find exact solution for pressvxre and velocity in star region
  starpu(&pm, &um, dl, vxl, pl, cl, dr, vxr, pr, cr, gamma);
 
  //sample the solution at the cell interface
  sample(pm, um, &ds, &us, &ps, dl, vxl, pl, cl, dr, vxr, pr, cr, gamma);
 
  // calculate the fluxes through the cell interface
  fluxes[0] = ds*us;
  fluxes[1] = ds*us*us+ps;
  if (us >= 0)
  {
    fluxes[2] = ds*us*vyl;
    fluxes[3] = ds*us*vzl;
    #ifdef DE
    fluxes[5] = ds*us*gel;
    #endif
    Es = (ps/(gamma - 1.0)) + 0.5*ds*(us*us + vyl*vyl + vzl*vzl);
  }
  else
  {
    fluxes[2] = ds*us*vyr;
    fluxes[3] = ds*us*vzr;
    #ifdef DE
    fluxes[5] = ds*us*ger;
    #endif
    Es = (ps/(gamma - 1.0)) + 0.5*ds*(us*us + vyr*vyr + vzr*vzr);
  }
  fluxes[4] = (Es+ps)*us;

}


Real guessp(Real dl, Real vxl, Real pl, Real cl, Real dr, Real vxr, Real pr, Real cr, Real gamma)
{
  // purpose:  to provide a guessed value for pressure
  //    pm in the Star Region. The choice is made
  //    according to adaptive Riemann solver using
  //    the PVRS, TRRS and TSRS approximate Riemann
  //    solvers. See Sect. 9.5 of Toro (1999)

  Real gl, gr, ppv, p_0;
  const Real TOL = 1.0e-6;
 
  // compute guess pressure from PVRS Riemann solver
  ppv = 0.5*(pl + pr) + 0.125*(vxl - vxr)*(dl + dr)*(cl + cr);
 
  if (ppv < 0.0) ppv = 0.0;
  // Two-Shock Riemann solver with PVRS as estimate
  gl = sqrt( (2.0/((gamma + 1.0)*dl)) / (((gamma - 1.0) / (gamma + 1.0))*pl + ppv) );
  gr = sqrt((2.0 / ((gamma + 1.0)*dr))/(((gamma - 1.0) / (gamma + 1.0))*pr + ppv));  
  p_0 = (gl*pl + gr*pr - (vxr-vxl))/(gl + gr);

  if (p_0 < 0.0) p_0 = TOL;

  return p_0;

}
 
 
void prefun(Real *f, Real *fd, Real p, Real *dk, Real *pk, Real *ck, Real gamma)
{
  // purpose:  to evaluate the pressure functions
  //    fl and fr in exact Riemann solver
  //    and their first derivatives
 
  Real ak, bk, pratio, qrt;
 
  if (p <= *pk) {
    // rarefaction wave
    pratio = p / *pk;
    *f = (2.0 / (gamma - 1.0)) * *ck * (pow(pratio, ((gamma - 1.0)/(2.0 * gamma))) - 1.0);
    *fd = (1.0/(*dk * *ck))*pow(pratio, -((gamma + 1.0)/(2.0 * gamma)));
  } 
  else 
  {
    // shock wave
    ak = (2.0 / (gamma + 1.0)) / *dk;
    bk = ((gamma - 1.0) / (gamma + 1.0)) * *pk;
    qrt = sqrt(ak/(bk + p));
    *f = (p - *pk)*qrt;
    *fd = (1.0 - 0.5*(p - *pk)/(bk + p))*qrt;
  }
}
 
 
void starpu(Real *p, Real *u, Real dl, Real vxl, Real pl, Real cl, Real dr, Real vxr, Real pr, Real cr, Real gamma)
{
  // purpose:  to compute the solution for pressure and
  //   velocity in the Star Region
 
  const int nriter = 20;
  const Real TOL = 1.0e-6;
  Real change, fl, fld, fr, frd, pold, pstart;
 
  //guessed value pstart is computed
  pstart = guessp(dl, vxl, pl, cl, dr, vxr, pr, cr, gamma);
  pold = pstart;

  int i = 1;
  for (i=0 ; i <= nriter; i++) {
    prefun(&fl, &fld, pold, &dl, &pl, &cl, gamma);
    prefun(&fr, &frd, pold, &dr, &pr, &cr, gamma);
    *p = pold - (fl + fr + vxr - vxl)/(fld + frd);
    change = 2.0*fabs((*p - pold)/(*p + pold));


    if (change <= TOL) break;
    if (*p < 0.0) *p = TOL;
    pold = *p;
  }

  if (i > nriter) {
    //printf("Divergence in Newton-Raphson iteration. p = %e\n", *p);
    //printf("%f %f %f %f %f %f\n", dl, vxl, pl, dr, vxr, pr);
    //exit(0);
  }

 
  // compute velocity in star region
  *u = 0.5*(vxl + vxr + fr - fl);

}
 
 
void sample(const Real pm, const Real vm, 
      Real *d, Real *u, Real *p,
      Real dl, Real vxl, Real pl, Real cl,
      Real dr, Real vxr, Real pr, Real cr, Real gamma)
{
  // purpose:  to sample the solution throughout the wave
  //   pattern. Pressure pm and velocity vm in the
  //   star region are known. Sampled
  //   values are d, u, p.
 
  Real c, cml, cmr, pml, pmr, sl, sr;
 
  if (vm >= 0) // sampling point lies to the left of the contact discontinuity
  {
    if (pm <= pl) // left rarefaction
    {    
      if (vxl - cl >= 0) // sampled point is in left data state
      {   
        *d = dl;
        *u = vxl;
        *p = pl;
      }
      else 
      {
        cml = cl*pow(pm/pl, ((gamma - 1.0)/(2.0 * gamma)));
        if (vm - cml < 0) // sampled point is in star left state
        {
          *d = dl*pow(pm/pl, 1.0/gamma);
          *u = vm;
          *p = pm;
        } 
        else // sampled point is inside left fan
        {
          *u = (2.0 / (gamma + 1.0))*(cl + ((gamma - 1.0) / 2.0)*vxl);
          c = (2.0 / (gamma + 1.0))*(cl + ((gamma - 1.0) / 2.0)*vxl);
          *d = dl*pow(c/cl, (2.0 / (gamma - 1.0)));
          *p = pl*pow(c/cl, (2.0 * gamma / (gamma - 1.0)));
        }
      }
    } 
    else // left shock
    { 
      pml = pm/pl;
      sl = vxl - cl*sqrt(((gamma + 1.0)/(2.0 * gamma))*pml + ((gamma - 1.0)/(2.0 * gamma)));
      if (sl >= 0) // sampled point is in left data state
      {
        *d = dl;
        *u = vxl;
        *p = pl;
      } 
      else // sampled point is in star left state
      { 
        *d = dl*(pml + ((gamma - 1.0) / (gamma + 1.0)))/(pml*((gamma - 1.0) / (gamma + 1.0)) + 1.0);
        *u = vm;
        *p = pm;
      }
    } 
  } 
  else // sampling point lies to the right of the contact discontinuity
  { 
    if (pm > pr) // right shock
    {
      pmr = pm/pr;
      sr = vxr + cr*sqrt(((gamma + 1.0)/(2.0 * gamma))*pmr + ((gamma - 1.0)/(2.0 * gamma)));
      if (sr <= 0) // sampled point is in right data state
      {
        *d = dr;
        *u = vxr;
        *p = pr;
      } 
      else // sampled point is in star right state
      { 
        *d = dr*(pmr + ((gamma - 1.0) / (gamma + 1.0)))/(pmr*((gamma - 1.0) / (gamma + 1.0)) + 1.0);
        *u = vm;
        *p = pm;
      }
    } 
    else // right rarefaction
    { 
      if (vxr + cr <= 0) // sampled point is in right data state
      { 
        *d = dr;
        *u = vxr;
        *p = pr;
      } 
      else 
      {
        cmr = cr*pow(pm/pr, ((gamma - 1.0)/(2.0 * gamma)));
        if (vm + cmr >= 0) // sampled point is in star right state
        {    
          *d = dr*pow(pm/pr, 1.0/gamma);
          *u = vm;
          *p = pm;
        } 
        else // sampled point is inside right fan
        {    
          *u = (2.0 / (gamma + 1.0))*(-cr + ((gamma - 1.0) / 2.0)*vxr);
          c = (2.0 / (gamma + 1.0))*(cr - ((gamma - 1.0) / 2.0)*vxr);
          *d = dr*pow(c/cr, (2.0 / (gamma - 1.0)));
          *p = pr*pow(c/cr, (2.0 * gamma / (gamma - 1.0)));
        }
      }
    } 
  }
}


