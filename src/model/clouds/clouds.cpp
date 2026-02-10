/*! \file clouds.cpp
 *  \brief Definitions of initial conditions for hydrostatic disks.
           Note that the grid is mapped to 1D as i + (x_dim)*j +
 (x_dim*y_dim)*k. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

#include "../../global/global.h"
#include "../../grid/grid3D.h"
#include "../../io/io.h"
#include "../../mpi/mpi_routines.h"
#include "../../utils/error_handling.h"
#include "../../utils/hydro_utilities.h"
#include "../../utils/math_utilities.h"
#include "../../utils/mhd_utilities.h"


/*! \fn void Clouds()
 *  \brief Bunch of clouds. */
void Grid3D::Clouds(Parameters P)
{
  int i, j, k, id;
  int istart, jstart, kstart, iend, jend, kend;
  Real x_pos, y_pos, z_pos;
  Real n_bg, n_cl;      // background and cloud number density
  Real rho_bg, rho_cl;  // background and cloud density
  Real vx_bg, vx_cl;    // background and cloud velocity
  Real vy_bg, vy_cl;
  Real vz_bg, vz_cl;
  Real T_bg, T_cl;       // background and cloud temperature
  Real p_bg, p_cl;       // background and cloud pressure
  Real mu   = 0.6;       // mean atomic weight
  int N_cl  = 1;         // number of clouds
  Real R_cl = 2.5;       // cloud radius in code units (kpc)
  Real cl_pos[N_cl][3];  // array of cloud positions
  Real r;

  // Multiple Cloud Setup
  // for (int nn=0; nn<N_cl; nn++) {
  //  cl_pos[nn][0] = (nn+1)*0.1*H.xdglobal+0.5*H.xdglobal;
  //  cl_pos[nn][1] = (nn%2*0.1+0.45)*H.ydglobal;
  //  cl_pos[nn][2] = 0.5*H.zdglobal;
  //  printf("Cloud positions: %f %f %f\n", cl_pos[nn][0], cl_pos[nn][1],
  //  cl_pos[nn][2]);
  //}

  // single centered cloud setup
  for (int nn = 0; nn < N_cl; nn++) {
    cl_pos[nn][0] = 0.5 * H.xdglobal;
    cl_pos[nn][1] = 0.5 * H.ydglobal;
    cl_pos[nn][2] = 0.5 * H.zdglobal;
    printf("Cloud positions: %f %f %f\n", cl_pos[nn][0], cl_pos[nn][1], cl_pos[nn][2]);
  }

  n_bg   = 1.68e-4;
  n_cl   = 5.4e-2;
  rho_bg = n_bg * mu * MP / DENSITY_UNIT;
  rho_cl = n_cl * mu * MP / DENSITY_UNIT;
  vx_bg  = 0.0;
  // vx_c  = -200*TIME_UNIT/KPC; // convert from km/s to kpc/kyr
  vx_cl = 0.0;
  vy_bg = vy_cl = 0.0;
  vz_bg = vz_cl = 0.0;
  T_bg          = 3e6;
  T_cl          = 1e4;
  p_bg          = n_bg * KB * T_bg / PRESSURE_UNIT;
  p_cl          = p_bg;

  istart = H.n_ghost;
  iend   = H.nx - H.n_ghost;
  if (H.ny > 1) {
    jstart = H.n_ghost;
    jend   = H.ny - H.n_ghost;
  } else {
    jstart = 0;
    jend   = H.ny;
  }
  if (H.nz > 1) {
    kstart = H.n_ghost;
    kend   = H.nz - H.n_ghost;
  } else {
    kstart = 0;
    kend   = H.nz;
  }

  // set initial values of conserved variables
  for (k = kstart; k < kend; k++) {
    for (j = jstart; j < jend; j++) {
      for (i = istart; i < iend; i++) {
        // get cell index
        id = i + j * H.nx + k * H.nx * H.ny;

        // get cell-centered position
        Get_Position(i, j, k, &x_pos, &y_pos, &z_pos);

        // set background state
        C.density[id]    = rho_bg;
        C.momentum_x[id] = rho_bg * vx_bg;
        C.momentum_y[id] = rho_bg * vy_bg;
        C.momentum_z[id] = rho_bg * vz_bg;
        C.Energy[id]     = p_bg / (gama - 1.0) + 0.5 * rho_bg * (vx_bg * vx_bg + vy_bg * vy_bg + vz_bg * vz_bg);
#ifdef DE
        C.GasEnergy[id] = p_bg / (gama - 1.0);
#endif
#ifdef SCALAR
  #ifdef DUST
        C.host[id + H.n_cells * grid_enum::dust_density] = 0.0;
  #endif
#endif
        // add clouds
        for (int nn = 0; nn < N_cl; nn++) {
          r = sqrt((x_pos - cl_pos[nn][0]) * (x_pos - cl_pos[nn][0]) +
                   (y_pos - cl_pos[nn][1]) * (y_pos - cl_pos[nn][1]) +
                   (z_pos - cl_pos[nn][2]) * (z_pos - cl_pos[nn][2]));
          if (r < R_cl) {
            C.density[id]    = rho_cl;
            C.momentum_x[id] = rho_cl * vx_cl;
            C.momentum_y[id] = rho_cl * vy_cl;
            C.momentum_z[id] = rho_cl * vz_cl;
            C.Energy[id]     = p_cl / (gama - 1.0) + 0.5 * rho_cl * (vx_cl * vx_cl + vy_cl * vy_cl + vz_cl * vz_cl);
#ifdef DE
            C.GasEnergy[id] = p_cl / (gama - 1.0);
#endif  // DE
#ifdef SCALAR
  #ifdef DUST
            C.host[id + H.n_cells * grid_enum::dust_density] = rho_cl * 1e-2;
  #endif  // DUST
#endif    // SCALAR
          }
        }
      }
    }
  }
}
