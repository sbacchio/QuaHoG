#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <mpi.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

#include <qhg_types.h>
#include <qhg_defs.h>
#include <qhg_spinor_linalg.h>
#include <qhg_spinor_field.h>

void
qhg_spinors_untwist_bc(qhg_spinor_field psi[], int n_spinors, double theta, int t_origin)
{
  qhg_lattice *lat = psi[0].lat;
  int Lt = lat->dims[0];
  int lt = lat->ldims[0];
  int t0 = lat->comms->proc_coords[0]*lt;  
  int lv3 = lat->lv3;
  for(int isp=0; isp<n_spinors; isp++) {
    for(int t=0; t<lt; t++) {
      int gt = t0+t;
      int dt = gt - t_origin;
      double phi = (M_PI*theta*dt)/Lt;
      _Complex double phase = cos(phi) + _Complex_I*sin(phi);
      for(int v=lv3*t; v<lv3*(t+1); v++) {
	_Complex double *q = &(psi[isp].field[v*NC*NS]);
	spinor_linalg_ax(phase, q);
      }
    }
  }
  return;
}

void
qhg_spinors_set_bc(qhg_spinor_field psi[], int n_spinors, _Complex double bc[ND])
{
  for(int isp=0; isp<n_spinors; isp++)
    qhg_spinor_field_set_bc(psi[isp], bc);
  
  return;
}
