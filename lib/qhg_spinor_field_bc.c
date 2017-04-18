#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <mpi.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

#include <qhg_types.h>
#include <qhg_defs.h>
#include <qhg_spinor_linalg.h>
#include <qhg_spinor_field.h>

void
qhg_spinor_twist_t_bc(qhg_spinor_field out, qhg_spinor_field in, double angle)
{
  qhg_lattice *lat = in.lat;
  int Lt = lat->dims[0];
  int lt = lat->ldims[0];
  int t0 = lat->comms->proc_coords[0]*lt;  
  unsigned long int lv3 = lat->lv3;
  for(int t=0; t<lt; t++) {
    int gt = t0+t;
    double phi = (M_PI*angle*gt)/(double)Lt;
    _Complex double phase = cos(phi) + _Complex_I*sin(phi);
    for(unsigned long int v=lv3*t; v<lv3*(t+1); v++) {
      _Complex double *p = &(out.field[v*NC*NS]);
      _Complex double *q = &(in.field[v*NC*NS]);
      memcpy(p, q, NC*NS*sizeof(_Complex double));
      spinor_linalg_ax(phase, p);
    }
  }
  return;
}

void
qhg_spinor_set_bc(qhg_spinor_field psi, enum qhg_fermion_bc_time bc)
{  
  psi.bc = bc;
  return;
}
