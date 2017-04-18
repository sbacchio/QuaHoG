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

void
qhg_spinors_ax(_Complex double a, qhg_spinor_field psi[], int n_spinors)
{
  qhg_lattice *lat = psi[0].lat;
  unsigned long int lvol = lat->lvol;
  for(int isp=0; isp<n_spinors; isp++) 
    for(unsigned long int v=0; v<lvol; v++) {
      _Complex double *q = &(psi[isp].field[v*NC*NS]);
      spinor_linalg_ax(a, q);
    }

  return;
}
