#ifndef _QHG_SU3_OPS_H
#define _QHG_SU3_OPS_H 1

#include <qhg_defs.h>
#include <qhg_types.h>
#include <qhg_su3_linalg.h>
#include <qhg_su3_mul.h>

#define GVMC(v, m, a, b) (b + NC*(a + NC*(m + ND*v)))

static void
su3_load(_Complex double A[NC][NC], qhg_gauge_field g, int v, int mu)
{
  for(int c0=0; c0<NC; c0++)
    for(int c1=0; c1<NC; c1++) {
      A[c0][c1] = g.field[GVMC(v, mu, c0, c1)];
    }
  return;
}

#endif /* _QHG_SU3_OPS_H */
