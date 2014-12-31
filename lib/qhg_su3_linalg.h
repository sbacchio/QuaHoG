#ifndef _QHG_SU3_LINALG_H
#define _QHG_SU3_LINALG_H 1
#include <complex.h>
#include <qhg_idx.h>
/*
 * Perform basic linear algebra on 3x3 matrices
 *
 */

/*
 * scale in-place an su3 matrix by a complex number
 */
static void
su3_linalg_au(_Complex double a, _Complex double *u)
{
  for(int i=0; i<NC*NC; i++)
    u[i] = a*u[i];

  return;
}

/*
 * u = v, with u, v 3x3
 */
static void
su3_linalg_ueqv(_Complex double *u, _Complex double *v)
{
  for(int i=0; i<NC*NC; i++)
    u[i] = v[i];

  return;
}

/*
 * u = v^+, with u, v 3x3
 */
static void
su3_linalg_ueqvd(_Complex double *u, _Complex double *v)
{
  for(int i=0; i<NC; i++)
    for(int j=0; j<NC; j++)    
      u[CC(i, j)] = conj(v[CC(j, i)]);

  return;
}

/*
 * u = u + v, with u, v 3x3
 */
static void
su3_linalg_upeqv(_Complex double *u, _Complex double *v)
{
  for(int i=0; i<NC*NC; i++)
    u[i] = u[i] + v[i];

  return;
}

/*
 * u_ij = 0 i,j = 1,..,3
 */
static void
su3_linalg_zero(_Complex double *u)
{
  for(int i=0; i<NC*NC; i++)
    u[i] = 0;
  
  return;
}

/*
 * Set u_ii = z_i, i=1,...,3
 */
static void
su3_linalg_diag(_Complex double *u, _Complex double z[NC])
{
  for(int i=0; i<NC*NC; i++)
    u[i] = 0;
  
  for(int i=0; i<NC; i++)
    u[CC(i, i)] = z[i];
  
  return;
}

/*
 * return trace(u), with u 3x3
 */
static _Complex double
su3_linalg_trace_u(_Complex double *u)
{
  _Complex double trace = 0;
  for(int i=0; i<NC; i++)
    trace += u[CC(i, i)];
  
  return trace;
}

/*
 * return det(u), with u 3x3
 */
static _Complex double
su3_linalg_det_u(_Complex double *u)
{
  _Complex double det = 0;
  det +=  u[CC(0,0)]*(u[CC(1,1)]*u[CC(2,2)] - u[CC(1,2)]*u[CC(2,1)]);
  det += -u[CC(0,1)]*(u[CC(1,0)]*u[CC(2,2)] - u[CC(1,2)]*u[CC(2,0)]);
  det +=  u[CC(0,2)]*(u[CC(1,0)]*u[CC(2,1)] - u[CC(1,1)]*u[CC(2,1)]);  
  return det;
}

#endif /* _QHG_SU3_LINALG_H */
