#ifndef _QHG_SPINOR_LINALG_H
#define _QHG_SPINOR_LINALG_H 1
#include <complex.h>
#include <qhg_idx.h>

/*
 * Perform basic linear algebra on spinors and/or spinor * su(3)
 */

/*
 * scale in-place an spinor by a complex number
 */
static void
spinor_linalg_ax(_Complex double a, _Complex double *x)
{
  for(int i=0; i<NC*NS; i++)
    x[i] = a*x[i];

  return;
}

/*
 * y += x, with x, y spinors
 */
static void
spinor_linalg_ypeqx(_Complex double *y, _Complex double *x)
{
  for(int i=0; i<NC*NS; i++)
    y[i] += x[i];

  return;
}

/*
 * x = 0, with x a spinor
 */
static void
spinor_linalg_zero(_Complex double *x)
{
  for(int i=0; i<NC*NS; i++)
    x[i] = 0;
  
  return;
}

/*
 * y += u*x, with y, x spinors and u an su(3)
 */
static void
spinor_linalg_ypequx(_Complex double *y, _Complex double *u, _Complex double *x)
{
  for(int s=0; s<NS; s++)
    for(int c0=0; c0<NC; c0++)
      for(int c1=0; c1<NC; c1++)
	y[CS(s, c0)] += u[CC(c0, c1)]*x[CS(s, c1)];
  
  return;
}

/*
 * y += u^+*x, with y, x spinors and u an su(3)
 */
static void
spinor_linalg_ypeqdx(_Complex double *y, _Complex double *u, _Complex double *x)
{
  for(int s=0; s<NS; s++)
    for(int c0=0; c0<NC; c0++)
      for(int c1=0; c1<NC; c1++)
	y[CS(s, c0)] += conj(u[CC(c1, c0)])*x[CS(s, c1)];
      
  return;
}
#endif /* _QHG_SU3_LINALG_H */
