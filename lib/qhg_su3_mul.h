#ifndef _QHG_SU3_MUL_H
#define _QHG_SU3_MUL_H 1
#include <complex.h>
#include <qhg_idx.h>
/*
 * Perform dense matrix mul 3x3 times 3x3
 *
 */

/*
 * w = u v
 */
static inline void
su3_mul_uu(_Complex double *w, _Complex double *u, _Complex double *v)
{
  for(int i=0; i<NC; i++)
    for(int j=0; j<NC; j++) {
      w[CC(i,j)] = 0;
      for(int k=0; k<NC; k++) {
	w[CC(i,j)] += u[CC(i,k)]*v[CC(k,j)];
      }
    }
  return;
}

/*
 * w = u v^\dagger
 */
static inline void
su3_mul_ud(_Complex double *w, _Complex double *u, _Complex double *v)
{
  for(int i=0; i<NC; i++)
    for(int j=0; j<NC; j++) {
      w[CC(i,j)] = 0;
      for(int k=0; k<NC; k++) {
	w[CC(i,j)] += u[CC(i,k)]*conj(v[CC(j,k)]);
      }
    }
  return;
}

/*
 * w = u^\dagger v
 */
static inline void
su3_mul_du(_Complex double *w, _Complex double *u, _Complex double *v)
{
  for(int i=0; i<NC; i++)
    for(int j=0; j<NC; j++) {
      w[CC(i,j)] = 0;
      for(int k=0; k<NC; k++) {
	w[CC(i,j)] += conj(u[CC(k,i)])*v[CC(k,j)];
      }
    }
  return;
}

/*
 * w = u^\dagger v^\dagger
 */
static inline void
su3_mul_dd(_Complex double *w, _Complex double *u, _Complex double *v)
{
  for(int i=0; i<NC; i++)
    for(int j=0; j<NC; j++) {
      w[CC(i,j)] = 0;
      for(int k=0; k<NC; k++) {
	w[CC(i,j)] += conj(u[CC(k,i)])*conj(v[CC(j,k)]);
      }
    }
  return;
}
#endif /* _QHG_SU3_MUL_H */
