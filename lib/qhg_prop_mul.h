#ifndef _QHG_PROP_MUL_H
#define _QHG_PROP_MUL_H 1
#include <complex.h>
#include <qhg_idx.h>

static void
prop_mul_gg(_Complex double Z[NS*NC][NS*NC], _Complex double X[NS*NC][NS*NC], _Complex double Y[NS*NC][NS*NC])
{
  for(int cs0=0; cs0<NC*NS; cs0++)
    for(int cs1=0; cs1<NC*NS; cs1++) {
      Z[cs0][cs1] = 0;
      for(int cs2=0; cs2<NC*NS; cs2++)
	Z[cs0][cs1] += X[cs0][cs2]*Y[cs2][cs1];
    }
  return;
}

static void
prop_mul_gd(_Complex double Z[NS*NC][NS*NC], _Complex double X[NS*NC][NS*NC], _Complex double Y[NS*NC][NS*NC])
{
  for(int cs0=0; cs0<NC*NS; cs0++)
    for(int cs1=0; cs1<NC*NS; cs1++) {
      Z[cs0][cs1] = 0;
      for(int cs2=0; cs2<NC*NS; cs2++)
	Z[cs0][cs1] += X[cs0][cs2]*conj(Y[cs1][cs2]);
    }
  return;
}

static void
prop_mul_dg(_Complex double Z[NS*NC][NS*NC], _Complex double X[NS*NC][NS*NC], _Complex double Y[NS*NC][NS*NC])
{
  for(int cs0=0; cs0<NC*NS; cs0++)
    for(int cs1=0; cs1<NC*NS; cs1++) {
      Z[cs0][cs1] = 0;
      for(int cs2=0; cs2<NC*NS; cs2++)
	Z[cs0][cs1] += conj(X[cs2][cs0])*Y[cs2][cs1];
    }
  return;
}

static void
prop_mul_dd(_Complex double Z[NS*NC][NS*NC], _Complex double X[NS*NC][NS*NC], _Complex double Y[NS*NC][NS*NC])
{
  for(int cs0=0; cs0<NC*NS; cs0++)
    for(int cs1=0; cs1<NC*NS; cs1++) {
      Z[cs0][cs1] = 0;
      for(int cs2=0; cs2<NC*NS; cs2++)
	Z[cs0][cs1] += conj(X[cs2][cs0])*conj(Y[cs1][cs2]);
    }
  return;
}

/*
 * This is not your standard matrix-matrix multiplication. Returns:
 * Z^{a b}_{mu nu} = sum_{k, c} X^{ac}_{mu nu} Y^{cb}_{k k}
 */
static void
prop_mul_gtr(_Complex double Z[NS*NC][NS*NC], _Complex double X[NS*NC][NS*NC], _Complex double Y[NS*NC][NS*NC])
{
  for(int sp0=0; sp0<NS; sp0++)
    for(int sp1=0; sp1<NS; sp1++)
      for(int col0=0; col0<NC; col0++)
	for(int col1=0; col1<NC; col1++) {
	  Z[CS(sp0, col0)][CS(sp1, col1)] = 0.0;
	  for(int sp2=0; sp2<NS; sp2++)
	    for(int col2=0; col2<NC; col2++)
	      Z[CS(sp0, col0)][CS(sp1, col1)] +=
		X[CS(sp0, col0)][CS(sp1, col2)]*Y[CS(sp2, col2)][CS(sp2, col1)];
    }
  return;
}

#endif /* _QHG_PROP_MUL_H */
