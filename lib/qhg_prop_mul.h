#ifndef _QHG_PROP_MUL_H
#define _QHG_PROP_MUL_H 1
#include <complex.h>

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

#endif /* _QHG_PROP_MUL_H */
