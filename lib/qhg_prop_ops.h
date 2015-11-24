#ifndef _QHG_PROP_OPS_H
#define _QHG_PROP_OPS_H 1

#define VSC(v, sc) (sc + NC*NS*v)
#include <qhg_idx.h>
#include <qhg_prop_mul.h>
#include <qhg_prop_mul_su3.h>
#include <qhg_prop_contract.h>
#include <qhg_prop_trace.h>

static void
prop_load(_Complex double A[NC*NS][NC*NS], qhg_spinor_field s[NC*NS], int v)
{
  for(int cs0=0; cs0<NS*NC; cs0++)
    for(int cs1=0; cs1<NS*NC; cs1++) {
      A[cs0][cs1] = s[cs1].field[VSC(v, cs0)];
    }
  return;
}

static void
prop_add_gg(_Complex double A[NC*NS][NC*NS], _Complex double B[NC*NS][NC*NS], _Complex double C[NC*NS][NC*NS])
{
  for(int cs0=0; cs0<NS*NC; cs0++)
    for(int cs1=0; cs1<NS*NC; cs1++) {
      A[cs0][cs1] = B[cs0][cs1] + C[cs0][cs1];
    }
  return;
}

static void
prop_sub_gg(_Complex double A[NC*NS][NC*NS], _Complex double B[NC*NS][NC*NS], _Complex double C[NC*NS][NC*NS])
{
  for(int cs0=0; cs0<NS*NC; cs0++)
    for(int cs1=0; cs1<NS*NC; cs1++) {
      A[cs0][cs1] = B[cs0][cs1] - C[cs0][cs1];
    }
  return;
}

static void
prop_peq_g(_Complex double A[NC*NS][NC*NS], _Complex double B[NC*NS][NC*NS])
{
  for(int cs0=0; cs0<NS*NC; cs0++)
    for(int cs1=0; cs1<NS*NC; cs1++) {
      A[cs0][cs1] += B[cs0][cs1];
    }
  return;
}

static void
prop_meq_g(_Complex double A[NC*NS][NC*NS], _Complex double B[NC*NS][NC*NS])
{
  for(int cs0=0; cs0<NS*NC; cs0++)
    for(int cs1=0; cs1<NS*NC; cs1++) {
      A[cs0][cs1] -= B[cs0][cs1];
    }
  return;
}

static void
prop_store(qhg_spinor_field s[NC*NS], int v, _Complex double A[NC*NS][NC*NS])
{
  for(int cs0=0; cs0<NS*NC; cs0++)
    for(int cs1=0; cs1<NS*NC; cs1++) {
      s[cs1].field[VSC(v, cs0)] = A[cs0][cs1];
    }
  return;
}

static void
prop_scale(_Complex double a, _Complex double A[NC*NS][NC*NS])
{
  for(int cs0=0; cs0<NS*NC; cs0++)
    for(int cs1=0; cs1<NS*NC; cs1++) {
      A[cs0][cs1] = a*A[cs0][cs1];
    }
  return;
}

static void
prop_zero(_Complex double A[NC*NS][NC*NS])
{
  for(int cs0=0; cs0<NS*NC; cs0++)
    for(int cs1=0; cs1<NS*NC; cs1++) {
      A[cs0][cs1] = 0.0;
    }
  return;
}

static void
prop_transpose(_Complex double A[NC*NS][NC*NS])
{
  for(int csp0=0; csp0<NC*NS-1; csp0++)
    for(int csp1=csp0+1; csp1<NC*NS; csp1++) {
      _Complex double swap = A[csp0][csp1];
      A[csp0][csp1] = A[csp1][csp0];
      A[csp1][csp0] = swap;
    }

  return;
}

static void
prop_color_transpose(_Complex double A[NC*NS][NC*NS])
{
  for(int sp0=0; sp0<NS; sp0++)
    for(int sp1=0; sp1<NS; sp1++)    
      for(int col0=0; col0<NC-1; col0++)
	for(int col1=col0+1; col1<NC; col1++) {
	  _Complex double swap = A[CS(sp0, col0)][CS(sp1, col1)];
	  A[CS(sp0, col0)][CS(sp1, col1)] = A[CS(sp0, col1)][CS(sp1, col0)];
	  A[CS(sp0, col1)][CS(sp1, col0)] = swap;
	}
  return;
}

static void
prop_color_trace(_Complex double T[NS][NS], _Complex double A[NC*NS][NC*NS])
{
  for(int sp0=0; sp0<NS; sp0++)
    for(int sp1=0; sp1<NS; sp1++) {
      T[sp0][sp1] = 0.0;
      for(int col0=0; col0<NC; col0++) {
	int cs0 = col0 + sp0*NC;
	int cs1 = col0 + sp1*NC;
	T[sp0][sp1] += A[cs0][cs1];
      }
    }
  return;
}

static void
prop_spinor_transpose(_Complex double A[NC*NS][NC*NS])
{
  for(int sp0=0; sp0<NS-1; sp0++)
    for(int sp1=sp0+1; sp1<NS; sp1++)    
      for(int col0=0; col0<NC; col0++)
	for(int col1=0; col1<NC; col1++) {
	  int cs0a = col0 + sp0*NC;
	  int cs1a = col1 + sp1*NC;
	  int cs0b = col0 + sp1*NC;
	  int cs1b = col1 + sp0*NC;
	  _Complex double swap = A[cs0a][cs1a];
	  A[cs0a][cs1a] =  A[cs0b][cs1b];
	  A[cs1b][cs1b] = swap;
	}
  return;
}

static void
prop_eye(_Complex double A[NC*NS][NC*NS])
{
  prop_zero(A);
  for(int cs0=0; cs0<NS*NC; cs0++) {
    A[cs0][cs0] = 1.;
  }
  return;
}

#endif /* _QHG_PROP_OPS_H */
