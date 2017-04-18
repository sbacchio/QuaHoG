#ifndef _QHG_PROP_OPS_H
#define _QHG_PROP_OPS_H 1

#define VSC(v, sc) (sc + NC*NS*v)

#include <qhg_prop_mul.h>
#include <qhg_prop_mul_su3.h>
#include <qhg_prop_contract.h>
#include <qhg_prop_trace.h>

static void
prop_load(_Complex double A[NC*NS][NC*NS], qhg_spinor_field s[NC*NS], unsigned long int v)
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
prop_store(qhg_spinor_field s[NC*NS], unsigned long int v, _Complex double A[NC*NS][NC*NS])
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
prop_eye(_Complex double A[NC*NS][NC*NS])
{
  prop_zero(A);
  for(int cs0=0; cs0<NS*NC; cs0++) {
    A[cs0][cs0] = 1.;
  }
  return;
}

#endif /* _QHG_PROP_OPS_H */
