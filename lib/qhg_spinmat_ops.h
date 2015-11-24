#ifndef _QHG_SPINMAT_OPS_H
#define _QHG_SPINMAT_OPS_H 1

static void
spinmat_peq_s(_Complex double A[NS][NS], _Complex double B[NS][NS])
{
  for(int sp0=0; sp0<NS; sp0++)
    for(int sp1=0; sp1<NS; sp1++) {
      A[sp0][sp1] += B[sp0][sp1];
    }
  return;
}

static void
spinmat_meq_s(_Complex double A[NS][NS], _Complex double B[NS][NS])
{
  for(int sp0=0; sp0<NS; sp0++)
    for(int sp1=0; sp1<NS; sp1++) {
      A[sp0][sp1] -= B[sp0][sp1];
    }
  return;
}

static void
spinmat_add_ss(_Complex double C[NS][NS], _Complex double A[NS][NS], _Complex double B[NS][NS])
{
  for(int sp0=0; sp0<NS; sp0++)
    for(int sp1=0; sp1<NS; sp1++) {
      C[sp0][sp1] = A[sp0][sp1] + B[sp0][sp1];
    }
  return;
}

static void
spinmat_sub_ss(_Complex double C[NS][NS], _Complex double A[NS][NS], _Complex double B[NS][NS])
{
  for(int sp0=0; sp0<NS; sp0++)
    for(int sp1=0; sp1<NS; sp1++) {
      C[sp0][sp1] = A[sp0][sp1] - B[sp0][sp1];
    }
  return;
}

static _Complex double
spinmat_trace(_Complex double A[NS][NS])
{
  _Complex double trace = 0;
  for(int sp0=0; sp0<NS; sp0++)
    for(int sp1=0; sp1<NS; sp1++) {
      trace += A[sp0][sp1];
    }
  return trace;
}

static void
spinmat_zero(_Complex double A[NS][NS])
{
  for(int sp0=0; sp0<NS; sp0++)
    for(int sp1=0; sp1<NS; sp1++) {
      A[sp0][sp1] = 0.0;
    }
  return;
}

#endif /* _QHG_SPINMAT_OPS_H */
