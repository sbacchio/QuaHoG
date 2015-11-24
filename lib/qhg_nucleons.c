#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <complex.h>
#include <qhg_idx.h>
#include <qhg_defs.h>
#include <qhg_types.h>
#include <qhg_alloc.h>
#include <qhg_correlator.h>
#include <qhg_prop_gammas.h>
#include <qhg_prop_ops.h>
#include <qhg_spinmat_ops.h>
#include <qhg_nucleon_defs.h>

static void
nucleon(_Complex double corr[4][NS][NS], _Complex double (*x)[NC*NS],
	_Complex double (*y)[NC*NS])
{
  _Complex double z[NS*NC][NS*NC];
  _Complex double aux[NS*NC][NS*NC];  
  _Complex double A[NS][NS], B[NS][NS];
      
  /*
   * \chi_1 - to - \chi_1
   */
  _Complex double Cg5xCg5[NS*NC][NS*NC];
  prop_Cg5_G(aux, x);
  prop_G_Cg5(Cg5xCg5, aux);
  prop_contract_02(z, Cg5xCg5, y);       
  prop_color_transpose(z);
  prop_mul_gtr(aux, y, z);
  prop_color_trace(A, aux);
  prop_mul_gg(aux, y, z);
  prop_color_trace(B, aux);
  spinmat_add_ss(corr[0], A, B);
  
  /*
   * \chi_1 - to - \chi_2
   */
  _Complex double yg5[NS*NC][NS*NC];
  prop_G_g5(yg5, y);
  _Complex double Cg5xC[NS*NC][NS*NC];
  prop_Cg5_G(aux, x);
  prop_G_C(Cg5xC, aux);
  prop_contract_13(z, Cg5xC, y);
  prop_transpose(z);
  prop_mul_gtr(aux, yg5, z);
  prop_color_trace(A, aux);
  prop_mul_gg(aux, z, yg5);
  prop_color_trace(B, aux);
  spinmat_add_ss(corr[1], A, B);

  /*
   * \chi_2 - to - \chi_1
   */
  _Complex double g5y[NS*NC][NS*NC];
  prop_g5_G(g5y, y);
  _Complex double CxCg5[NS*NC][NS*NC];
  prop_G_Cg5(aux, x);
  prop_C_G(CxCg5, aux);
  prop_contract_02(z, CxCg5, y);
  prop_color_transpose(z);
  prop_mul_gtr(aux, g5y, z);
  prop_color_trace(A, aux);
  prop_mul_gg(aux, g5y, z);
  prop_color_trace(B, aux);
  spinmat_add_ss(corr[2], A, B);

  /*
   * \chi_2 - to - \chi_2
   */
  _Complex double g5yg5[NS*NC][NS*NC];
  prop_g5_G(g5yg5, yg5);
  _Complex double CxC[NS*NC][NS*NC];
  _Complex double w[NS*NC][NS*NC];
  prop_G_C(aux, x);
  prop_C_G(CxC, aux);
  prop_contract_02(z, CxC, y);
  prop_contract_02(w, CxC, yg5);
  prop_color_transpose(z);
  prop_color_transpose(w);
  prop_mul_gtr(aux, g5yg5, z);
  prop_color_trace(A, aux);
  prop_mul_gg(aux, g5y, w);
  prop_color_trace(B, aux);
  spinmat_add_ss(corr[3], A, B);
  
  return;
}

qhg_correlator
qhg_nucleons(qhg_spinor_field sp_u[NS*NC], qhg_spinor_field sp_d[NS*NC],
	     int source_coords[ND])
{
  qhg_lattice *lat = sp_u[0].lat;
  qhg_correlator corr = qhg_correlator_init(SITE_SIZE, lat);  
  for(int i=0; i<ND; i++)
    corr.origin[i] = source_coords[i];
  int tsrc = corr.origin[0];  
  int lvol = lat->lvol;
  int lv3 = lat->lv3;
  int lt = lat->ldims[0];
  int t0 = lat->ldims[0]*lat->comms->proc_coords[0];

#ifdef QHG_OMP
#pragma omp parallel for
#endif
  for(int v=0; v<lvol; v++) {
    _Complex double U[NS*NC][NS*NC];
    _Complex double D[NS*NC][NS*NC];    
    _Complex double C[NS*NC][NS*NC];    
    _Complex double W[NS*NC][NS*NC];
    _Complex double V[NS*NC][NS*NC];        

    prop_load(U, sp_u, v);
    prop_load(D, sp_d, v);    

    int t = v/lv3;
    int gt = t + t0;
    if(gt < tsrc) {
      prop_scale(sp_u[0].bc[0], U);
      prop_scale(sp_d[0].bc[0], D);
    }
    
    _Complex double (*P[2])[NC*NS] = {U, D};
    for(int iflav=0; iflav<NFLAV; iflav++) {
      _Complex double c[4][NS][NS];
      nucleon(c, P[iflav], P[(iflav+1)%2]);
      for(int i=0; i<4; i++)
	for(int s0=0; s0<NS; s0++)
	  for(int s1=0; s1<NS; s1++) {
	    corr.C[NIDX(v, iflav, i, s0, s1)] = c[i][s0][s1];
	  }
    }
  }
  
  corr.mom_list = NULL;
  return corr;
}
