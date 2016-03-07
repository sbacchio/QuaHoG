#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <complex.h>
#include <qhg_idx.h>
#include <qhg_defs.h>
#include <qhg_types.h>
#include <qhg_alloc.h>
#include <qhg_correlator.h>
#include <qhg_prop_ops.h>
#include <qhg_prop_gammas.h>
#include <qhg_meson_udsc_defs.h>

static void
prop_bc(enum qhg_fermion_bc_time bc, _Complex double (*p)[NC*NS])
{
  switch(bc) {
  case PERIODIC:
    break;
  case ANTIPERIODIC:
    prop_scale(-1, p);
    break;
  }
  return;
}

qhg_correlator
qhg_mesons_udsc(qhg_spinor_field light[2][NS*NC], qhg_spinor_field strange[2][NS*NC], qhg_spinor_field charm[2][NS*NC],
		int source_coords[ND])
{
  qhg_lattice *lat = light[0][0].lat;
  qhg_correlator corr = qhg_correlator_init(NCHANNELS, lat);
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
    _Complex double Q[NFLAV][NS*NC][NS*NC];
    _Complex double C[NS*NC][NS*NC];    
    _Complex double W[NS*NC][NS*NC];
    _Complex double V[NS*NC][NS*NC];        

    prop_load(Q[0], light[0], v);
    prop_load(Q[1], light[1], v);
    
    prop_load(Q[2], strange[0], v);    
    prop_load(Q[3], strange[1], v);
    
    prop_load(Q[4], charm[0], v);
    prop_load(Q[5], charm[1], v);

    int t = v/lv3;
    int gt = t + t0;
    if(gt < tsrc) {
      prop_bc(light[0][0].bc, Q[0]);
      prop_bc(light[1][0].bc, Q[1]);
      prop_bc(strange[0][0].bc, Q[2]);
      prop_bc(strange[1][0].bc, Q[3]);
      prop_bc(charm[0][0].bc, Q[4]);
      prop_bc(charm[1][0].bc, Q[5]);
    }
       
    for(int igamma=0; igamma<NGAMMAS; igamma++) {
      for(int iflav0=0; iflav0<NFLAV; iflav0++) 
	for(int iflav1=0; iflav1<NFLAV; iflav1++) {
	  /*
	   * Sign flips are for consistency with libqcd
	   */
	  switch(igamma) {
	  case 0: /* 1 */
	    prop_mul_gd(C, Q[iflav0], Q[iflav1]);
	    corr.C[VGF(v, igamma, iflav0, iflav1)] = prop_trace(C);
	    break;
	  case 1: /* g5 */
	    prop_g5_G(W, Q[iflav0]);
	    prop_G_g5(V, W);	  
	    prop_mul_gd(C, V, Q[iflav1]);
	    corr.C[VGF(v, igamma, iflav0, iflav1)] = -prop_trace(C);
	    break;
	  case 2: /* gx */
	    prop_gx_G(W, Q[iflav0]);
	    prop_G_gx(V, W);	  
	    prop_mul_gd(C, V, Q[iflav1]);
	    corr.C[VGF(v, igamma, iflav0, iflav1)] = -prop_trace(C);
	    break;
	  case 3: /* gy */
	    prop_gy_G(W, Q[iflav0]);
	    prop_G_gy(V, W);	  
	    prop_mul_gd(C, V, Q[iflav1]);
	    corr.C[VGF(v, igamma, iflav0, iflav1)] = -prop_trace(C);
	    break;
	  case 4: /* gz */
	    prop_gz_G(W, Q[iflav0]);
	    prop_G_gz(V, W);	  
	    prop_mul_gd(C, V, Q[iflav1]);
	    corr.C[VGF(v, igamma, iflav0, iflav1)] = -prop_trace(C);
	    break;
	  case 5: /* gt */
	    prop_gt_G(W, Q[iflav0]);
	    prop_G_gt(V, W);	  
	    prop_mul_gd(C, V, Q[iflav1]);
	    corr.C[VGF(v, igamma, iflav0, iflav1)] = prop_trace(C);
	    break;
	  case 6: /* g5gx */
	    prop_g5gx_G(W, Q[iflav0]);
	    prop_G_g5gx(V, W);	  
	    prop_mul_gd(C, V, Q[iflav1]);
	    corr.C[VGF(v, igamma, iflav0, iflav1)] = -prop_trace(C);
	    break;
	  case 7: /* g5gy */
	    prop_g5gy_G(W, Q[iflav0]);
	    prop_G_g5gy(V, W);	  
	    prop_mul_gd(C, V, Q[iflav1]);
	    corr.C[VGF(v, igamma, iflav0, iflav1)] = -prop_trace(C);
	    break;
	  case 8: /* g5gz */
	    prop_g5gz_G(W, Q[iflav0]);
	    prop_G_g5gz(V, W);	  
	    prop_mul_gd(C, V, Q[iflav1]);
	    corr.C[VGF(v, igamma, iflav0, iflav1)] = -prop_trace(C);
	    break;
	  case 9: /* g5gt */
	    prop_g5gt_G(W, Q[iflav0]);
	    prop_G_g5gt(V, W);	  
	    prop_mul_gd(C, V, Q[iflav1]);
	    corr.C[VGF(v, igamma, iflav0, iflav1)] = prop_trace(C);
	    break;	  
	  }
	}
    }
  }
  
  corr.mom_list = NULL;
  return corr;
}
  
