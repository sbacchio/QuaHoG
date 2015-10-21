#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <complex.h>
#include <qhg_idx.h>
#include <qhg_defs.h>
#include <qhg_types.h>
#include <qhg_alloc.h>
#include <qhg_correlator.h>
#include <qhg_xchange_spinor.h>
#include <qhg_prop_gammas.h>
#include <qhg_prop_ops.h>
#include <qhg_su3_ops.h>
#include <qhg_io_utils.h>
#include <qhg_nn_thrp_defs.h>

#ifdef QHG_ACC
#include <openacc.h>
#endif

qhg_correlator
qhg_nn_thrp(qhg_spinor_field fwd[NS*NC], qhg_spinor_field bwd[NS*NC], qhg_gauge_field gf,
	    int source_coords[ND], qhg_thrp_nn_sink_params thrp_sink)
{  
  qhg_lattice *lat = fwd[0].lat;
  qhg_correlator corr = qhg_correlator_init(SITE_SIZE, lat);
  int lvol = lat->lvol;
  int **nn = lat->nn;
  for(int i=0; i<ND; i++)
    corr.origin[i] = source_coords[i];
  int tsrc = corr.origin[0];  
  int Lt = lat->dims[0];
  int lv3 = lat->lv3;
  int lt = lat->ldims[0];
  int t0 = lat->ldims[0]*lat->comms->proc_coords[0];  

  for(int i=0; i<NS*NC; i++) {
    qhg_xchange_spinor(bwd[i]);
    qhg_xchange_spinor(fwd[i]);
  }

  int *bvol = lat->bvol;
  int **evol = lat->evol;
  int lvol_bvol = lvol;
  for(int i=0; i<ND; i++)
    lvol_bvol += 2*bvol[i];
  
  int lvol_evol = lvol;
  for(int i=0; i<ND; i++) {
    lvol_evol += 2*bvol[i];
    for(int j=i+1; j<ND; j++)
      lvol_evol += 4*evol[i][j];
  }

  _Complex double F[NS*NC][NS*NC];
  _Complex double B[NS*NC][NS*NC];
  _Complex double T[NS*NC][NS*NC];
  _Complex double gF[NS*NC][NS*NC];          

  _Complex double gmFP[NS*NC][NS*NC];
  _Complex double gpFM[NS*NC][NS*NC];
  _Complex double gmF[NS*NC][NS*NC];
  _Complex double gpF[NS*NC][NS*NC];
  
  _Complex double FP[NS*NC][NS*NC];
  _Complex double FM[NS*NC][NS*NC];

  _Complex double gmuCnu[NC*NS][NC*NS];
  _Complex double gnuCmu[NC*NS][NC*NS];
  _Complex double g5simunuCku[NC*NS][NC*NS];
  _Complex double g5simukuCnu[NC*NS][NC*NS];
  
  _Complex double BP[NS*NC][NS*NC];
  _Complex double BM[NS*NC][NS*NC];
  
  _Complex double A0[NS*NC][NS*NC];
  _Complex double A1[NS*NC][NS*NC];
  _Complex double A2[NS*NC][NS*NC];
  _Complex double A3[NS*NC][NS*NC];

  _Complex double U0[NC][NC];
  _Complex double Um[NC][NC];

  _Complex double C[ND][NS*NC][NS*NC];  

  _Complex double fwd_tbc = fwd[0].bc[0];
  _Complex double bwd_tbc = bwd[0].bc[0];
  _Complex double *corrC = corr.C;

#ifndef QHG_ACC
#ifdef QHG_OMP
#pragma omp parallel for						\
  private(F, B, T, gF, gmFP, gpFM, gmF, gpF, FP, FM, BP, BM,		\
	  A0, A1, A2, A3, U0, Um, gmuCnu, gnuCmu, g5simukuCnu, g5simunuCku, \
	  C)
#endif
#endif  
#ifdef QHG_ACC
  qhg_spinor_field *fwd_d;
  qhg_spinor_field *bwd_d;
  fwd_d = (qhg_spinor_field *)acc_pcopyin(fwd, sizeof(qhg_spinor_field)*NC*NS);
  bwd_d = (qhg_spinor_field *)acc_pcopyin(bwd, sizeof(qhg_spinor_field)*NC*NS);
  for(int cs=0; cs<NC*NS; cs++) {
    _Complex double *field_d;
    field_d = (_Complex double *) acc_pcreate(fwd[cs].field, sizeof(_Complex double)*lvol_bvol*NC*NS);
    acc_update_device(fwd[cs].field, sizeof(_Complex double)*lvol_bvol*NC*NS);
    acc_memcpy_to_device(&fwd_d[cs].field, &field_d, sizeof(_Complex double *));
  }
  for(int cs=0; cs<NC*NS; cs++) {
    _Complex double *field_d;
    field_d = (_Complex double *) acc_pcreate(bwd[cs].field, sizeof(_Complex double)*lvol_bvol*NC*NS);
    acc_update_device(bwd[cs].field, sizeof(_Complex double)*lvol_bvol*NC*NS);
    acc_memcpy_to_device(&bwd_d[cs].field, &field_d, sizeof(_Complex double *));
  }
  
#pragma acc parallel loop copyout(corrC[0:(lvol*SITE_SIZE)])	\
  private(F, B, T, gF, fwd_tbc, bwd_tbc, tsrc)			\
  present(fwd,bwd)
#endif
  for(int v=0; v<lvol; v++) {
    prop_load(F, fwd, v);
    prop_load(B, bwd, v);
        
    int t = v/lv3;
    int gt = t + t0;
    if(gt < tsrc) {
      prop_scale(fwd_tbc, F);
      prop_scale(bwd_tbc, B);
    }

    int i = 0;
    prop_1_G(gF, F);
        
    prop_mul_gg(T, gF, B);    
    corrC[TIDX(v, i)] = prop_trace(T);

    /* /\* Local *\/ */
    /* for(int i=0; i<NLOC; i++) { */
    /*   switch(i) { */
    /* 	/\* Scalar *\/ */
    /*   case one: */
    /* 	prop_1_G(gF, F); */
    /* 	break; */
    /* 	/\* Pseudo-scalar *\/ */
    /*   case g5: */
    /* 	prop_g5_G(gF, F); */
    /* 	break; */
    /* 	/\* Vector *\/ */
    /*   case g0: */
    /* 	prop_g0_G(gF, F); */
    /* 	break; */
    /*   case gx: */
    /* 	prop_gx_G(gF, F); */
    /* 	break; */
    /*   case gy: */
    /* 	prop_gy_G(gF, F); */
    /* 	break; */
    /*   case gz: */
    /* 	prop_gz_G(gF, F); */
    /* 	break; */
    /* 	/\* Axial *\/ */
    /*   case g5g0: */
    /* 	prop_g5g0_G(gF, F); */
    /* 	break; */
    /*   case g5gx: */
    /* 	prop_g5gx_G(gF, F); */
    /* 	break; */
    /*   case g5gy: */
    /* 	prop_g5gy_G(gF, F); */
    /* 	break; */
    /*   case g5gz: */
    /* 	prop_g5gz_G(gF, F); */
    /* 	break; */
    /* 	/\* Tensor *\/ */
    /*   case g5si0x: */
    /* 	prop_g5si0x_G(gF, F); */
    /* 	break; */
    /*   case g5si0y: */
    /* 	prop_g5si0y_G(gF, F); */
    /* 	break; */
    /*   case g5si0z: */
    /* 	prop_g5si0z_G(gF, F); */
    /* 	break; */
    /*   case g5sixy: */
    /* 	prop_g5sixy_G(gF, F); */
    /* 	break; */
    /*   case g5sixz: */
    /* 	prop_g5sixz_G(gF, F); */
    /* 	break; */
    /*   case g5siyz: */
    /* 	prop_g5siyz_G(gF, F); */
    /* 	break; */
    /*   } */
    /*   prop_mul_gg(T, gF, B); */
    /*   corr.C[TIDX(v, i)] = prop_trace(T); */
    /* } */

    /* /\* Noether *\/ */
    /* for(int i=NLOC; i<NLOC+NNOE; i++) { */
    /*   int mu = i - NLOC; */
    /*   int vp = nn[mu][v]; */
    /*   int vm = nn[mu+ND][v]; */
      
    /*   su3_load(U0, gf, v, mu); */
    /*   su3_load(Um, gf, vm, mu); */
      
    /*   prop_load(FP, fwd, vp); */
    /*   prop_load(FM, fwd, vm); */
    /*   prop_load(BP, bwd, vp); */
    /*   prop_load(BM, bwd, vm); */

    /*   if(mu == 0) { */
    /* 	int t = v/lv3; */
    /* 	int gt0 = t + t0;			/\* Global absolute insertion time *\/ */
    /* 	int gdt = ((gt0 - tsrc) + Lt) % Lt;     /\* Global sink-source separation  *\/ */
    /* 	if(gdt + tsrc + 1 >= Lt) { */
    /* 	  prop_scale(fwd[0].bc[0], FP); */
    /* 	  prop_scale(bwd[0].bc[0], BP); */
    /* 	} */
    /*   } */

    /*   if(mu == 0) { */
    /* 	int t = v/lv3; */
    /* 	int gt0 = t + t0;			/\* Global absolute insertion time *\/ */
    /* 	int gdt = ((gt0 - tsrc) + Lt) % Lt;     /\* Global sink-source separation  *\/ */
    /* 	if(gdt + tsrc - 1 < 0 || gdt + tsrc - 1 >= Lt) { */
    /* 	  prop_scale(fwd[0].bc[0], FM); */
    /* 	  prop_scale(bwd[0].bc[0], BM); */
    /* 	} */
    /*   } */
      
    /*   switch(i) { */
    /*   case noe_g0: */
    /* 	prop_1pg0_G(gpFM, FM); */
    /* 	prop_1mg0_G(gmFP, FP); */
    /* 	prop_1pg0_G(gpF, F); */
    /* 	prop_1mg0_G(gmF, F); */
    /* 	break; */
    /*   case noe_gx: */
    /* 	prop_1pgx_G(gpFM, FM); */
    /* 	prop_1mgx_G(gmFP, FP); */
    /* 	prop_1pgx_G(gpF, F); */
    /* 	prop_1mgx_G(gmF, F); */
    /* 	break; */
    /*   case noe_gy: */
    /* 	prop_1pgy_G(gpFM, FM); */
    /* 	prop_1mgy_G(gmFP, FP); */
    /* 	prop_1pgy_G(gpF, F); */
    /* 	prop_1mgy_G(gmF, F); */
    /* 	break; */
    /*   case noe_gz: */
    /* 	prop_1pgz_G(gpFM, FM); */
    /* 	prop_1mgz_G(gmFP, FP); */
    /* 	prop_1pgz_G(gpF, F); */
    /* 	prop_1mgz_G(gmF, F); */
    /* 	break; */
    /*   } */
    /*   prop_mul_su3_D_G(T, gpF, U0); */
    /*   prop_mul_gg(A0, T, BP); */

    /*   prop_mul_su3_D_G(T, gpFM, Um); */
    /*   prop_mul_gg(A1, T, B); */

    /*   prop_mul_su3_U_G(T, gmFP, U0); */
    /*   prop_mul_gg(A2, T, B); */
      
    /*   prop_mul_su3_U_G(T, gmF, Um); */
    /*   prop_mul_gg(A3, T, BM); */

    /*   prop_add_gg(T, A0, A1); */
    /*   prop_meq_g(T, A2); */
    /*   prop_meq_g(T, A3); */
    /*   corr.C[TIDX(v, i)] = prop_trace(T)*0.25; */
    /* } */

    /* /\* Compute gamma-less derivative *\/ */
    /* for(int mu=0; mu<ND; mu++) {       */
    /*   int vp = nn[mu][v]; */
    /*   int vm = nn[mu+ND][v]; */

    /*   su3_load(U0, gf, v, mu); */
    /*   su3_load(Um, gf, vm, mu); */
      
    /*   prop_load(FP, fwd, vp); */
    /*   prop_load(FM, fwd, vm); */
    /*   prop_load(BP, bwd, vp); */
    /*   prop_load(BM, bwd, vm); */


    /*   if(mu == 0) { */
    /* 	int t = v/lv3; */
    /* 	int gt0 = t + t0;			/\* Global absolute insertion time *\/ */
    /* 	int gdt = ((gt0 - tsrc) + Lt) % Lt;     /\* Global sink-source separation  *\/ */
    /* 	if(gdt + tsrc + 1 >= Lt) { */
    /* 	  prop_scale(fwd[0].bc[0], FP); */
    /* 	  prop_scale(bwd[0].bc[0], BP); */
    /* 	} */
    /*   } */

    /*   if(mu == 0) { */
    /* 	int t = v/lv3; */
    /* 	int gt0 = t + t0;			/\* Global absolute insertion time *\/ */
    /* 	int gdt = ((gt0 - tsrc) + Lt) % Lt;     /\* Global sink-source separation  *\/ */
    /* 	if(gdt + tsrc - 1 < 0 || gdt + tsrc - 1 >= Lt) { */
    /* 	  prop_scale(fwd[0].bc[0], FM); */
    /* 	  prop_scale(bwd[0].bc[0], BM); */
    /* 	} */
    /*   } */
      
    /*   prop_mul_su3_U_G(T, FP, U0); */
    /*   prop_mul_gg(A0, T, B);       */

    /*   prop_mul_su3_D_G(T, FM, Um); */
    /*   prop_mul_gg(A1, T, B);       */

    /*   prop_mul_su3_D_G(T, F, U0); */
    /*   prop_mul_gg(A2, T, BP);       */

    /*   prop_mul_su3_U_G(T, F, Um); */
    /*   prop_mul_gg(A3, T, BM);       */
    /*   prop_sub_gg(C[mu], A0, A1); */
    /*   prop_meq_g(C[mu], A2); */
    /*   prop_peq_g(C[mu], A3); */
    /* } */

    /* /\* Vector derivative *\/ */
    /* for(int i=NLOC+NNOE; i<NLOC+NNOE+NVDER; i++) { */
    /*   int mu, nu; */
    /*   switch(i) { */
    /*   case der_g0_0: */
    /* 	mu = 0; */
    /* 	nu = 0; */
    /* 	prop_g0_G(gnuCmu, C[mu]);	 */
    /* 	prop_g0_G(gmuCnu, C[nu]); */
    /* 	break; */
    /*   case der_gx_0: */
    /* 	mu = 0; */
    /* 	nu = 1; */
    /* 	prop_gx_G(gnuCmu, C[mu]);	 */
    /* 	prop_g0_G(gmuCnu, C[nu]); */
    /* 	break; */
    /*   case der_gy_0: */
    /* 	mu = 0; */
    /* 	nu = 2; */
    /* 	prop_gy_G(gnuCmu, C[mu]);	 */
    /* 	prop_g0_G(gmuCnu, C[nu]); */
    /* 	break; */
    /*   case der_gz_0: */
    /* 	mu = 0; */
    /* 	nu = 3; */
    /* 	prop_gz_G(gnuCmu, C[mu]);	 */
    /* 	prop_g0_G(gmuCnu, C[nu]); */
    /* 	break; */
    /* 	// */
    /*   case der_gx_x: */
    /* 	mu = 1; */
    /* 	nu = 1; */
    /* 	prop_gx_G(gnuCmu, C[mu]);	 */
    /* 	prop_gx_G(gmuCnu, C[nu]); */
    /* 	break; */
    /*   case der_gy_x: */
    /* 	mu = 1; */
    /* 	nu = 2; */
    /* 	prop_gy_G(gnuCmu, C[mu]);	 */
    /* 	prop_gx_G(gmuCnu, C[nu]); */
    /* 	break; */
    /*   case der_gz_x: */
    /* 	mu = 1; */
    /* 	nu = 3; */
    /* 	prop_gz_G(gnuCmu, C[mu]);	 */
    /* 	prop_gx_G(gmuCnu, C[nu]); */
    /* 	break; */
    /* 	// */
    /*   case der_gy_y: */
    /* 	mu = 2; */
    /* 	nu = 2; */
    /* 	prop_gy_G(gnuCmu, C[mu]);	 */
    /* 	prop_gy_G(gmuCnu, C[nu]); */
    /* 	break; */
    /*   case der_gz_y: */
    /* 	mu = 2; */
    /* 	nu = 3; */
    /* 	prop_gz_G(gnuCmu, C[mu]);	 */
    /* 	prop_gy_G(gmuCnu, C[nu]); */
    /* 	break; */
    /* 	// */
    /*   case der_gz_z: */
    /* 	mu = 3; */
    /* 	nu = 3; */
    /* 	prop_gz_G(gnuCmu, C[mu]);	 */
    /* 	prop_gz_G(gmuCnu, C[nu]); */
    /* 	break;	   */
    /*   } */
      
    /*   prop_add_gg(T, gmuCnu, gnuCmu); */
    /*   corr.C[TIDX(v, i)] = prop_trace(T)*0.125; */
    /* } */

    /* /\* Axial derivative, symmetric *\/ */
    /* for(int i=NLOC+NNOE+NVDER; i<NLOC+NNOE+NVDER+NADSY; i++) { */
    /*   int mu, nu; */
    /*   switch(i) { */
    /*   case der_g5g0_0_sym: */
    /* 	mu = 0; */
    /* 	nu = 0; */
    /* 	prop_g5g0_G(gnuCmu, C[mu]);	 */
    /* 	prop_g5g0_G(gmuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5gx_0_sym: */
    /* 	mu = 0; */
    /* 	nu = 1; */
    /* 	prop_g5gx_G(gnuCmu, C[mu]);	 */
    /* 	prop_g5g0_G(gmuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5gy_0_sym: */
    /* 	mu = 0; */
    /* 	nu = 2; */
    /* 	prop_g5gy_G(gnuCmu, C[mu]);	 */
    /* 	prop_g5g0_G(gmuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5gz_0_sym: */
    /* 	mu = 0; */
    /* 	nu = 3; */
    /* 	prop_g5gz_G(gnuCmu, C[mu]);	 */
    /* 	prop_g5g0_G(gmuCnu, C[nu]); */
    /* 	break; */
    /* 	// */
    /*   case der_g5gx_x_sym: */
    /* 	mu = 1; */
    /* 	nu = 1; */
    /* 	prop_g5gx_G(gnuCmu, C[mu]);	 */
    /* 	prop_g5gx_G(gmuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5gy_x_sym: */
    /* 	mu = 1; */
    /* 	nu = 2; */
    /* 	prop_g5gy_G(gnuCmu, C[mu]);	 */
    /* 	prop_g5gx_G(gmuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5gz_x_sym: */
    /* 	mu = 1; */
    /* 	nu = 3; */
    /* 	prop_g5gz_G(gnuCmu, C[mu]);	 */
    /* 	prop_g5gx_G(gmuCnu, C[nu]); */
    /* 	break; */
    /* 	// */
    /*   case der_g5gy_y_sym: */
    /* 	mu = 2; */
    /* 	nu = 2; */
    /* 	prop_g5gy_G(gnuCmu, C[mu]);	 */
    /* 	prop_g5gy_G(gmuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5gz_y_sym: */
    /* 	mu = 2; */
    /* 	nu = 3; */
    /* 	prop_g5gz_G(gnuCmu, C[mu]);	 */
    /* 	prop_g5gy_G(gmuCnu, C[nu]); */
    /* 	break; */
    /* 	// */
    /*   case der_g5gz_z_sym: */
    /* 	mu = 3; */
    /* 	nu = 3; */
    /* 	prop_g5gz_G(gnuCmu, C[mu]);	 */
    /* 	prop_g5gz_G(gmuCnu, C[nu]); */
    /* 	break;	   */
    /*   } */
      
    /*   prop_add_gg(T, gmuCnu, gnuCmu); */
    /*   corr.C[TIDX(v, i)] = prop_trace(T)*0.125; */
    /* } */

    /* /\* Axial derivative, anti-symmetric *\/ */
    /* for(int i=NLOC+NNOE+NVDER+NADSY; i<NLOC+NNOE+NVDER+NADSY+NADAS; i++) { */
    /*   int mu, nu; */
    /*   switch(i) { */
    /*   case der_g5gx_0_asy: */
    /* 	mu = 0; */
    /* 	nu = 1; */
    /* 	prop_g5gx_G(gnuCmu, C[mu]);	 */
    /* 	prop_g5g0_G(gmuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5gy_0_asy: */
    /* 	mu = 0; */
    /* 	nu = 2; */
    /* 	prop_g5gy_G(gnuCmu, C[mu]);	 */
    /* 	prop_g5g0_G(gmuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5gz_0_asy: */
    /* 	mu = 0; */
    /* 	nu = 3; */
    /* 	prop_g5gz_G(gnuCmu, C[mu]);	 */
    /* 	prop_g5g0_G(gmuCnu, C[nu]); */
    /* 	break; */
    /* 	// */
    /*   case der_g5gy_x_asy: */
    /* 	mu = 1; */
    /* 	nu = 2; */
    /* 	prop_g5gy_G(gnuCmu, C[mu]);	 */
    /* 	prop_g5gx_G(gmuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5gz_x_asy: */
    /* 	mu = 1; */
    /* 	nu = 3; */
    /* 	prop_g5gz_G(gnuCmu, C[mu]);	 */
    /* 	prop_g5gx_G(gmuCnu, C[nu]); */
    /* 	break; */
    /* 	// */
    /*   case der_g5gz_y_asy: */
    /* 	mu = 2; */
    /* 	nu = 3; */
    /* 	prop_g5gz_G(gnuCmu, C[mu]);	 */
    /* 	prop_g5gy_G(gmuCnu, C[nu]); */
    /* 	break; */
    /*   } */
      
    /*   prop_sub_gg(T, gnuCmu, gmuCnu); */
    /*   corr.C[TIDX(v, i)] = prop_trace(T)*0.125; */
    /* } */

    /* /\* Tensor derivative, symmetric *\/ */
    /* for(int i=NLOC+NNOE+NVDER+NADSY+NADAS; i<NLOC+NNOE+NVDER+NADSY+NADAS+NTDER; i++) { */
    /*   int mu, nu, ku; */
    /*   switch(i) { */
    /*   case der_g5si00_x: */
    /* 	mu = 0; */
    /* 	nu = 0; */
    /* 	ku = 1; */
    /* 	prop_zero(g5simunuCku); */
    /* 	prop_g5si0x_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5si00_y: */
    /* 	mu = 0; */
    /* 	nu = 0; */
    /* 	ku = 2; */
    /* 	prop_zero(g5simunuCku); */
    /* 	prop_g5si0y_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5si00_z: */
    /* 	mu = 0; */
    /* 	nu = 0; */
    /* 	ku = 3; */
    /* 	prop_zero(g5simunuCku); */
    /* 	prop_g5si0z_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5si0x_x: */
    /* 	mu = 0; */
    /* 	nu = 1; */
    /* 	ku = 1; */
    /* 	prop_g5si0x_G(g5simunuCku, C[ku]); */
    /* 	prop_g5si0x_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5si0x_y: */
    /* 	mu = 0; */
    /* 	nu = 1; */
    /* 	ku = 2; */
    /* 	prop_g5si0x_G(g5simunuCku, C[ku]); */
    /* 	prop_g5si0y_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5si0x_z: */
    /* 	mu = 0; */
    /* 	nu = 1; */
    /* 	ku = 3; */
    /* 	prop_g5si0x_G(g5simunuCku, C[ku]); */
    /* 	prop_g5si0z_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5si0y_y: */
    /* 	mu = 0; */
    /* 	nu = 2; */
    /* 	ku = 2; */
    /* 	prop_g5si0y_G(g5simunuCku, C[ku]); */
    /* 	prop_g5si0y_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5si0y_z: */
    /* 	mu = 0; */
    /* 	nu = 2; */
    /* 	ku = 3; */
    /* 	prop_g5si0y_G(g5simunuCku, C[ku]); */
    /* 	prop_g5si0z_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5si0z_z: */
    /* 	mu = 0; */
    /* 	nu = 3; */
    /* 	ku = 3; */
    /* 	prop_g5si0z_G(g5simunuCku, C[ku]); */
    /* 	prop_g5si0z_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /* 	// */
    /*   case der_g5six0_0: */
    /* 	mu = 1; */
    /* 	nu = 0; */
    /* 	ku = 0; */
    /* 	prop_g5six0_G(g5simunuCku, C[ku]); */
    /* 	prop_g5six0_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5six0_x: */
    /* 	mu = 1; */
    /* 	nu = 0; */
    /* 	ku = 1; */
    /* 	prop_g5six0_G(g5simunuCku, C[ku]); */
    /* 	prop_zero(g5simukuCnu); */
    /* 	break; */
    /*   case der_g5six0_y: */
    /* 	mu = 1; */
    /* 	nu = 0; */
    /* 	ku = 2; */
    /* 	prop_g5six0_G(g5simunuCku, C[ku]); */
    /* 	prop_g5sixy_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5six0_z: */
    /* 	mu = 1; */
    /* 	nu = 0; */
    /* 	ku = 3; */
    /* 	prop_g5six0_G(g5simunuCku, C[ku]); */
    /* 	prop_g5sixz_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5sixx_y: */
    /* 	mu = 1; */
    /* 	nu = 1; */
    /* 	ku = 2; */
    /* 	prop_zero(g5simunuCku); */
    /* 	prop_g5sixy_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5sixx_z: */
    /* 	mu = 1; */
    /* 	nu = 1; */
    /* 	ku = 3; */
    /* 	prop_zero(g5simunuCku); */
    /* 	prop_g5sixz_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5sixy_y: */
    /* 	mu = 1; */
    /* 	nu = 2; */
    /* 	ku = 2; */
    /* 	prop_g5sixy_G(g5simunuCku, C[ku]); */
    /* 	prop_g5sixy_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5sixy_z: */
    /* 	mu = 1; */
    /* 	nu = 2; */
    /* 	ku = 3; */
    /* 	prop_g5sixy_G(g5simunuCku, C[ku]); */
    /* 	prop_g5sixz_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5sixz_z: */
    /* 	mu = 1; */
    /* 	nu = 3; */
    /* 	ku = 3; */
    /* 	prop_g5sixz_G(g5simunuCku, C[ku]); */
    /* 	prop_g5sixz_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /* 	// */
    /*   case der_g5siy0_0: */
    /* 	mu = 2; */
    /* 	nu = 0; */
    /* 	ku = 0; */
    /* 	prop_g5siy0_G(g5simunuCku, C[ku]); */
    /* 	prop_g5siy0_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5siy0_x: */
    /* 	mu = 2; */
    /* 	nu = 0; */
    /* 	ku = 1; */
    /* 	prop_g5siy0_G(g5simunuCku, C[ku]); */
    /* 	prop_g5siyx_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5siy0_y: */
    /* 	mu = 2; */
    /* 	nu = 0; */
    /* 	ku = 2; */
    /* 	prop_g5siy0_G(g5simunuCku, C[ku]); */
    /* 	prop_zero(g5simukuCnu); */
    /* 	break; */
    /*   case der_g5siy0_z: */
    /* 	mu = 2; */
    /* 	nu = 0; */
    /* 	ku = 3; */
    /* 	prop_g5siy0_G(g5simunuCku, C[ku]); */
    /* 	prop_g5siyz_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5siyx_x: */
    /* 	mu = 2; */
    /* 	nu = 1; */
    /* 	ku = 1; */
    /* 	prop_g5siyx_G(g5simunuCku, C[ku]); */
    /* 	prop_g5siyx_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5siyx_y: */
    /* 	mu = 2; */
    /* 	nu = 1; */
    /* 	ku = 2; */
    /* 	prop_g5siyx_G(g5simunuCku, C[ku]); */
    /* 	prop_zero(g5simukuCnu); */
    /* 	break; */
    /*   case der_g5siyx_z: */
    /* 	mu = 2; */
    /* 	nu = 1; */
    /* 	ku = 3; */
    /* 	prop_g5siyx_G(g5simunuCku, C[ku]); */
    /* 	prop_g5siyz_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5siyy_z: */
    /* 	mu = 2; */
    /* 	nu = 2; */
    /* 	ku = 3; */
    /* 	prop_zero(g5simunuCku); */
    /* 	prop_g5siyz_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5siyz_z: */
    /* 	mu = 2; */
    /* 	nu = 3; */
    /* 	ku = 3; */
    /* 	prop_g5siyz_G(g5simunuCku, C[ku]); */
    /* 	prop_g5siyz_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /* 	// */
    /*   case der_g5siz0_0: */
    /* 	mu = 3; */
    /* 	nu = 0; */
    /* 	ku = 0; */
    /* 	prop_g5siz0_G(g5simunuCku, C[ku]); */
    /* 	prop_g5siz0_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5siz0_x: */
    /* 	mu = 3; */
    /* 	nu = 0; */
    /* 	ku = 1; */
    /* 	prop_g5siz0_G(g5simunuCku, C[ku]); */
    /* 	prop_g5sizx_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5siz0_y: */
    /* 	mu = 3; */
    /* 	nu = 0; */
    /* 	ku = 2; */
    /* 	prop_g5siz0_G(g5simunuCku, C[ku]); */
    /* 	prop_g5sizy_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5siz0_z: */
    /* 	mu = 3; */
    /* 	nu = 0; */
    /* 	ku = 3; */
    /* 	prop_g5siz0_G(g5simunuCku, C[ku]); */
    /* 	prop_zero(g5simukuCnu); */
    /* 	break; */
    /*   case der_g5sizx_x: */
    /* 	mu = 3; */
    /* 	nu = 1; */
    /* 	ku = 1; */
    /* 	prop_g5sizx_G(g5simunuCku, C[ku]); */
    /* 	prop_g5sizx_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5sizx_y: */
    /* 	mu = 3; */
    /* 	nu = 1; */
    /* 	ku = 2; */
    /* 	prop_g5sizx_G(g5simunuCku, C[ku]); */
    /* 	prop_g5sizy_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5sizx_z: */
    /* 	mu = 3; */
    /* 	nu = 1; */
    /* 	ku = 3; */
    /* 	prop_g5sizx_G(g5simunuCku, C[ku]); */
    /* 	prop_zero(g5simukuCnu); */
    /* 	break; */
    /*   case der_g5sizy_y: */
    /* 	mu = 3; */
    /* 	nu = 2; */
    /* 	ku = 2; */
    /* 	prop_g5sizy_G(g5simunuCku, C[ku]); */
    /* 	prop_g5sizy_G(g5simukuCnu, C[nu]); */
    /* 	break; */
    /*   case der_g5sizy_z: */
    /* 	mu = 3; */
    /* 	nu = 2; */
    /* 	ku = 3; */
    /* 	prop_g5sizy_G(g5simunuCku, C[ku]); */
    /* 	prop_zero(g5simukuCnu); */
    /* 	break; */
    /*   } */
      
    /*   prop_add_gg(T, g5simunuCku, g5simukuCnu); */
    /*   corr.C[TIDX(v, i)] = prop_trace(T)*0.125; */
    /* } */
  }

  corr.cutoff[0] = (source_coords[0] + thrp_sink.dt) % Lt;
  corr.mom_list = NULL;

  return corr;
}
