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
#include <lines_types.h>
#include <lines_utils.h>

#define NLOC 16			// number of local operator directions
#define NNOE 4			// number of noether operator directions
#define NVDER ((4*(4+1))/2)	// number of vector-derivative directions
#define NADSY ((4*(4+1))/2)	// number of axial-derivative directions (symmetric)
#define NADAS ((4*(4-1))/2)	// number of axial-derivative directions (anti-symmetric)
#define NTDER (4*((4+1)*4)/2-4) // number of tensor-derivative directions [ (si_{mu,nu} D_ku + si_{mu,ku} D_nu) with nu <= ku and ommit ku == mu == nu which is zero]
#define NCHAN (NLOC+NNOE+NVDER+NADSY+NADAS+NTDER)
#define SITE_SIZE (NCHAN)
#define TIDX(v, ch) (ch + v*NCHAN)

static char chan_tags[NCHAN][256] = {
  "=loc:1=\0", "=loc:g5=\0",
  "=loc:g0=\0", "=loc:gx=\0", "=loc:gy=\0", "=loc:gz=\0",
  "=loc:g5g0=\0", "=loc:g5gx=\0", "=loc:g5gy=\0", "=loc:g5gz=\0",
  "=loc:g5si0x=\0", "=loc:g5si0y=\0", "=loc:g5si0z=\0", "=loc:g5sixy=\0", "=loc:g5sixz=\0", "=loc:g5siyz=\0",
  //
  "=noe:g0=\0", "=noe:gx=\0", "=noe:gy=\0", "=noe:gz=\0",
  //
  "=der:g0D0:sym=\0", "=der:gxD0:sym=\0", "=der:gyD0:sym=\0", "=der:gzD0:sym=\0",
  "=der:gxDx:sym=\0", "=der:gyDx:sym=\0", "=der:gzDx:sym=\0",
  "=der:gyDy:sym=\0", "=der:gzDy:sym=\0",
  "=der:gzDz:sym=\0",
  //
  "=der:g5g0D0:sym=\0", "=der:g5gxD0:sym=\0", "=der:g5gyD0:sym=\0", "=der:g5gzD0:sym=\0",
  "=der:g5gxDx:sym=\0", "=der:g5gyDx:sym=\0", "=der:g5gzDx:sym=\0",
  "=der:g5gyDy:sym=\0", "=der:g5gzDy:sym=\0",
  "=der:g5gzDz:sym=\0",
  //
  "=der:g5gxD0:asy=\0", "=der:g5gyD0:asy=\0", "=der:g5gzD0:asy=\0",
  "=der:g5gyDx:asy=\0", "=der:g5gzDx:asy=\0",
  "=der:g5gzDy:asy=\0",
  //
  "=der:g5si00Dx:sym=\0", "=der:g5si00Dy:sym=\0", "=der:g5si00Dz:sym=\0",
  "=der:g5si0xDx:sym=\0", "=der:g5si0xDy:sym=\0", "=der:g5si0xDz:sym=\0",
  "=der:g5si0yDy:sym=\0", "=der:g5si0yDz:sym=\0",
  "=der:g5si0zDz:sym=\0",
  //
  "=der:g5six0D0:sym=\0", "=der:g5six0Dx:sym=\0", "=der:g5six0Dy:sym=\0", "=der:g5six0Dz:sym=\0",
  "=der:g5sixxDy:sym=\0", "=der:g5sixxDz:sym=\0",
  "=der:g5sixyDy:sym=\0", "=der:g5sixyDz:sym=\0",
  "=der:g5sixzDz:sym=\0",
  //
  "=der:g5siy0D0:sym=\0", "=der:g5siy0Dx:sym=\0", "=der:g5siy0Dy:sym=\0", "=der:g5siy0Dz:sym=\0",
  "=der:g5siyxDx:sym=\0", "=der:g5siyxDy:sym=\0", "=der:g5siyxDz:sym=\0",
  "=der:g5siyyDz:sym=\0",
  "=der:g5siyzDz:sym=\0",  
  //
  "=der:g5siz0D0:sym=\0", "=der:g5siz0Dx:sym=\0", "=der:g5siz0Dy:sym=\0", "=der:g5siz0Dz:sym=\0",
  "=der:g5sizxDx:sym=\0", "=der:g5sizxDy:sym=\0", "=der:g5sizxDz:sym=\0",
  "=der:g5sizyDy:sym=\0", "=der:g5sizyDz:sym=\0",
};

static enum {
  one, g5, g0, gx, gy, gz,
  g5g0, g5gx, g5gy, g5gz,
  g5si0x, g5si0y, g5si0z, g5sixy, g5sixz, g5siyz,
  noe_g0, noe_gx, noe_gy, noe_gz,
  der_g0_0, der_gx_0, der_gy_0, der_gz_0,
  der_gx_x, der_gy_x, der_gz_x,   
  der_gy_y, der_gz_y,
  der_gz_z,     
  der_g5g0_0_sym, der_g5gx_0_sym, der_g5gy_0_sym, der_g5gz_0_sym,
  der_g5gx_x_sym, der_g5gy_x_sym, der_g5gz_x_sym,   
  der_g5gy_y_sym, der_g5gz_y_sym,
  der_g5gz_z_sym,     
  der_g5gx_0_asy, der_g5gy_0_asy, der_g5gz_0_asy,
  der_g5gy_x_asy, der_g5gz_x_asy,   
  der_g5gz_y_asy,

  der_g5si00_x, der_g5si00_y, der_g5si00_z,
  der_g5si0x_x, der_g5si0x_y, der_g5si0x_z,
  der_g5si0y_y, der_g5si0y_z,
  der_g5si0z_z,

  der_g5six0_0, der_g5six0_x, der_g5six0_y, der_g5six0_z,
  der_g5sixx_y, der_g5sixx_z,
  der_g5sixy_y, der_g5sixy_z,
  der_g5sixz_z,

  der_g5siy0_0, der_g5siy0_x, der_g5siy0_y, der_g5siy0_z,
  der_g5siyx_x, der_g5siyx_y, der_g5siyx_z,
  der_g5siyy_z,
  der_g5siyz_z,  

  der_g5siz0_0, der_g5siz0_x, der_g5siz0_y, der_g5siz0_z,
  der_g5sizx_x, der_g5sizx_y, der_g5sizx_z,
  der_g5sizy_y, der_g5sizy_z,
} channels;

qhg_correlator
qhg_nn_thrp(qhg_spinor_field fwd[NS*NC], qhg_spinor_field bwd[NS*NC], qhg_gauge_field gf,
	    int source_coords[ND], qhg_thrp_nn_sink_params thrp_sink)
{  
  qhg_lattice *lat = fwd[0].lat;
  qhg_correlator corr = qhg_correlator_init(SITE_SIZE, lat);
  int lvol = lat->lvol;
  int **nn = lat->nn;

  for(int i=0; i<NS*NC; i++) {
    qhg_xchange_spinor(bwd[i]);
    qhg_xchange_spinor(fwd[i]);
  }
  
  for(int v=0; v<lvol; v++) {
    _Complex double F[NS*NC][NS*NC];
    _Complex double B[NS*NC][NS*NC];
    _Complex double T[NS*NC][NS*NC];
    prop_load(F, fwd, v);
    prop_load(B, bwd, v);

    /* Local */
    for(int i=0; i<NLOC; i++) {
      _Complex double gF[NS*NC][NS*NC];          
      switch(i) {
	/* Scalar */
      case one:
	prop_1_G(gF, F);
	break;
	/* Pseudo-scalar */
      case g5:
	prop_g5_G(gF, F);
	break;
	/* Vector */
      case g0:
	prop_g0_G(gF, F);
	break;
      case gx:
	prop_gx_G(gF, F);
	break;
      case gy:
	prop_gy_G(gF, F);
	break;
      case gz:
	prop_gz_G(gF, F);
	break;
	/* Axial */
      case g5g0:
	prop_g5g0_G(gF, F);
	break;
      case g5gx:
	prop_g5gx_G(gF, F);
	break;
      case g5gy:
	prop_g5gy_G(gF, F);
	break;
      case g5gz:
	prop_g5gz_G(gF, F);
	break;
	/* Tensor */
      case g5si0x:
	prop_g5si0x_G(gF, F);
	break;
      case g5si0y:
	prop_g5si0y_G(gF, F);
	break;
      case g5si0z:
	prop_g5si0z_G(gF, F);
	break;
      case g5sixy:
	prop_g5sixy_G(gF, F);
	break;
      case g5sixz:
	prop_g5sixz_G(gF, F);
	break;
      case g5siyz:
	prop_g5siyz_G(gF, F);
	break;
      }
      prop_mul_gg(T, gF, B);
      corr.C[TIDX(v, i)] = prop_trace(T);
    }

    /* Noether */
    for(int i=NLOC; i<NLOC+NNOE; i++) {
      _Complex double gmFP[NS*NC][NS*NC];
      _Complex double gpFM[NS*NC][NS*NC];
      _Complex double gmF[NS*NC][NS*NC];
      _Complex double gpF[NS*NC][NS*NC];
      
      _Complex double FP[NS*NC][NS*NC];
      _Complex double FM[NS*NC][NS*NC];

      _Complex double BP[NS*NC][NS*NC];
      _Complex double BM[NS*NC][NS*NC];

      _Complex double A0[NS*NC][NS*NC];
      _Complex double A1[NS*NC][NS*NC];
      _Complex double A2[NS*NC][NS*NC];
      _Complex double A3[NS*NC][NS*NC];
      int mu = i - NLOC;
      int vp = nn[mu][v];
      int vm = nn[mu+ND][v];
      
      _Complex double U0[NC][NC];
      _Complex double Um[NC][NC];
      su3_load(U0, gf, v, mu);
      su3_load(Um, gf, vm, mu);
      
      prop_load(FP, fwd, vp);
      prop_load(FM, fwd, vm);
      prop_load(BP, bwd, vp);
      prop_load(BM, bwd, vm);
      switch(i) {
      case noe_g0:
	prop_1pg0_G(gpFM, FM);
	prop_1mg0_G(gmFP, FP);
	prop_1pg0_G(gpF, F);
	prop_1mg0_G(gmF, F);
	break;
      case noe_gx:
	prop_1pgx_G(gpFM, FM);
	prop_1mgx_G(gmFP, FP);
	prop_1pgx_G(gpF, F);
	prop_1mgx_G(gmF, F);
	break;
      case noe_gy:
	prop_1pgy_G(gpFM, FM);
	prop_1mgy_G(gmFP, FP);
	prop_1pgy_G(gpF, F);
	prop_1mgy_G(gmF, F);
	break;
      case noe_gz:
	prop_1pgz_G(gpFM, FM);
	prop_1mgz_G(gmFP, FP);
	prop_1pgz_G(gpF, F);
	prop_1mgz_G(gmF, F);
	break;
      }
      prop_mul_su3_D_G(T, gpF, U0);
      prop_mul_gg(A0, T, BP);

      prop_mul_su3_D_G(T, gpFM, Um);
      prop_mul_gg(A1, T, B);

      prop_mul_su3_U_G(T, gmFP, U0);
      prop_mul_gg(A2, T, B);
      
      prop_mul_su3_U_G(T, gmF, Um);
      prop_mul_gg(A3, T, BM);

      prop_add_gg(T, A0, A1);
      prop_meq_g(T, A2);
      prop_meq_g(T, A3);
      corr.C[TIDX(v, i)] = prop_trace(T)*0.25;
    }

    /* Compute gamma-less derivative */
    _Complex double C[ND][NS*NC][NS*NC];    
    for(int mu=0; mu<ND; mu++) {
      _Complex double FP[NS*NC][NS*NC];
      _Complex double FM[NS*NC][NS*NC];

      _Complex double BP[NS*NC][NS*NC];
      _Complex double BM[NS*NC][NS*NC];

      _Complex double A0[NS*NC][NS*NC];
      _Complex double A1[NS*NC][NS*NC];
      _Complex double A2[NS*NC][NS*NC];
      _Complex double A3[NS*NC][NS*NC];
      
      int vp = nn[mu][v];
      int vm = nn[mu+ND][v];
      
      _Complex double U0[NC][NC];
      _Complex double Um[NC][NC];
      su3_load(U0, gf, v, mu);
      su3_load(Um, gf, vm, mu);
      
      prop_load(FP, fwd, vp);
      prop_load(FM, fwd, vm);
      prop_load(BP, bwd, vp);
      prop_load(BM, bwd, vm);

      prop_mul_su3_U_G(T, FP, U0);
      prop_mul_gg(A0, T, B);      

      prop_mul_su3_D_G(T, FM, Um);
      prop_mul_gg(A1, T, B);      

      prop_mul_su3_D_G(T, F, U0);
      prop_mul_gg(A2, T, BP);      

      prop_mul_su3_U_G(T, F, Um);
      prop_mul_gg(A3, T, BM);      
      prop_sub_gg(C[mu], A0, A1);
      prop_meq_g(C[mu], A2);
      prop_peq_g(C[mu], A3);
    }

    /* Vector derivative */
    for(int i=NLOC+NNOE; i<NLOC+NNOE+NVDER; i++) {
      _Complex double gmuCnu[NC*NS][NC*NS];
      _Complex double gnuCmu[NC*NS][NC*NS];
      int mu, nu;
      switch(i) {
      case der_g0_0:
	mu = 0;
	nu = 0;
	prop_g0_G(gnuCmu, C[mu]);	
	prop_g0_G(gmuCnu, C[nu]);
	break;
      case der_gx_0:
	mu = 0;
	nu = 1;
	prop_gx_G(gnuCmu, C[mu]);	
	prop_g0_G(gmuCnu, C[nu]);
	break;
      case der_gy_0:
	mu = 0;
	nu = 2;
	prop_gy_G(gnuCmu, C[mu]);	
	prop_g0_G(gmuCnu, C[nu]);
	break;
      case der_gz_0:
	mu = 0;
	nu = 3;
	prop_gz_G(gnuCmu, C[mu]);	
	prop_g0_G(gmuCnu, C[nu]);
	break;
	//
      case der_gx_x:
	mu = 1;
	nu = 1;
	prop_gx_G(gnuCmu, C[mu]);	
	prop_gx_G(gmuCnu, C[nu]);
	break;
      case der_gy_x:
	mu = 1;
	nu = 2;
	prop_gy_G(gnuCmu, C[mu]);	
	prop_gx_G(gmuCnu, C[nu]);
	break;
      case der_gz_x:
	mu = 1;
	nu = 3;
	prop_gz_G(gnuCmu, C[mu]);	
	prop_gx_G(gmuCnu, C[nu]);
	break;
	//
      case der_gy_y:
	mu = 2;
	nu = 2;
	prop_gy_G(gnuCmu, C[mu]);	
	prop_gy_G(gmuCnu, C[nu]);
	break;
      case der_gz_y:
	mu = 2;
	nu = 3;
	prop_gz_G(gnuCmu, C[mu]);	
	prop_gy_G(gmuCnu, C[nu]);
	break;
	//
      case der_gz_z:
	mu = 3;
	nu = 3;
	prop_gz_G(gnuCmu, C[mu]);	
	prop_gz_G(gmuCnu, C[nu]);
	break;	  
      }
      
      prop_add_gg(T, gmuCnu, gnuCmu);
      corr.C[TIDX(v, i)] = prop_trace(T)*0.125;
    }

    /* Axial derivative, symmetric */
    for(int i=NLOC+NNOE+NVDER; i<NLOC+NNOE+NVDER+NADSY; i++) {
      _Complex double gmuCnu[NC*NS][NC*NS];
      _Complex double gnuCmu[NC*NS][NC*NS];
      int mu, nu;
      switch(i) {
      case der_g5g0_0_sym:
	mu = 0;
	nu = 0;
	prop_g5g0_G(gnuCmu, C[mu]);	
	prop_g5g0_G(gmuCnu, C[nu]);
	break;
      case der_g5gx_0_sym:
	mu = 0;
	nu = 1;
	prop_g5gx_G(gnuCmu, C[mu]);	
	prop_g5g0_G(gmuCnu, C[nu]);
	break;
      case der_g5gy_0_sym:
	mu = 0;
	nu = 2;
	prop_g5gy_G(gnuCmu, C[mu]);	
	prop_g5g0_G(gmuCnu, C[nu]);
	break;
      case der_g5gz_0_sym:
	mu = 0;
	nu = 3;
	prop_g5gz_G(gnuCmu, C[mu]);	
	prop_g5g0_G(gmuCnu, C[nu]);
	break;
	//
      case der_g5gx_x_sym:
	mu = 1;
	nu = 1;
	prop_g5gx_G(gnuCmu, C[mu]);	
	prop_g5gx_G(gmuCnu, C[nu]);
	break;
      case der_g5gy_x_sym:
	mu = 1;
	nu = 2;
	prop_g5gy_G(gnuCmu, C[mu]);	
	prop_g5gx_G(gmuCnu, C[nu]);
	break;
      case der_g5gz_x_sym:
	mu = 1;
	nu = 3;
	prop_g5gz_G(gnuCmu, C[mu]);	
	prop_g5gx_G(gmuCnu, C[nu]);
	break;
	//
      case der_g5gy_y_sym:
	mu = 2;
	nu = 2;
	prop_g5gy_G(gnuCmu, C[mu]);	
	prop_g5gy_G(gmuCnu, C[nu]);
	break;
      case der_g5gz_y_sym:
	mu = 2;
	nu = 3;
	prop_g5gz_G(gnuCmu, C[mu]);	
	prop_g5gy_G(gmuCnu, C[nu]);
	break;
	//
      case der_g5gz_z_sym:
	mu = 3;
	nu = 3;
	prop_g5gz_G(gnuCmu, C[mu]);	
	prop_g5gz_G(gmuCnu, C[nu]);
	break;	  
      }
      
      prop_add_gg(T, gmuCnu, gnuCmu);
      corr.C[TIDX(v, i)] = prop_trace(T)*0.125;
    }

    /* Axial derivative, anti-symmetric */
    for(int i=NLOC+NNOE+NVDER+NADSY; i<NLOC+NNOE+NVDER+NADSY+NADAS; i++) {
      _Complex double gmuCnu[NC*NS][NC*NS];
      _Complex double gnuCmu[NC*NS][NC*NS];
      int mu, nu;
      switch(i) {
      case der_g5gx_0_asy:
	mu = 0;
	nu = 1;
	prop_g5gx_G(gnuCmu, C[mu]);	
	prop_g5g0_G(gmuCnu, C[nu]);
	break;
      case der_g5gy_0_asy:
	mu = 0;
	nu = 2;
	prop_g5gy_G(gnuCmu, C[mu]);	
	prop_g5g0_G(gmuCnu, C[nu]);
	break;
      case der_g5gz_0_asy:
	mu = 0;
	nu = 3;
	prop_g5gz_G(gnuCmu, C[mu]);	
	prop_g5g0_G(gmuCnu, C[nu]);
	break;
	//
      case der_g5gy_x_asy:
	mu = 1;
	nu = 2;
	prop_g5gy_G(gnuCmu, C[mu]);	
	prop_g5gx_G(gmuCnu, C[nu]);
	break;
      case der_g5gz_x_asy:
	mu = 1;
	nu = 3;
	prop_g5gz_G(gnuCmu, C[mu]);	
	prop_g5gx_G(gmuCnu, C[nu]);
	break;
	//
      case der_g5gz_y_asy:
	mu = 2;
	nu = 3;
	prop_g5gz_G(gnuCmu, C[mu]);	
	prop_g5gy_G(gmuCnu, C[nu]);
	break;
      }
      
      prop_sub_gg(T, gnuCmu, gmuCnu);
      corr.C[TIDX(v, i)] = prop_trace(T)*0.125;
    }

    /* Tensor derivative, symmetric */
    for(int i=NLOC+NNOE+NVDER+NADSY+NADAS; i<NLOC+NNOE+NVDER+NADSY+NADAS+NTDER; i++) {
      _Complex double g5simunuCku[NC*NS][NC*NS];
      _Complex double g5simukuCnu[NC*NS][NC*NS];
      int mu, nu, ku;
      switch(i) {
      case der_g5si00_x:
	mu = 0;
	nu = 0;
	ku = 1;
	prop_zero(g5simunuCku);
	prop_g5si0x_G(g5simukuCnu, C[nu]);
	break;
      case der_g5si00_y:
	mu = 0;
	nu = 0;
	ku = 2;
	prop_zero(g5simunuCku);
	prop_g5si0y_G(g5simukuCnu, C[nu]);
	break;
      case der_g5si00_z:
	mu = 0;
	nu = 0;
	ku = 3;
	prop_zero(g5simunuCku);
	prop_g5si0z_G(g5simukuCnu, C[nu]);
	break;
      case der_g5si0x_x:
	mu = 0;
	nu = 1;
	ku = 1;
	prop_g5si0x_G(g5simunuCku, C[ku]);
	prop_g5si0x_G(g5simukuCnu, C[nu]);
	break;
      case der_g5si0x_y:
	mu = 0;
	nu = 1;
	ku = 2;
	prop_g5si0x_G(g5simunuCku, C[ku]);
	prop_g5si0y_G(g5simukuCnu, C[nu]);
	break;
      case der_g5si0x_z:
	mu = 0;
	nu = 1;
	ku = 3;
	prop_g5si0x_G(g5simunuCku, C[ku]);
	prop_g5si0z_G(g5simukuCnu, C[nu]);
	break;
      case der_g5si0y_y:
	mu = 0;
	nu = 2;
	ku = 2;
	prop_g5si0y_G(g5simunuCku, C[ku]);
	prop_g5si0y_G(g5simukuCnu, C[nu]);
	break;
      case der_g5si0y_z:
	mu = 0;
	nu = 2;
	ku = 3;
	prop_g5si0y_G(g5simunuCku, C[ku]);
	prop_g5si0z_G(g5simukuCnu, C[nu]);
	break;
      case der_g5si0z_z:
	mu = 0;
	nu = 3;
	ku = 3;
	prop_g5si0z_G(g5simunuCku, C[ku]);
	prop_g5si0z_G(g5simukuCnu, C[nu]);
	break;
	//
      case der_g5six0_0:
	mu = 1;
	nu = 0;
	ku = 0;
	prop_g5six0_G(g5simunuCku, C[ku]);
	prop_g5six0_G(g5simukuCnu, C[nu]);
	break;
      case der_g5six0_x:
	mu = 1;
	nu = 0;
	ku = 1;
	prop_g5six0_G(g5simunuCku, C[ku]);
	prop_zero(g5simukuCnu);
	break;
      case der_g5six0_y:
	mu = 1;
	nu = 0;
	ku = 2;
	prop_g5six0_G(g5simunuCku, C[ku]);
	prop_g5sixy_G(g5simukuCnu, C[nu]);
	break;
      case der_g5six0_z:
	mu = 1;
	nu = 0;
	ku = 3;
	prop_g5six0_G(g5simunuCku, C[ku]);
	prop_g5sixz_G(g5simukuCnu, C[nu]);
	break;
      case der_g5sixx_y:
	mu = 1;
	nu = 1;
	ku = 2;
	prop_zero(g5simunuCku);
	prop_g5sixy_G(g5simukuCnu, C[nu]);
	break;
      case der_g5sixx_z:
	mu = 1;
	nu = 1;
	ku = 3;
	prop_zero(g5simunuCku);
	prop_g5sixz_G(g5simukuCnu, C[nu]);
	break;
      case der_g5sixy_y:
	mu = 1;
	nu = 2;
	ku = 2;
	prop_g5sixy_G(g5simunuCku, C[ku]);
	prop_g5sixy_G(g5simukuCnu, C[nu]);
	break;
      case der_g5sixy_z:
	mu = 1;
	nu = 2;
	ku = 3;
	prop_g5sixy_G(g5simunuCku, C[ku]);
	prop_g5sixz_G(g5simukuCnu, C[nu]);
	break;
      case der_g5sixz_z:
	mu = 1;
	nu = 3;
	ku = 3;
	prop_g5sixz_G(g5simunuCku, C[ku]);
	prop_g5sixz_G(g5simukuCnu, C[nu]);
	break;
	//
      case der_g5siy0_0:
	mu = 2;
	nu = 0;
	ku = 0;
	prop_g5siy0_G(g5simunuCku, C[ku]);
	prop_g5siy0_G(g5simukuCnu, C[nu]);
	break;
      case der_g5siy0_x:
	mu = 2;
	nu = 0;
	ku = 1;
	prop_g5siy0_G(g5simunuCku, C[ku]);
	prop_g5siyx_G(g5simukuCnu, C[nu]);
	break;
      case der_g5siy0_y:
	mu = 2;
	nu = 0;
	ku = 2;
	prop_g5siy0_G(g5simunuCku, C[ku]);
	prop_zero(g5simukuCnu);
	break;
      case der_g5siy0_z:
	mu = 2;
	nu = 0;
	ku = 3;
	prop_g5siy0_G(g5simunuCku, C[ku]);
	prop_g5siyz_G(g5simukuCnu, C[nu]);
	break;
      case der_g5siyx_x:
	mu = 2;
	nu = 1;
	ku = 1;
	prop_g5siyx_G(g5simunuCku, C[ku]);
	prop_g5siyx_G(g5simukuCnu, C[nu]);
	break;
      case der_g5siyx_y:
	mu = 2;
	nu = 1;
	ku = 2;
	prop_g5siyx_G(g5simunuCku, C[ku]);
	prop_zero(g5simukuCnu);
	break;
      case der_g5siyx_z:
	mu = 2;
	nu = 1;
	ku = 3;
	prop_g5siyx_G(g5simunuCku, C[ku]);
	prop_g5siyz_G(g5simukuCnu, C[nu]);
	break;
      case der_g5siyy_z:
	mu = 2;
	nu = 2;
	ku = 3;
	prop_zero(g5simunuCku);
	prop_g5siyz_G(g5simukuCnu, C[nu]);
	break;
      case der_g5siyz_z:
	mu = 2;
	nu = 3;
	ku = 3;
	prop_g5siyz_G(g5simunuCku, C[ku]);
	prop_g5siyz_G(g5simukuCnu, C[nu]);
	break;
	//
      case der_g5siz0_0:
	mu = 3;
	nu = 0;
	ku = 0;
	prop_g5siz0_G(g5simunuCku, C[ku]);
	prop_g5siz0_G(g5simukuCnu, C[nu]);
	break;
      case der_g5siz0_x:
	mu = 3;
	nu = 0;
	ku = 1;
	prop_g5siz0_G(g5simunuCku, C[ku]);
	prop_g5sizx_G(g5simukuCnu, C[nu]);
	break;
      case der_g5siz0_y:
	mu = 3;
	nu = 0;
	ku = 2;
	prop_g5siz0_G(g5simunuCku, C[ku]);
	prop_g5sizy_G(g5simukuCnu, C[nu]);
	break;
      case der_g5siz0_z:
	mu = 3;
	nu = 0;
	ku = 3;
	prop_g5siz0_G(g5simunuCku, C[ku]);
	prop_zero(g5simukuCnu);
	break;
      case der_g5sizx_x:
	mu = 3;
	nu = 1;
	ku = 1;
	prop_g5sizx_G(g5simunuCku, C[ku]);
	prop_g5sizx_G(g5simukuCnu, C[nu]);
	break;
      case der_g5sizx_y:
	mu = 3;
	nu = 1;
	ku = 2;
	prop_g5sizx_G(g5simunuCku, C[ku]);
	prop_g5sizy_G(g5simukuCnu, C[nu]);
	break;
      case der_g5sizx_z:
	mu = 3;
	nu = 1;
	ku = 3;
	prop_g5sizx_G(g5simunuCku, C[ku]);
	prop_zero(g5simukuCnu);
	break;
      case der_g5sizy_y:
	mu = 3;
	nu = 2;
	ku = 2;
	prop_g5sizy_G(g5simunuCku, C[ku]);
	prop_g5sizy_G(g5simukuCnu, C[nu]);
	break;
      case der_g5sizy_z:
	mu = 3;
	nu = 2;
	ku = 3;
	prop_g5sizy_G(g5simunuCku, C[ku]);
	prop_zero(g5simukuCnu);
	break;
      }
      
      prop_add_gg(T, g5simunuCku, g5simukuCnu);
      corr.C[TIDX(v, i)] = prop_trace(T)*0.125;
    }  
  }
    
  for(int i=0; i<ND; i++)
    corr.origin[i] = source_coords[i];
  int T = lat->dims[0];
  corr.cutoff[0] = (source_coords[0] + thrp_sink.dt) % T;
  corr.mom_list = NULL;
  return corr;
}

void
qhg_write_nn_thrp(char fname[], qhg_correlator corr)
{
  qhg_lattice *lat = corr.lat;
  int proc_id = lat->comms->proc_id;
  int am_io_proc = proc_id == 0 ? 1 : 0;
  int site_size = corr.site_size;
  int nm = corr.mom_list->n_mom_vecs;
  int lvol = lat->lvol;
  int (*mv)[3] = corr.mom_list->mom_vecs;
  int *pd = lat->comms->proc_dims;
  int np = lat->comms->nprocs;
  int *ld = lat->ldims;
  int *d = lat->dims;
  int nlines[np];
  for(int i=0; i<np; i++)
    nlines[i] = 0;

  int T = lat->dims[0];
  int L = lat->dims[1];
  int ts = corr.cutoff[0];
  int t0 = corr.origin[0];
  int dt = (T+ts-t0) % T;
  
  lines lines_loc = lines_new(dt*nm*NCHAN);
  for(int ta=t0; ta<=dt+t0; ta++) {
    int t = ta % T;
    for(int n=0; n<nm; n++) {
      /* Momentum vector */
      int *k = mv[n];
      /* Global coordinates of momentum vector */
      int x[] = {t, (k[0]+d[1]) % d[1], (k[1]+d[2]) % d[2], (k[2]+d[3]) % d[3]};
      /* Coordinates of process on which x[] resides */
      int pc[] = {x[0]/ld[0], x[1]/ld[1], x[2]/ld[2], x[3]/ld[3]};
      /* Lexico index of pc[] */
      int ip = IDX(pc, pd);
      /* Local coordinates of x[] within process ip */
      int lc[] = {x[0]%ld[0], x[1]%ld[1], x[2]%ld[2], x[3]%ld[3]};
      /* Lexico index of lc[] */
      int ix = IDX(lc,ld);
      _Complex double *c = corr.C;
      line li[NCHAN];
      /* Shift to time relative to source */
      for(int ich=0; ich<NCHAN; ich++) {
	li[ich].n = ich + NCHAN*(n + nm*(ta-t0));
	sprintf(li[ich].c, "%4d %+d %+d %+d  %+e %+e %s\n",
		ta-t0, k[0], k[1], k[2],
		creal(c[TIDX(ix, ich)]), cimag(c[TIDX(ix, ich)]),
		chan_tags[ich]);
      }
      nlines[ip] += NCHAN;
      if(proc_id == ip)
	lines_loc = lines_append(lines_loc, li, NCHAN);
      
    }
  }
  
  int nl = 0;
  int displs[np];
  for(int i=0; i<np; i++)
    displs[i] = 0;
  
  for(int i=0; i<np; i++) {
    nl += nlines[i];
    nlines[i] *= sizeof(line);
    for(int j=0; j<i; j++)
      displs[i] += nlines[j];
  }
  
  lines lines_glob = lines_new(nl);
  MPI_Gatherv(lines_loc.l, nlines[proc_id], MPI_BYTE,
	      lines_glob.l, nlines, displs, MPI_BYTE, 0, lat->comms->comm);
  lines_glob.cur = lines_glob.len;
  lines_glob = lines_sorted(lines_glob, NCHAN);

  if(am_io_proc) {
    FILE *fp = qhg_fopen(fname, "w");
    for(int i=0; i<nl; i++)
      fprintf(fp, "%s", lines_glob.l[i].c);
    fclose(fp);
  }
  
  lines_del(lines_loc);
  lines_del(lines_glob);  
  return;
}
