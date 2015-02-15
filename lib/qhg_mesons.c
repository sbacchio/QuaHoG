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

#define NGAMMAS 10			/* 1,gx,...,gt,g5,gxg5,...,gtg5 */
#define NFLAVS 2			/* up/down */
#define NCHANNELS ((NGAMMAS)*(NFLAVS))
#define VGF(v, g, f) ((v)*NCHANNELS + (g)*NFLAVS + (f))
#define VSC(v, cs) ((v)*NC*NS + (cs))

static char gamma_tags[NGAMMAS][256] = {"   =1=\0",
					"  =g5=\0",
					"  =gx=\0",
					"  =gy=\0",
					"  =gz=\0",
					"  =gt=\0",
					"=g5gx=\0",
					"=g5gy=\0",
					"=g5gz=\0",
					"=g5gt=\0"};

static char flav_tags[NFLAVS][256] = {"up\0",
				      "dn\0"};

qhg_correlator
qhg_mesons(qhg_spinor_field sp_u[NS*NC], qhg_spinor_field sp_d[NS*NC], int source_coords[ND])
{
  qhg_lattice *lat = sp_u[0].lat;
  qhg_correlator corr = qhg_correlator_init(NCHANNELS, lat);
  int lvol = lat->lvol;
  for(int v=0; v<lvol; v++) {
    _Complex double U[NS*NC][NS*NC];
    _Complex double D[NS*NC][NS*NC];    
    _Complex double C[NS*NC][NS*NC];    
    _Complex double W[NS*NC][NS*NC];
    _Complex double V[NS*NC][NS*NC];        
    for(int cs0=0; cs0<NS*NC; cs0++)
      for(int cs1=0; cs1<NS*NC; cs1++) {
	U[cs0][cs1] = sp_u[cs1].field[VSC(v, cs0)];
	D[cs0][cs1] = sp_d[cs1].field[VSC(v, cs0)];	
      }
    for(int igamma=0; igamma<NGAMMAS; igamma++) {
      _Complex double (*P[2])[NC*NS] = {U, D};
      for(int iflav=0; iflav<NFLAVS; iflav++) {
	/*
	 * Sign flips are for consistency with libqcd
	 */
	switch(igamma) {
	case 0: /* 1 */
	  prop_mul_gd(C, P[iflav], P[iflav]);
	  corr.C[VGF(v, igamma, iflav)] = prop_trace(C);
	  break;
	case 1: /* g5 */
	  prop_g5_G(W, P[iflav]);
	  prop_G_g5(V, W);	  
	  prop_mul_gd(C, V, P[iflav]);
	  corr.C[VGF(v, igamma, iflav)] = -prop_trace(C);
	  break;
	case 2: /* gx */
	  prop_gx_G(W, P[iflav]);
	  prop_G_gx(V, W);	  
	  prop_mul_gd(C, V, P[iflav]);
	  corr.C[VGF(v, igamma, iflav)] = -prop_trace(C);
	  break;
	case 3: /* gy */
	  prop_gy_G(W, P[iflav]);
	  prop_G_gy(V, W);	  
	  prop_mul_gd(C, V, P[iflav]);
	  corr.C[VGF(v, igamma, iflav)] = -prop_trace(C);
	  break;
	case 4: /* gz */
	  prop_gz_G(W, P[iflav]);
	  prop_G_gz(V, W);	  
	  prop_mul_gd(C, V, P[iflav]);
	  corr.C[VGF(v, igamma, iflav)] = -prop_trace(C);
	  break;
	case 5: /* gt */
	  prop_gt_G(W, P[iflav]);
	  prop_G_gt(V, W);	  
	  prop_mul_gd(C, V, P[iflav]);
	  corr.C[VGF(v, igamma, iflav)] = -prop_trace(C);
	  break;
	case 6: /* g5gx */
	  prop_g5gx_G(W, P[iflav]);
	  prop_G_g5gx(V, W);	  
	  prop_mul_gd(C, V, P[iflav]);
	  corr.C[VGF(v, igamma, iflav)] = -prop_trace(C);
	  break;
	case 7: /* g5gy */
	  prop_g5gy_G(W, P[iflav]);
	  prop_G_g5gy(V, W);	  
	  prop_mul_gd(C, V, P[iflav]);
	  corr.C[VGF(v, igamma, iflav)] = -prop_trace(C);
	  break;
	case 8: /* g5gz */
	  prop_g5gz_G(W, P[iflav]);
	  prop_G_g5gz(V, W);	  
	  prop_mul_gd(C, V, P[iflav]);
	  corr.C[VGF(v, igamma, iflav)] = -prop_trace(C);
	  break;
	case 9: /* g5gt */
	  prop_g5gt_G(W, P[iflav]);
	  prop_G_g5gt(V, W);	  
	  prop_mul_gd(C, V, P[iflav]);
	  corr.C[VGF(v, igamma, iflav)] = prop_trace(C);
	  break;	  
	}
      }
    }
  }
  
  for(int i=0; i<ND; i++)
    corr.origin[i] = source_coords[i];
  corr.mom_list = NULL;
  return corr;
}

void
qhg_write_mesons(char fname[], qhg_correlator corr)
{
  qhg_lattice *lat = corr.lat;
  int am_io_proc = lat->comms->proc_id == 0 ? 1 : 0;
  int site_size = corr.site_size;
  int nm = corr.mom_list->n_mom_vecs;
  int lvol = lat->lvol;
  int (*mv)[3] = corr.mom_list->mom_vecs;
  int *pd = lat->comms->proc_dims;
  int *ld = lat->ldims;
  int *d = lat->dims;

  /* Gather the correlator into the 0th process */  
  size_t nelems = site_size*lat->vol;
  size_t nelems_loc = site_size*lat->lvol;  
  _Complex double *data = NULL;
  if(am_io_proc) {
    data = qhg_alloc(sizeof(_Complex double)*nelems);  
    memset(data, '\0', sizeof(_Complex double)*nelems);
  }
  MPI_Gather(corr.C, 2*nelems_loc, MPI_DOUBLE, data, 2*nelems_loc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  /* Loop over momentum vectors, calculate which process the momentum
     resides on, and print */
  if(am_io_proc) {
    FILE *fp = fopen(fname, "w");
    for(int tt=0; tt<d[0]; tt++) {    
      int t = (d[0] + tt + corr.origin[0]) % d[0];
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
	_Complex double *c = &data[ip*nelems_loc];
	for(int ig=0; ig<NGAMMAS; ig++)
	  for(int ifl=0; ifl<NFLAVS; ifl++) {	  
	    fprintf(fp, "%4d %+d %+d %+d %+e %+e \t%s\t%s\n",
		    tt, k[0], k[1], k[2],
		    creal(c[VGF(ix, ig, ifl)]), cimag(c[VGF(ix, ig, ifl)]),
		    gamma_tags[ig],flav_tags[ifl]);
	  }
      }
    }
    fclose(fp);
  }

  if(am_io_proc)
    free(data);

  MPI_Barrier(lat->comms->comm);
  return;
}
