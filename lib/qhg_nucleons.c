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

#define NFLAV 2
#define NCOMP (4)
#define NCHAN (4)
#define SITE_SIZE (NFLAV*NCOMP*NCOMP*NCHAN)
#define NIDX(v, f, chi, si0, si1) (si1 + NCOMP*(si0 + NCOMP*(chi + NCHAN*(f+NFLAV*v))))
#define VSC(v, sc) (sc + NC*NS*v)

static char flav_tags[NFLAV][256] = {"ppm\0",
				     "pmm\0"};

static char chan_tags[NCHAN][256] = {"1-1\0",
				     "1-2\0",
				     "2-1\0",
				     "2-2\0"};


qhg_correlator
qhg_nucleons(qhg_spinor_field sp_u[NS*NC], qhg_spinor_field sp_d[NS*NC], int source_coords[ND])
{
  qhg_lattice *lat = sp_u[0].lat;
  qhg_correlator corr = qhg_correlator_init(SITE_SIZE, lat);
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
    _Complex double (*P[2])[NC*NS] = {U, D};
    for(int iflav=0; iflav<NFLAV; iflav++) {
      _Complex double (*y)[NC*NS] = P[(iflav+0)%2];
      _Complex double (*x)[NC*NS] = P[(iflav+1)%2];
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

      for(int s0=0; s0<NS; s0++)
	for(int s1=0; s1<NS; s1++){
	  A[s0][s1] = 0.;
	  B[s0][s1] = 0.;
	}
	  
      for(int c1=0; c1<NC; c1++)
	for(int c0=0; c0<NC; c0++)
	  for(int s0=0; s0<NS; s0++)
	    for(int s1=0; s1<NS; s1++)
	      for(int s2=0; s2<NS; s2++) {
		A[s0][s1] += y[CS(s0,c0)][CS(s1,c1)]*z[CS(s2,c0)][CS(s2,c1)];
		B[s0][s1] += y[CS(s0,c0)][CS(s2,c1)]*z[CS(s2,c0)][CS(s1,c1)];		
	      }
      
      for(int s0=0; s0<NS; s0++)
	for(int s1=0; s1<NS; s1++) {
	  corr.C[NIDX(v, iflav, 0, s0, s1)] = A[s0][s1] + B[s0][s1];
	}

      /*
       * \chi_1 - to - \chi_2
       */
      _Complex double yg5[NS*NC][NS*NC];
      prop_G_g5(yg5, y);
      _Complex double Cg5xC[NS*NC][NS*NC];
      prop_Cg5_G(aux, x);
      prop_G_C(Cg5xC, aux);
      prop_contract_13(z, Cg5xC, y);

      for(int s0=0; s0<NS; s0++)
	for(int s1=0; s1<NS; s1++){
	  A[s0][s1] = 0.;
	  B[s0][s1] = 0.;
	}
	  
      for(int c1=0; c1<NC; c1++)
	for(int c0=0; c0<NC; c0++)
	  for(int s0=0; s0<NS; s0++)
	    for(int s1=0; s1<NS; s1++)
	      for(int s2=0; s2<NS; s2++) {
		A[s0][s1] += yg5[CS(s0,c0)][CS(s1,c1)]*z[CS(s2,c0)][CS(s2,c1)];
		B[s0][s1] += yg5[CS(s2,c0)][CS(s1,c1)]*z[CS(s2,c0)][CS(s0,c1)];		
	      }
      
      for(int s0=0; s0<NS; s0++)
	for(int s1=0; s1<NS; s1++) {
	  corr.C[NIDX(v, iflav, 1, s0, s1)] = A[s0][s1] + B[s0][s1];
	}

      /*
       * \chi_2 - to - \chi_1
       */
      _Complex double g5y[NS*NC][NS*NC];
      prop_g5_G(g5y, y);
      _Complex double CxCg5[NS*NC][NS*NC];
      prop_G_Cg5(aux, x);
      prop_C_G(CxCg5, aux);
      prop_contract_02(z, CxCg5, y);

      for(int s0=0; s0<NS; s0++)
	for(int s1=0; s1<NS; s1++){
	  A[s0][s1] = 0.;
	  B[s0][s1] = 0.;
	}
	  
      for(int c1=0; c1<NC; c1++)
	for(int c0=0; c0<NC; c0++)
	  for(int s0=0; s0<NS; s0++)
	    for(int s1=0; s1<NS; s1++)
	      for(int s2=0; s2<NS; s2++) {
		A[s0][s1] += g5y[CS(s0,c0)][CS(s1,c1)]*z[CS(s2,c0)][CS(s2,c1)];
		B[s0][s1] += g5y[CS(s0,c0)][CS(s2,c1)]*z[CS(s2,c0)][CS(s1,c1)];		
	      }
      
      for(int s0=0; s0<NS; s0++)
	for(int s1=0; s1<NS; s1++) {
	  corr.C[NIDX(v, iflav, 2, s0, s1)] = A[s0][s1] + B[s0][s1];
	}


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

      for(int s0=0; s0<NS; s0++)
	for(int s1=0; s1<NS; s1++){
	  A[s0][s1] = 0.;
	  B[s0][s1] = 0.;
	}
	  
      for(int c1=0; c1<NC; c1++)
	for(int c0=0; c0<NC; c0++)
	  for(int s0=0; s0<NS; s0++)
	    for(int s1=0; s1<NS; s1++)
	      for(int s2=0; s2<NS; s2++) {
		A[s0][s1] += g5yg5[CS(s0,c0)][CS(s1,c1)]*z[CS(s2,c0)][CS(s2,c1)];
		B[s0][s1] += g5y[CS(s0,c0)][CS(s2,c1)]*w[CS(s2,c0)][CS(s1,c1)];		
	      }
      
      for(int s0=0; s0<NS; s0++)
	for(int s1=0; s1<NS; s1++) {
	  corr.C[NIDX(v, iflav, 3, s0, s1)] = A[s0][s1] + B[s0][s1];
	}
    }
  }

  
  for(int i=0; i<ND; i++)
    corr.origin[i] = source_coords[i];
  corr.mom_list = NULL;
  return corr;
}

void
qhg_write_nucleons(char fname[], qhg_correlator corr)
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
	for(int ifl=0; ifl<NFLAV; ifl++)
	  for(int ich=0; ich<NCHAN; ich++)
	    for(int si0=0; si0<NCOMP; si0++) {
	      fprintf(fp, "%4d %+d %+d %+d  %+e %+e  %+e %+e  %+e %+e  %+e %+e  %s %s\n",
		      tt, k[0], k[1], k[2],
		      creal(c[NIDX(ix, ifl, ich, si0, 0)]), cimag(c[NIDX(ix, ifl, ich, si0, 0)]),
		      creal(c[NIDX(ix, ifl, ich, si0, 1)]), cimag(c[NIDX(ix, ifl, ich, si0, 1)]),
		      creal(c[NIDX(ix, ifl, ich, si0, 2)]), cimag(c[NIDX(ix, ifl, ich, si0, 2)]),
		      creal(c[NIDX(ix, ifl, ich, si0, 3)]), cimag(c[NIDX(ix, ifl, ich, si0, 3)]),
		      flav_tags[ifl],
		      chan_tags[ich]);
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
