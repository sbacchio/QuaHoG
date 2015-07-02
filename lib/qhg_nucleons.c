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
#include <qhg_io_utils.h>
#include <lines_types.h>
#include <lines_utils.h>

#define NFLAV 2
#define NCOMP (4)
#define NCHAN (4)
#define SITE_SIZE (NFLAV*NCOMP*NCOMP*NCHAN)
#define NIDX(v, f, chi, si0, si1) (si1 + NCOMP*(si0 + NCOMP*(chi + NCHAN*(f+NFLAV*v))))

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
  for(int i=0; i<ND; i++)
    corr.origin[i] = source_coords[i];
  int tsrc = corr.origin[0];  
  int lvol = lat->lvol;
  int lv3 = lat->lv3;
  int lt = lat->ldims[0];
  int t0 = lat->ldims[0]*lat->comms->proc_coords[0];
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
  
  lines lines_loc = lines_new(d[0]*nm*NFLAV*NCHAN*NCOMP*NCOMP);
  for(int t=0; t<d[0]; t++) {    
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
      line li[NFLAV*NCHAN*NCOMP];
      /* Shift to time relative to source */
      int tt = (d[0] + t - corr.origin[0]) % d[0];      
      for(int ifl=0; ifl<NFLAV; ifl++)
	for(int ich=0; ich<NCHAN; ich++)
	  for(int si0=0; si0<NCOMP; si0++) {
	    int j = si0 + NCOMP*(ich + NCHAN*ifl);
	    li[j].n = j + NFLAV*NCHAN*NCOMP*(n + nm*tt);
	    sprintf(li[j].c, "%4d %+d %+d %+d  %+e %+e  %+e %+e  %+e %+e  %+e %+e  %s %s\n",
		    tt, k[0], k[1], k[2],
		    creal(c[NIDX(ix, ifl, ich, si0, 0)]), cimag(c[NIDX(ix, ifl, ich, si0, 0)]),
		    creal(c[NIDX(ix, ifl, ich, si0, 1)]), cimag(c[NIDX(ix, ifl, ich, si0, 1)]),
		    creal(c[NIDX(ix, ifl, ich, si0, 2)]), cimag(c[NIDX(ix, ifl, ich, si0, 2)]),
		    creal(c[NIDX(ix, ifl, ich, si0, 3)]), cimag(c[NIDX(ix, ifl, ich, si0, 3)]),
		    flav_tags[ifl],
		    chan_tags[ich]);
	  }
      nlines[ip] += NFLAV*NCHAN*NCOMP;
      if(proc_id == ip)
	lines_loc = lines_append(lines_loc, li, NFLAV*NCHAN*NCOMP);
      
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
  lines_glob = lines_sorted(lines_glob, NFLAV*NCHAN*NCOMP);

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
