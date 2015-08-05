#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>

#include <complex.h>
#include <qhg_io_utils.h>
#include <qhg_idx.h>
#include <qhg_defs.h>
#include <qhg_types.h>
#include <qhg_alloc.h>
#include <qhg_meson_udsc_defs.h>
#include <qhg_correlator.h>
#include <qhg_correlator_shift.h>
#include <lines_types.h>
#include <lines_utils.h>

void
qhg_write_mom_mesons_udsc(char fname[], qhg_correlator corr)
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
  lines lines_loc = lines_new(d[0]*nm*NGAMMAS*NFLAV*NFLAV);
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
      line li[NGAMMAS*NFLAV*NFLAV];
      /* Shift to time relative to source */
      int tt = (d[0] + t - corr.origin[0]) % d[0];
      for(int ig=0; ig<NGAMMAS; ig++)
        for(int ifl0=0; ifl0<NFLAV; ifl0++)
	  for(int ifl1=0; ifl1<NFLAV; ifl1++) {
	    int j = ifl1 + ifl0*NFLAV + ig*NFLAV*NFLAV;
	    li[j].n = j + NGAMMAS*NFLAV*NFLAV*(n + nm*tt);
	    sprintf(li[j].c, "%4d %+d %+d %+d %+e %+e \t%s\t%s%s\n",
		    tt, k[0], k[1], k[2],
		    creal(c[VGF(ix, ig, ifl0, ifl1)]), cimag(c[VGF(ix, ig, ifl0, ifl1)]),
		    gamma_tags[ig],flav_tags[ifl0],flav_tags[ifl1]);
	  }
      nlines[ip] += NGAMMAS*NFLAV*NFLAV;
      if(proc_id == ip)
        lines_loc = lines_append(lines_loc, li, NGAMMAS*NFLAV*NFLAV);
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
  lines_glob = lines_sorted(lines_glob, NGAMMAS*NFLAV*NFLAV);

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
