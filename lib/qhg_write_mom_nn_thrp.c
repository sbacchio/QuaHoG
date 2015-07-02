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
#include <qhg_nn_thrp_defs.h>
#include <qhg_correlator.h>
#include <qhg_correlator_shift.h>
#include <lines_types.h>
#include <lines_utils.h>

void
qhg_write_mom_nn_thrp(char fname[], qhg_thrp_correlator corr_thrp)
{
  qhg_correlator corr = corr_thrp.corr;
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
  int t0 = corr.origin[0];
  int dt = corr_thrp.dt;
  int ts = corr.cutoff[0];

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
