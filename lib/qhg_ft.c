#include <mpi.h>
#include <complex.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <qhg_defs.h>
#include <qhg_types.h>
#include <qhg_idx.h>
#include <qhg_alloc.h>
#include <qhg_mom_list.h>
#include <qhg_correlator.h>

static int
get_ft_sign(char dir[], qhg_lattice *lat)
{
  int am_io_proc = lat->comms->proc_id == 0 ? 1 : 0;
  int s = 0;
  if(strcmp(dir, "fwd") == 0) {
    s = -1;
  } else if(strcmp(dir, "bwd") == 0) {
    s = +1;
  } else {
    if(am_io_proc)
      fprintf(stderr, "Bad direction argument in %s; should be one of %s or %s\n",
	      __func__, "fwd", "bwd");
    MPI_Abort(lat->comms->comm, 3);
  }
  return s;
}

qhg_correlator
qhg_ft(qhg_correlator corr_x, qhg_mom_list *mom_list, char direction[])
{
  qhg_lattice *lat = corr_x.lat;
  int sign = get_ft_sign(direction, lat);

  int *ldims = lat->ldims;
  int lt = ldims[0];
  int lx = ldims[1];
  int ly = ldims[2];
  int lz = ldims[3];  

  int Lx = lat->dims[1];
  int Ly = lat->dims[2];
  int Lz = lat->dims[3];  

  int x0 = corr_x.origin[1];
  int y0 = corr_x.origin[2];
  int z0 = corr_x.origin[3];  

  int *pdims = lat->comms->proc_dims;
  int lvol = lat->lvol;
  int site_size = corr_x.site_size;  

  qhg_correlator corr_p = qhg_correlator_init(site_size, lat);
  memset(corr_p.C, '\0', site_size*lvol*sizeof(_Complex double));
  
  /* Split the communicator so that all processes of the same
     time-slices have a separate communicator. This way we can do a
     reduction over the new communicator. */
  MPI_Comm tcomm;
  int *pcoords = lat->comms->proc_coords;
  int proc_id = lat->comms->proc_id;
  MPI_Comm_split(lat->comms->comm, pcoords[0], proc_id, &tcomm);  

  /* The rank within the tcomm communicator */
  int s_proc_id;
  MPI_Comm_rank(tcomm, &s_proc_id);  

  _Complex double *ftacc = qhg_alloc(lt*site_size*sizeof(_Complex double));
  _Complex double *recv = qhg_alloc(lt*site_size*sizeof(_Complex double));
  int nmom = mom_list->n_mom_vecs;
  int (*moms)[3] = mom_list->mom_vecs;
  for(int m=0; m<nmom; m++) {
    int *k = moms[m];
    memset(ftacc, '\0', lt*site_size*sizeof(_Complex double));
    for(int x=0; x<lx; x++)
      for(int y=0; y<ly; y++)
	for(int z=0; z<lz; z++) {
	  int xx = pcoords[1]*lx + x;
	  int yy = pcoords[2]*ly + y;
	  int zz = pcoords[3]*lz + z;
	  double ph =
	    (double)(Lx + xx - x0)*k[0]/Lx +
	    (double)(Ly + yy - y0)*k[1]/Ly +
	    (double)(Lz + zz - z0)*k[2]/Lz;
	  ph = ph*2*M_PI;
	  _Complex double ex = cos(ph) + sign * _Complex_I * sin(ph);
	  for(int t=0; t<lt; t++) {
	    int co[] = {t,x,y,z};
	    int v = IDX(co, ldims);
	    for(int s=0; s<site_size; s++)
	      ftacc[t*site_size + s] += ex*corr_x.C[v*site_size + s];
	  }
	}
    int kx = (Lx+k[0]) % Lx;
    int ky = (Ly+k[1]) % Ly;
    int kz = (Lz+k[2]) % Lz;

    /* Process coords of kx,ky,kz */
    int pco[] = {kx/lx, ky/ly, kz/lz};
    int x_proc_id = pco[2] + pdims[3]*(pco[1] + pdims[2]*pco[0]);

    /* Reduce on process which holds kx, ky, kz */
    MPI_Reduce(ftacc, recv, 2*lt*site_size, MPI_DOUBLE, MPI_SUM, x_proc_id, tcomm);

    /* Local coords of kx, ky, kz */
    int kco[] = {kx % lx, ky % ly, kz % lz};
    int vk = kco[2] + lz*(kco[1] + ly*kco[0]);
    int lv3 = lat->lv3;
    if(s_proc_id == x_proc_id)
      for(int t=0; t<lt; t++)
	for(int s=0; s<site_size; s++)
	  corr_p.C[(t*lv3 + vk)*site_size + s] = recv[t*site_size + s];
  }

  free(ftacc);
  free(recv);

  corr_p.mom_list = mom_list;
  corr_p.origin[0] = corr_x.origin[0];
  corr_p.origin[1] = 0;
  corr_p.origin[2] = 0;
  corr_p.origin[3] = 0;  
  return corr_p;
}
