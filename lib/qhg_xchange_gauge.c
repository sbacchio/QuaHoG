#include <stdlib.h>
#include <string.h>

#include <mpi.h>
#include <complex.h>
#include <tmLQCD.h>
#include <qhg_defs.h>
#include <qhg_idx.h>
#include <qhg_types.h>
#include <qhg_alloc.h>
#include <qhg_xchange_utils.h>

static void
get_boundary(_Complex double *bnd, int dir, _Complex double *field, qhg_lattice *lat)
{
  int sign = dir / ND == 0 ? +1 : -1;
  dir = dir % ND;
  int bvol = lat->bvol[dir];
  int *ldims = lat->ldims;
  int *procs = lat->comms->proc_dims;
  if(procs[dir] > 1) {
    size_t site_size = NC*NC*ND*sizeof(_Complex double);
    MPI_Datatype bnd_collect = get_slice(ND, ldims, dir, sign == +1 ? 0 : ldims[dir]-1, site_size);
    MPI_Comm comm = lat->comms->comm;
    int proc_id = lat->comms->proc_id;
    int fproc = lat->comms->neigh_proc[dir];
    int bproc = lat->comms->neigh_proc[dir+ND];
    if(sign > 0) {
      MPI_Sendrecv(field, 1, bnd_collect, bproc, bproc,
		   bnd, bvol*site_size, MPI_BYTE, fproc, proc_id, comm,
		   MPI_STATUS_IGNORE);
    } else {
      MPI_Sendrecv(field, 1, bnd_collect, fproc, fproc,
		   bnd, bvol*site_size, MPI_BYTE, bproc, proc_id, comm,
		   MPI_STATUS_IGNORE);
    }
    MPI_Type_free(&bnd_collect);
  }
  return;
}

static void
get_boundary_edge(qhg_gauge_field gf, int dir)
{
  qhg_lattice *lat = gf.lat;
  int d0 = dir % ND;
  int bvol = lat->bvol[d0];
  int *ldims = lat->ldims;
  int *procs = lat->comms->proc_dims;
  int *evol = lat->evol[d0];
  if(procs[d0] > 1) {
    size_t site_size = NC*NC*ND*sizeof(_Complex double);
    int n_edge = 0;
    for(int i=ND-1; i>d0; i--)
      if(procs[i] > 1)
	n_edge++;

    int blen_arr[1+2*n_edge];
    MPI_Datatype type_arr[1+2*n_edge];
    MPI_Aint displ_arr[1+2*n_edge];

    blen_arr[0] = 1;
    type_arr[0] = get_slice(ND, ldims, d0, (dir < ND) ? 0 : (ldims[d0]-1), site_size);
    displ_arr[0] = 0;

    int k = 1;
    for(int i=ND-1; i>d0; i--)
      if(procs[i] > 1) {
	int bdims[ND];
	for(int j=0; j<ND; j++)
	  bdims[j] = ldims[j];

	bdims[i] = 1;

	blen_arr[k] = 1;
	blen_arr[k+1] = 1;	
	type_arr[k] = get_slice(ND, bdims, d0, (dir < ND) ? 0 : (ldims[d0]-1), site_size);
	type_arr[k+1] = get_slice(ND, bdims, d0, (dir < ND) ? 0 : (ldims[d0]-1), site_size);
	displ_arr[k] = (size_t)gf.bnd[i]-(size_t)gf.field;
	displ_arr[k+1] = (size_t)gf.bnd[i+ND]-(size_t)gf.field;
	k+=2;
      }

    MPI_Datatype bnd_edge_collect;
    MPI_Type_create_struct(1+2*n_edge, blen_arr, displ_arr, type_arr, &bnd_edge_collect);
    MPI_Type_commit(&bnd_edge_collect);

    for(int i=0; i<1+2*n_edge; i++)
      MPI_Type_free(&type_arr[i]);

    blen_arr[0] = bvol*site_size;
    type_arr[0] = MPI_BYTE;
    displ_arr[0] = 0;

    k = 1;
    for(int i=ND-1; i>d0; i--)
      if(procs[i] > 1) {
	blen_arr[k] = evol[i]*site_size;
	blen_arr[k+1] = evol[i]*site_size;	
	type_arr[k] = MPI_BYTE;
	type_arr[k+1] = MPI_BYTE;	
	displ_arr[k] = (size_t)gf.edge[dir][i] - (size_t)gf.bnd[dir];
	displ_arr[k+1] = (size_t)gf.edge[dir][i+ND] - (size_t)gf.bnd[dir];
	k+=2;
      }

    MPI_Datatype bnd_edge_put;
    MPI_Type_create_struct(1+2*n_edge, blen_arr, displ_arr, type_arr, &bnd_edge_put);
    MPI_Type_commit(&bnd_edge_put);
       
    MPI_Comm comm = lat->comms->comm;
    int proc_id = lat->comms->proc_id;
    int to = lat->comms->neigh_proc[(dir+ND) % (2*ND)];
    int from = lat->comms->neigh_proc[dir];
    MPI_Sendrecv(gf.field, 1, bnd_edge_collect, to, to,
		 gf.bnd[dir], 1, bnd_edge_put, from, proc_id, comm,
		 MPI_STATUS_IGNORE);
    MPI_Type_free(&bnd_edge_collect);
    MPI_Type_free(&bnd_edge_put);    
  }
  return;
}

void
qhg_xchange_gauge(qhg_gauge_field gf)
{ 
  for(int i=ND-1; i>=0; i--) {
    get_boundary_edge(gf, i);
    get_boundary_edge(gf, i+ND);
  }
  
  return;
}
