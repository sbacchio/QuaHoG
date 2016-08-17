#include <stdlib.h>
#include <string.h>

#include <mpi.h>
#include <complex.h>
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
    size_t site_size = NC*NS*sizeof(_Complex double);
    MPI_Datatype bnd_collect = get_slice(ND, ldims, dir, sign == +1 ?  0 : ldims[dir]-1, site_size);
    MPI_Comm comm = lat->comms->comm;
    int proc_id = lat->comms->proc_id;
    int fproc = lat->comms->neigh_proc[dir];
    int bproc = lat->comms->neigh_proc[dir+ND];
    if(sign > 0) {
      MPI_Request req;
      /* MPI_Sendrecv(field, 1, bnd_collect, bproc, bproc, */
      /* 		   bnd, bvol*site_size, MPI_BYTE, fproc, proc_id, comm, */
      /* 		   MPI_STATUS_IGNORE); */
      MPI_Irecv(bnd, bvol*site_size, MPI_BYTE, fproc, proc_id, comm, &req);
      MPI_Send(field, 1, bnd_collect, bproc, bproc, comm);
      MPI_Wait(&req, MPI_STATUS_IGNORE);
    } else {
      MPI_Request req;
      /* MPI_Sendrecv(field, 1, bnd_collect, fproc, fproc, */
      /* 		   bnd, bvol*site_size, MPI_BYTE, bproc, proc_id, comm, */
      /* 		   MPI_STATUS_IGNORE); */
      MPI_Irecv(bnd, bvol*site_size, MPI_BYTE, bproc, proc_id, comm, &req);
      MPI_Send(field, 1, bnd_collect, fproc, fproc, comm);
      MPI_Wait(&req, MPI_STATUS_IGNORE);
    }
    MPI_Type_free(&bnd_collect);
  }
  return;
}

void
qhg_xchange_spinor(qhg_spinor_field sp)
{
  for(int i=0; i<2*ND; i++)
    get_boundary(sp.bnd[i], i, sp.field, sp.lat);
  
  return;
}
