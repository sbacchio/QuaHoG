#include <stdio.h>
#include <qhg_alloc.h>
#include <qhg_defs.h>
#include <qhg_idx.h>
#include <qhg_types.h>
#ifdef TMLQCD
#include <tmLQCD.h>
#endif
#ifdef QHG_OMP
#include <omp.h>
#endif

qhg_comms *
qhg_comms_init(int *proc_geom)
{
  MPI_Init(NULL, NULL);

  qhg_comms *comms = qhg_alloc(sizeof(qhg_comms));
  comms->nthreads = 1;
  
#ifdef TMLQCD
  int verbosity = 3;
  tmLQCD_invert_init(0, NULL, verbosity);

  tmLQCD_mpi_params *mp = qhg_alloc(sizeof(tmLQCD_mpi_params));
  int ret = 0;
  ret = tmLQCD_get_mpi_params(mp);  
  
  int procs[ND] = {mp->nproc_t, mp->nproc_x, mp->nproc_y, mp->nproc_z};
  for(int dir=0; dir<ND; dir++) {
    comms->proc_dims[dir] = procs[dir];
    comms->proc_coords[dir] = mp->proc_coords[dir];
  }

  comms->proc_id = mp->cart_id;
  comms->nprocs = mp->nproc;

  if(proc_geom != NULL) {
    if(proc_geom[0] != comms->proc_dims[0] ||
       proc_geom[1] != comms->proc_dims[1] ||
       proc_geom[2] != comms->proc_dims[2] ||
       proc_geom[3] != comms->proc_dims[3]) {
      if(comms->proc_id == 0) {
	fprintf(stderr, " process geometry missmatch with tmLQCD\n");	
      }
      MPI_Abort(MPI_COMM_WORLD, 5);
    }
  }
  
#ifdef QHG_OMP
  comms->nthreads = mp->omp_num_threads;
#endif
  free(mp);  
#else /* not TMLQCD */
  MPI_Comm_rank(MPI_COMM_WORLD, &comms->proc_id);
  MPI_Comm_size(MPI_COMM_WORLD, &comms->nprocs);
  if(proc_geom == NULL) {
    /* Use MPI's Dims_create() if proc_geom == NULL */
    for(int i=0; i<ND; i++)
      comms->proc_dims[i] = 0;
    MPI_Dims_create(comms->nprocs, ND, comms->proc_dims);
  } else {
    /* If proc_geom != NULL, expect an array of 4 ints, the
       process-grid dimensions */
    int n = 1;
    for(int dir=0; dir<ND; dir++) {
      comms->proc_dims[dir] = proc_geom[dir];
      n *= proc_geom[dir];
    }

    if(n != comms->nprocs) {
      if(comms->proc_id == 0)
	fprintf(stderr, " nprocs != prod(proc_dims), quitting\n");
      MPI_Abort(MPI_COMM_WORLD, 4);
    }
  }
  
  int periods[] = {1,1,1,1};
  MPI_Cart_create(MPI_COMM_WORLD, ND, comms->proc_dims, periods, 0, &comms->comm);
  MPI_Cart_coords(comms->comm, comms->proc_id, ND, comms->proc_coords);
#ifdef QHG_OMP
  comms->nthreads = omp_get_num_threads();
#endif
  int *procs = comms->proc_dims;
#endif
  
  for(int dir=0; dir<ND; dir++) {
    int c[ND];

    /* Neighber in +dir */
    for(int i=0; i<ND; i++)
      c[i] = comms->proc_coords[i];
    c[dir] = (c[dir] + 1) % procs[dir];
    int rp = IDX(c,procs);

    /* Neighber in -dir */
    for(int i=0; i<ND; i++)
      c[i] = comms->proc_coords[i];
    c[dir] = (procs[dir] + c[dir] - 1) % procs[dir];
    int rm = IDX(c,procs);

    comms->neigh_proc[dir] = rp;
    comms->neigh_proc[dir+ND] = rm;
  }
  comms->comm = MPI_COMM_WORLD;
  return comms;
}

void
qhg_comms_finalize(qhg_comms *comms)
{
#ifdef TMLQCD
  tmLQCD_finalise();
#endif
  free(comms);
  comms = NULL;
  MPI_Finalize();
  return;
}
