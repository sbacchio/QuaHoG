#include <stdio.h>
#include <qhg_alloc.h>
#include <qhg_defs.h>
#include <qhg_idx.h>
#include <qhg_types.h>
#include <tmLQCD.h>
#ifdef QHG_OMP
#include <omp.h>
#endif

qhg_comms *
qhg_comms_init(qhg_lattice *lat)
{
  MPI_Init(NULL, NULL);

  int verbosity = 3;
  tmLQCD_invert_init(0, NULL, verbosity);

  tmLQCD_lat_params *lp = qhg_alloc(sizeof(tmLQCD_lat_params));
  tmLQCD_mpi_params *mp = qhg_alloc(sizeof(tmLQCD_mpi_params));
  int ret = 0;
  ret = tmLQCD_get_lat_params(lp);
  ret = tmLQCD_get_mpi_params(mp);

  int ldims[ND] = {lp->T, lp->LX, lp->LY, lp->LZ};
  int procs[ND] = {mp->nproc_t, mp->nproc_x, mp->nproc_y, mp->nproc_z};
  int gdims[ND], lvol=1, gvol=1;
  for(int i=0; i<ND; i++) {
    gdims[i] = ldims[i]*procs[i];
    lvol *= ldims[i];
    gvol *= gdims[i];
  }

  for(int dir=0; dir<ND; dir++)
    if(gdims[dir] != lat->dims[dir]) {
      fprintf(stderr, " Dimension missmatch with tmLQCD\n");
      fprintf(stderr, " tmLQCD	: %d*%d %d*%d %d*%d %d*%d\n",
	      lp->T,procs[0],
	      lp->LX,procs[1],
	      lp->LY,procs[2],
	      lp->LZ,procs[3]);
      fprintf(stderr, " QHG	: %d %d %d %d\n",
	      lat->dims[0],
	      lat->dims[1],
	      lat->dims[2],
	      lat->dims[3]);
      exit(1);
    }
  qhg_comms *comms = qhg_alloc(sizeof(qhg_comms));
  for(int dir=0; dir<ND; dir++) {
    comms->proc_dims[dir] = procs[dir];
    comms->proc_coords[dir] = mp->proc_coords[dir];
  }
  comms->proc_id = mp->cart_id;
  comms->nprocs = mp->nproc;

#ifdef QHG_OMP
  comms->nthreads = mp->omp_num_threads;
  omp_set_num_threads(comms->nthreads);
#pragma omp parallel
  {
    int ithr = omp_get_thread_num();
    if(ithr == 0 && comms->proc_id == 0)
      printf(" Numb. threads = %d\n", comms->nthreads);
  }
#else
  comms->nthreads = 1;
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
  tmLQCD_finalise();
  free(comms);
  comms = NULL;
  MPI_Finalize();
  return;
}
