#include <stdio.h>
#include <qhg_alloc.h>
#include <qhg_defs.h>
#include <qhg_idx.h>
#include <qhg_types.h>
#include <tmLQCD.h>

qhg_comms *
qhg_comms_init(qhg_lattice *lat)
{
  MPI_Init(NULL, NULL);

  int verbosity = 3;
  tmLQCD_invert_init(0, NULL, verbosity);

  tmLQCD_lat_params lp;
  tmLQCD_mpi_params mp;   
  tmLQCD_get_lat_params(&lp);
  tmLQCD_get_mpi_params(&mp);  
  
  int ldims[ND] = {lp.T, lp.LX, lp.LY, lp.LZ};  
  int procs[ND] = {mp.nproc_t, mp.nproc_x, mp.nproc_y, mp.nproc_z};
  int gdims[ND], lvol=1, gvol=1;
  for(int i=0; i<ND; i++) {
    gdims[i] = ldims[i]*procs[i];
    lvol *= ldims[i];
    gvol *= gdims[i];
  }
  
  for(int dir=0; dir<ND; dir++)
    if(gdims[dir] != lat->dims[dir]) {
      fprintf(stderr, " Dimension missmatch with tmLQCD\n");
      exit(1);
    }
  qhg_comms *comms = qhg_alloc(sizeof(qhg_comms));  
  for(int dir=0; dir<ND; dir++) {
    comms->proc_dims[dir] = procs[dir];
    comms->proc_coords[dir] = mp.proc_coords[dir];
  }
  comms->proc_id = mp.cart_id;
  comms->nprocs = mp.nproc;
   
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

