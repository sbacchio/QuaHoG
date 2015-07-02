#ifndef _QHG_TYPES_H
#define _QHG_TYPES_H 1

#include <mpi.h>
#include <qhg_defs.h>


typedef struct {
  int proc_id;
  int nprocs;
  int nthreads;
  int proc_dims[ND];
  int proc_coords[ND];
  int neigh_proc[2*ND];
  MPI_Comm comm;
} qhg_comms;

typedef struct {
  /* Global dimensions */
  int dims[ND];
  int vol;
  int v3;
  /* Local dimensions */
  int ldims[ND];
  int lvol;
  int lv3;
  /* Volume of boundaries */
  int bvol[ND];
  /* Volume of edges */
  int evol[ND][ND];
  /* Nearest neighbor indexing */
  int *nn[2*ND];
  /* Pointer to comms struct */
  qhg_comms *comms;
} qhg_lattice;

typedef struct {
  _Complex double *field;
  _Complex double *bnd[2*ND];
  _Complex double *edge[2*ND][2*ND];
  qhg_lattice *lat;
} qhg_gauge_field;

typedef struct {
  _Complex double *field;
  _Complex double *bnd[2*ND];
  _Complex double *edge[2*ND][2*ND];
  _Complex double *bc;
  qhg_lattice *lat;
} qhg_spinor_field;

typedef struct {
  int (*mom_vecs)[3];
  int n_mom_vecs;
  int max_mom_sq;
} qhg_mom_list;
  
typedef struct {
  _Complex double *C;
  int site_size;
  int *origin;
  /* So far cutoff is not used */
  int cutoff[ND];
  qhg_mom_list *mom_list;
  qhg_lattice *lat;  
} qhg_correlator;

enum projector {
  P0,
  P3,
  P4,
  P5,
  P6
};

enum flavor {
  up,
  dn
};

typedef struct {
  qhg_correlator corr;
  int dt;
  enum flavor flav;
  enum projector proj;
} qhg_thrp_correlator;

typedef struct {
  enum projector proj;
  int dt;
} qhg_thrp_nn_sink_params;

#endif /* _QHG_TYPES_H */
