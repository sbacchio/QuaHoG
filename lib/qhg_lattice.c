#include <stdio.h>

#include <qhg_alloc.h>
#include <qhg_defs.h>
#include <qhg_idx.h>
#include <qhg_types.h>
#include <qhg_comms.h>

#ifdef TMLQCD
#include <tmLQCD.h>
#endif 

#define D(sign, dir) ((sign)*ND + (dir))

qhg_lattice *
qhg_lattice_init(int dims[ND], qhg_comms *comms)
{
  qhg_lattice *lat = qhg_alloc(sizeof(qhg_lattice));

  int err = 0;
  for(int d=0; d<ND; d++)
    err += dims[d] % comms->proc_dims[d];

  if(err) {
    if(comms->proc_id == 0) {
      fprintf(stderr, " Lattice dimensions not divisible by processes\n");
      fprintf(stderr, " (%d %d %d %d) / (%d %d %d %d)\n", 
					   dims[0], dims[1], dims[2], dims[3],
					   comms->proc_dims[0], comms->proc_dims[1], comms->proc_dims[2], comms->proc_dims[3]);
    }
    MPI_Abort(MPI_COMM_WORLD, 6);
  }
  /* Global dimensions */
  for(int d=0; d<ND; d++)
    lat->dims[d] = dims[d];
  lat->vol = 1;
  for(int d=0; d<ND; d++)
    lat->vol *= dims[d];
  lat->v3 = lat->vol / dims[0];

  /* Local dimensions */
  int ldims[ND];
  for(int dir=0; dir<ND; dir++)
    ldims[dir] = lat->dims[dir]/comms->proc_dims[dir];
  int lvol = lat->vol/comms->nprocs;
  int lv3 = lvol/ldims[0];

  int *procs = comms->proc_dims;  
  int par_dir[ND];
  /* par_dir[d] == 1 if direction d is parallelised, else 0 */
  for(int d=0; d<ND; d++)
    par_dir[d] = procs[d] == 1 ? 0 : 1;

#ifdef TMLQCD
  tmLQCD_lat_params *lp = qhg_alloc(sizeof(tmLQCD_lat_params));
  int ret = tmLQCD_get_lat_params(lp);
  int gdims[ND];
  int lp_ldims[] = {
    lp->T, lp->LX, lp->LY, lp->LZ
  };
  for(int i=0; i<ND; i++) {
    gdims[i] = lp_ldims[i]*procs[i];
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
  free(lp);
#endif
  
  int bvol[ND];  
  /* 
     The number of sites of the boundary in each direction. For
     directions which are not parallelised bvol[dir] == 0;
  */
  for(int dir=0; dir<ND; dir++)
    bvol[dir] = par_dir[dir] ? lvol/ldims[dir] : 0;

  /* 
     Maps sites to their nearest neighbor. Maps sites in the bulk to
     the boundaries and in the boundaries to the edges.
  */
  int nnv = lvol;
  for(int i=0; i<ND; i++)
    nnv += 2*bvol[i];

  int *nn[2*ND];
  for(int d=0; d<2*ND; d++)
    nn[d] = qhg_alloc(sizeof(int)*nnv);

  /* 
     First initialize assuming periodic boundaries within local
     lattice
  */ 
  for(int dir=0; dir<ND; dir++) {
    /* 
       Neighbors in +dir 
    */
    for(int v=0; v<lvol; v++) {
      int x[ND] = CO(v, ldims);
      x[dir] = (x[dir] + 1) % ldims[dir];
      nn[dir][v] = IDX(x, ldims);
    }
    /* 
       Neighbors in -dir 
    */
    for(int v=0; v<lvol; v++) {
      int x[ND] = CO(v, ldims);
      x[dir] = (ldims[dir] + x[dir] - 1) % ldims[dir];
      nn[dir+ND][v] = IDX(x, ldims);
    }
  }

  /* 
     The number of sites of the edges in each direction pair
     (d0,d1). For a pair of directions to have edges, both should be
     parallelised, otherwise evol[d0][d1] = evol[d1][d0] == 0;
  */
  int evol[ND][ND];  
  for(int d0=0; d0<ND; d0++) {
    for(int d1=0; d1<ND; d1++)
      evol[d0][d1] = (par_dir[d0]*par_dir[d1]) ? lvol/ldims[d0]/ldims[d1] : 0;
    /* the diagonal is always zero */
    evol[d0][d0] = 0;
  }
  
  int v_offset = lvol;
  for(int d0=0; d0<ND; d0++)
    for(int s0=0; s0<2; s0++)    
      if(par_dir[d0]) {
	/* 
	   bdims is the 4D-dimensional extent of the boundary,
	   i.e. with dimension d0 equal to 1
	*/
	int bdims[ND];
	for(int i=0; i<ND; i++)
	  bdims[i] = ldims[i];
	bdims[d0] = 1;      
      
	/* 
	   Index the sites which reside on the boundaries as neighbors of
	   sites in the bulk (interior). When allocating spinor- and
	   gauge-fields these are made to reside in memory adjecent to
	   lvol. First come the forward boundaries for all directions, then
	   the backwards.
	*/ 
	for(int v=0; v<bvol[d0]; v++) {
	  int xb[ND] = CO(v, bdims);
	  int x[ND];
	  for(int i=0; i<ND; i++)
	    x[i] = xb[i];
	  /* 
	     When d < ND, the boundary is neighbor to sites with
	     coordinate ldims[d0] - 1. 
	     
	     When d >= ND, the boundary is neighbor to sites with
	     coordinate 0.
	  */
	  x[d0] = s0 == 0 ? ldims[d0]-1 : 0;
	  int vv = IDX(x, ldims);
	  nn[D(s0, d0)][vv] = v + v_offset;
	}	

	/* 
	   Index the sites which reside on the boundaries as neighbors
	   of sites on the boundaries, for the perpendicular
	   directions to d0.
	*/ 
	for(int v=0; v<bvol[d0]; v++) {
	  int xb[ND] = CO(v, bdims);
	  int x[ND];
	  for(int d1=0; d1<ND; d1++)
	    for(int s1=0; s1<2; s1++) {
	      if(d1 == d0)
		continue;
	      for(int i=0; i<ND; i++)
		x[i] = xb[i];

	      x[d1] = (bdims[d1] + x[d1] + 1 - 2*s1) % bdims[d1];
	      int vv = IDX(x, bdims);
	      nn[D(s1, d1)][v + v_offset] = vv + v_offset;
	    }
	}
	v_offset += bvol[d0]; 
      }

  /* 
     Index the sites which reside on the edges as neighbors of sites
     on the boundary. When allocating spinor- and gauge-fields these
     are made to reside in memory adjecent to the boundaries.
  */ 
  for(int d0=0; d0<ND; d0++)
    for(int s0=0; s0<2; s0++)    
      for(int d1=d0+1; d1<ND; d1++)
	for(int s1=0; s1<2; s1++)
	  if(par_dir[d0]*par_dir[d1]) {
	    /* 
	       4D-dimensions of the boundary, i.e. with dimension d0 = d1 = 1
	    */
	    int edims[ND];
	    for(int i=0; i<ND; i++)
	      edims[i] = ldims[i];
	    edims[d0] = 1;
	    edims[d1] = 1;  	
      
	    for(int v=0; v<evol[d0][d1]; v++) {
	      int xb[ND] = CO(v, edims);
	      int x[ND];
	      for(int i=0; i<ND; i++)
		x[i] = xb[i];
	      x[d0] = s0 == 0 ? ldims[d0] - 1 : 0;
	      x[d1] = s1 == 0 ? ldims[d1] - 1 : 0;
	      int vv = IDX(x, ldims);
	      nn[D(s1, d1)][nn[D(s0, d0)][vv]] = v+v_offset;
	      /* The hop-order commutes */
	      nn[D(s0, d0)][nn[D(s1, d1)][vv]] = v+v_offset;
	    }
	    v_offset += evol[d0][d1];
	  }
  
  for(int dir=0; dir<2*ND; dir++)
    lat->nn[dir] = nn[dir];

  for(int dir=0; dir<ND; dir++) {
    lat->ldims[dir] = ldims[dir];
    lat->bvol[dir] = bvol[dir];
    for(int d=0; d<ND; d++)
      lat->evol[dir][d] = evol[dir][d];
  }
  
  lat->lvol = lvol;
  lat->lv3 = lv3;  
  lat->comms = comms;
  return lat;
}

void
qhg_lattice_finalize(qhg_lattice *lat)
{
  free(lat);
  return;
}

