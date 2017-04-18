#include <string.h>

#include <qhg_alloc.h>
#include <qhg_defs.h>
#include <qhg_types.h>

#define D(sign, dir) ((sign)*ND + (dir))

qhg_gauge_field
qhg_gauge_field_init(qhg_lattice *lat)
{
  unsigned long int lvol = lat->lvol;
  unsigned long int *bvol = lat->bvol;
  unsigned long int (*evol)[ND] = lat->evol;  
  unsigned long int ext_vol = lvol;
  for(int i=0; i<ND; i++) {
    ext_vol += 2*bvol[i];    
    for(int j=i+1; j<ND; j++)
      ext_vol += 4*evol[i][j];
  }
  
  qhg_gauge_field gf;
  gf.field = qhg_alloc(ext_vol*NC*NC*ND*sizeof(_Complex double));

  unsigned long int v_offset = lvol;
  for(int d0=0; d0<ND; d0++)
    for(int s0=0; s0<2; s0++)    
      if(bvol[d0] != 0) {
	gf.bnd[D(s0, d0)] = &gf.field[v_offset*NC*NC*ND];
	v_offset += bvol[d0];
      } else {
	gf.bnd[D(s0, d0)] = NULL;	
      }
	
  for(int d0=0; d0<ND; d0++)
    for(int s0=0; s0<2; s0++)
      for(int d1=d0+1; d1<ND; d1++)
	if(evol[d0][d1] != 0) {
	  for(int s1=0; s1<2; s1++) {
	    gf.edge[D(s0, d0)][D(s1, d1)] = &gf.field[v_offset*NC*NC*ND];
	    v_offset += evol[d0][d1];
	  }
	} else {
	  for(int s1=0; s1<2; s1++)	      
	    gf.edge[D(s0, d0)][D(s1, d1)] = NULL;	
	}
    
  gf.lat = lat;
  return gf;
}

void
qhg_gauge_field_finalize(qhg_gauge_field gf)
{
  free(gf.field);
  for(int dir=0; dir<2*ND; dir++)
    gf.bnd[dir] = NULL;
  gf.lat = NULL;
  return;
}

void
qhg_gauge_field_copy(qhg_gauge_field y, qhg_gauge_field x)
{
  y.lat = x.lat;
  memcpy(y.field, x.field, x.lat->lvol*NC*NC*ND*sizeof(_Complex double));
  return;
}
