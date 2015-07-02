#include <string.h>

#include <qhg_alloc.h>
#include <qhg_defs.h>
#include <qhg_types.h>

#define D(sign, dir) ((sign)*ND + (dir))

qhg_spinor_field
qhg_spinor_field_init(qhg_lattice *lat)
{
  int lvol = lat->lvol;
  int *bvol = lat->bvol;
  int ext_vol = lvol;
  for(int i=0; i<ND; i++)
    ext_vol += 2*bvol[i];
  
  qhg_spinor_field sp;
  sp.field = qhg_alloc(ext_vol*NS*NC*sizeof(_Complex double));


  int v_offset = lvol;
  for(int d0=0; d0<ND; d0++)
    for(int s0=0; s0<2; s0++)    
      if(bvol[d0] != 0) {
	sp.bnd[D(s0, d0)] = &sp.field[v_offset*NC*NS];
	v_offset += bvol[d0];
      } else {
	sp.bnd[D(s0, d0)] = NULL;	
      }
    
  sp.lat = lat;
  sp.bc = qhg_alloc(sizeof(_Complex double)*ND);
  for(int d=0; d<ND; d++)
    sp.bc[d] = 1.;
  return sp;
}

void
qhg_spinor_field_set_bc(qhg_spinor_field sp, _Complex double bc[])
{
  for(int d=0; d<ND; d++)
    sp.bc[d] = bc[d];
  return;
}

void
qhg_spinor_field_finalize(qhg_spinor_field sp)
{
  free(sp.field);
  for(int dir=0; dir<2*ND; dir++)
    sp.bnd[dir] = NULL;
  sp.lat = NULL;
  free(sp.bc);
  return;
}

void
qhg_spinor_field_copy(qhg_spinor_field y, qhg_spinor_field x)
{
  y.lat = x.lat;
  memcpy(y.bc, x.bc, sizeof(_Complex double)*ND);
  memcpy(y.field, x.field, x.lat->lvol*NC*NS*sizeof(_Complex double));
  return;
}
