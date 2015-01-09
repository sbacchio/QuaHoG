#include <stdlib.h>
#include <complex.h>
#include <qhg_defs.h>
#include <qhg_types.h>
#include <qhg_correlator.h>
#include <qhg_prop_ops.h>

#define NCHANNELS 2 /* +/- charged pions */
#define VC(v, c) ((v)*NCHANNELS + (c))
#define VSC(v, cs) ((v)*NC*NS + (cs))

qhg_correlator
qhg_mesons(qhg_spinor_field sp_u[NS*NC], qhg_spinor_field sp_d[NS*NC], int source_coords[ND])
{
  qhg_lattice *lat = sp_u[0].lat;
  qhg_correlator corr = qhg_correlator_init(NCHANNELS, lat);
  int lvol = lat->lvol;
  for(int v=0; v<lvol; v++) {
    _Complex double U[NS*NC][NS*NC];
    _Complex double D[NS*NC][NS*NC];    
    _Complex double C[NS*NC][NS*NC];    
    for(int cs0=0; cs0<NS*NC; cs0++)
      for(int cs1=0; cs1<NS*NC; cs1++) {
	U[cs0][cs1] = sp_u[cs1].field[VSC(v, cs0)];
	D[cs0][cs1] = sp_d[cs1].field[VSC(v, cs0)];	
      }
    prop_mul_gd(C, U, U);
    corr.C[VC(v, 0)] = prop_trace(C);

    prop_mul_gd(C, D, D);    
    corr.C[VC(v, 1)] = prop_trace(C);
  }
  
  for(int i=0; i<ND; i++)
    corr.origin[i] = source_coords[i];
  corr.mom_list = NULL;
  return corr;
}
