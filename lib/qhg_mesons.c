#include <stdlib.h>
#include <complex.h>
#include <qhg_defs.h>
#include <qhg_types.h>
#include <qhg_correlator.h>

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
    corr.C[VC(v, 0)] = 0.;
    corr.C[VC(v, 1)] = 0.;    
    for(int cs0=0; cs0<NS*NC; cs0++)
      for(int cs1=0; cs1<NS*NC; cs1++) {
	_Complex double u = sp_u[cs0].field[VSC(v, cs1)];
	_Complex double d = sp_d[cs0].field[VSC(v, cs1)];	
	corr.C[VC(v, 0)] += creal(u)*creal(u) + cimag(u)*cimag(u);
	corr.C[VC(v, 1)] += creal(d)*creal(d) + cimag(d)*cimag(d);	
      }
  }
  for(int i=0; i<ND; i++)
    corr.origin[i] = source_coords[i];
  corr.mom_list = NULL;
  return corr;
}
