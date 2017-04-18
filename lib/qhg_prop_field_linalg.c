#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <mpi.h>
#include <math.h>

#include <qhg_types.h>
#include <qhg_defs.h>
#include <qhg_prop_gammas.h>

#define VSC(v, sc) (sc + NC*NS*v)

void
qhg_prop_field_Gdag(qhg_spinor_field psi[NC*NS])
{
  qhg_lattice *lat = psi[0].lat;
  unsigned long int lvol = lat->lvol;
  for(unsigned long int v=0; v<lvol; v++) {
    _Complex double G[NS*NC][NS*NC];
    _Complex double D[NS*NC][NS*NC];
    for(int cs0=0; cs0<NS*NC; cs0++)
      for(int cs1=0; cs1<NS*NC; cs1++) {
	G[cs0][cs1] = psi[cs1].field[VSC(v, cs0)];
      }

    for(int cs0=0; cs0<NS*NC; cs0++)
      for(int cs1=0; cs1<NS*NC; cs1++) {
	D[cs0][cs1] = conj(G[cs1][cs0]);
      }
    
    for(int cs0=0; cs0<NS*NC; cs0++)
      for(int cs1=0; cs1<NS*NC; cs1++) {
	psi[cs1].field[VSC(v, cs0)] = D[cs0][cs1];
      }
  }  
  return;
}

void
qhg_prop_field_g5_G(qhg_spinor_field psi[NC*NS])
{
  qhg_lattice *lat = psi[0].lat;
  unsigned long int lvol = lat->lvol;
  for(unsigned long int v=0; v<lvol; v++) {
    _Complex double G[NS*NC][NS*NC];
    _Complex double g5G[NS*NC][NS*NC];
    for(int cs0=0; cs0<NS*NC; cs0++)
      for(int cs1=0; cs1<NS*NC; cs1++) {
	G[cs0][cs1] = psi[cs1].field[VSC(v, cs0)];
      }

    prop_g5_G(g5G, G);
    
    for(int cs0=0; cs0<NS*NC; cs0++)
      for(int cs1=0; cs1<NS*NC; cs1++) {
	psi[cs1].field[VSC(v, cs0)] = g5G[cs0][cs1];
      }
  }  
  return;
}
