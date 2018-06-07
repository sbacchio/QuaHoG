#include <string.h>
#include <complex.h>
#include <qhg_alloc.h>
#include <qhg_defs.h>
#include <qhg_idx.h>
#include <qhg_types.h>

qhg_fast_spinor_field
qhg_fast_spinor_field_init(qhg_lattice *lat, enum qhg_fermion_bc_time bc)
{
  unsigned long int lvol = lat->lvol;
  unsigned long int lv3 = lat->lv3;
  int lt = lat->ldims[0];
  
  qhg_fast_spinor_field sp;
  sp.alloc = qhg_alloc(lvol*NS*NS*NC*NC*2*sizeof(afloat));
  sp.field = qhg_alloc(lt*NS*NS*NC*NC*2*sizeof(afloat *));
  sp.lat = lat;
  sp.bc = bc;

  afloat *tmp = sp.alloc;

  for(int t=0; t<lt; t++)
    for(int s0=0; s0<NS; s0++)
      for(int s1=0; s1<NS; s1++)
        for(int c0=0; c0<NC; c0++)
          for(int c1=0; c1<NC; c1++)
            for(int ri=0; ri<2; ri++) {
              sp.field[t][s0][s1][c0][c1][ri] = tmp;
              tmp+=lv3;
            }

  return sp;
}

void
qhg_fast_spinor_field_finalize(qhg_fast_spinor_field sp)
{
  free(sp.alloc);
  free(sp.field);
  sp.lat = NULL;
  return;
}

void
qhg_fast_spinor_field_copy(qhg_fast_spinor_field y, qhg_fast_spinor_field x)
{
  y.lat = x.lat;
  y.bc = x.bc;
  memcpy(y.alloc, x.alloc, x.lat->lvol*NS*NS*NC*NC*2*sizeof(afloat));
  return;
}

void
qhg_fast_spinor_field_import(qhg_fast_spinor_field y, qhg_spinor_field x[NS*NC])
{
  y.lat = x[0].lat;
  y.bc = x[0].bc;
  unsigned long int lvol = y.lat->lvol;
  unsigned long int lv3 = y.lat->lv3;

  for(int cs0=0; cs0<NS*NC; cs0++)
    for(unsigned long int v=0; v<lvol; v++)
      for(int cs1=0; cs1<NC*NS; cs1++) {
        y.field[v/lv3][CS2S(cs0)][CS2S(cs1)][CS2C(cs0)][CS2C(cs1)][0][v%lv3] = creal(x[cs0].field[(v*NC*NS + cs1)]);
        y.field[v/lv3][CS2S(cs0)][CS2S(cs1)][CS2C(cs0)][CS2C(cs1)][1][v%lv3] = cimag(x[cs0].field[(v*NC*NS + cs1)]);
      }
  return;
}

#define SQUARE(a) (a)*(a)

double
qhg_fast_spinor_field_normsq(qhg_fast_spinor_field y)
{
  int lt = y.lat->ldims[0];
  unsigned long int lv3 = y.lat->lv3;

  double norm = 0;

  for(int t=0; t<lt; t++)
    for(int s0=0; s0<NS; s0++)
      for(int s1=0; s1<NS; s1++)
        for(int c0=0; c0<NC; c0++)
          for(int c1=0; c1<NC; c1++)
            for(int ri=0; ri<2; ri++)
              for(unsigned long int v=0; v<lv3; v++) {
                norm += SQUARE(y.field[t][s0][s1][c0][c1][ri][v]);
              }

  MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_SUM, y.lat->comms->comm);
  return norm;
}

