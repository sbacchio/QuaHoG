#include <string.h>

#include <qhg_alloc.h>
#include <qhg_defs.h>
#include <qhg_types.h>

qhg_correlator
qhg_correlator_init(int site_size, qhg_lattice *lat)
{
  int lvol = lat->lvol;  
  qhg_correlator corr;
  corr.C = qhg_alloc(lvol*site_size*sizeof(_Complex double));
  memset(corr.C, '\0', lvol*site_size*sizeof(_Complex double));
  corr.lat = lat;
  corr.site_size = site_size;
  return corr;
}

void
qhg_correlator_finalize(qhg_correlator corr)
{
  free(corr.C);
  corr.lat = NULL;
  return;
}
