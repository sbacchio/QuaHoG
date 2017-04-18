#include <string.h>

#include <qhg_alloc.h>
#include <qhg_defs.h>
#include <qhg_types.h>

qhg_correlator
qhg_correlator_init(size_t site_size, qhg_lattice *lat)
{
  unsigned long int lvol = lat->lvol;  
  qhg_correlator corr;
  corr.C = qhg_alloc(lvol*site_size*sizeof(_Complex double));
  memset(corr.C, '\0', lvol*site_size*sizeof(_Complex double));
  corr.lat = lat;
  corr.site_size = site_size;
  corr.origin = qhg_alloc(ND*sizeof(int));
  for(int i=0; i<ND; i++)
    corr.origin[i] = 0;
  return corr;
}

qhg_correlator
qhg_correlator_copy(qhg_correlator x)
{
  qhg_lattice *lat = x.lat;
  unsigned long int lvol = lat->lvol;  
  size_t site_size = x.site_size;
  qhg_correlator y = qhg_correlator_init(site_size, lat);  
  memcpy(y.C, x.C, lvol*site_size*sizeof(_Complex double));
  for(int i=0; i<ND; i++)
    y.origin[i] = x.origin[i];
  return y;
}

void
qhg_correlator_finalize(qhg_correlator corr)
{
  free(corr.C);
  corr.lat = NULL;
  free(corr.origin);
  corr.origin = NULL;
  return;
}
