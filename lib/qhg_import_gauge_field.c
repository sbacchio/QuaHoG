#include <stdlib.h>
#include <string.h>

#include <string.h>
#include <complex.h>
#include <qhg_defs.h>
#include <qhg_types.h>

void
qhg_import_gauge_field(qhg_gauge_field gf, _Complex double *ptr)
{
  unsigned long int lvol = gf.lat->lvol;
  memcpy(gf.field, ptr, sizeof(_Complex double)*NC*NC*ND*lvol);
  return;
}
