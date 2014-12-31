#ifndef _QHG_GAUGE_FIELD_H
#define _QHG_GAUGE_FIELD_H 1

#include <qhg_types.h>

qhg_gauge_field qhg_gauge_field_init(qhg_lattice *);
void qhg_gauge_field_finalize(qhg_gauge_field );
void qhg_gauge_field_copy(qhg_gauge_field, qhg_gauge_field);

#endif /* _QHG_GAUGE_FIELD_H */
