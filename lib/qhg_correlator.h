#ifndef _QHG_CORRELATOR_H
#define _QHG_CORRELATOR_H 1

#include <qhg_types.h>

qhg_correlator qhg_correlator_init(int, qhg_lattice *);
qhg_correlator qhg_correlator_copy(qhg_correlator);
void qhg_correlator_finalize(qhg_correlator);

#endif /* _QHG_CORRELATOR_H */
