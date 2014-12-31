#ifndef _QHG_COMMS_H
#define _QHG_COMMS_H 1

#include <qhg_types.h>

qhg_comms *qhg_comms_init(qhg_lattice *);
void qhg_comms_finalize(qhg_comms *);

#endif /* _QHG_COMMS_H */
