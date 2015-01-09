#ifndef _QHG_SPINOR_FIELD_IO_H
#define _QHG_SPINOR_FIELD_IO_H 1

#include <qhg_types.h>

void qhg_write_spinors(char [], int, qhg_spinor_field *);
void qhg_read_spinors(qhg_spinor_field *, int, char []);

#endif /* _QHG_SPINOR_FIELD_IO_H */
