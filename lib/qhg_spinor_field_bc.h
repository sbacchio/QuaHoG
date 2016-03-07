#ifndef _QHG_SPINOR_FIELD_BC_H
#define _QHG_SPINOR_FIELD_BC_H 1
#include <complex.h>
#include <qhg_types.h>

void qhg_spinor_twist_t_bc(qhg_spinor_field , qhg_spinor_field , double);
void qhg_spinors_set_bc(qhg_spinor_field *, int, _Complex double []);

#endif /* _QHG_SPINOR_FIELD_BC_H */
