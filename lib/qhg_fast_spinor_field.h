#ifndef _QHG_FAST_SPINOR_FIELD_H
#define _QHG_FAST_SPINOR_FIELD_H 1

#include <qhg_types.h>

qhg_fast_spinor_field qhg_fast_spinor_field_init(qhg_lattice *, enum qhg_fermion_bc_time);
void qhg_fast_spinor_field_finalize(qhg_fast_spinor_field );
void qhg_fast_spinor_field_copy(qhg_fast_spinor_field, qhg_fast_spinor_field);
void qhg_fast_spinor_field_import(qhg_fast_spinor_field, qhg_spinor_field[NS*NC]);
double qhg_fast_spinor_field_normsq(qhg_fast_spinor_field);

#endif /* _QHG_FAST_SPINOR_FIELD_H */
