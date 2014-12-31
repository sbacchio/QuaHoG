#ifndef _QHG_GAUSS_SMEAR
#define _QHG_GAUSS_SMEAR 1

#include <qhg_types.h>

void qhg_gauss_smear(qhg_spinor_field, qhg_spinor_field, qhg_gauge_field, double, int);
void qhg_gauss_smear_iter(qhg_spinor_field, qhg_spinor_field, qhg_gauge_field, double);

#endif /* _QHG_GAUSS_SMEAR */
