#ifndef _QHG_APE_SMEAR
#define _QHG_APE_SMEAR 1

#include <qhg_types.h>

void qhg_ape_smear(qhg_gauge_field, qhg_gauge_field, double, int);
void qhg_ape_smear_3d(qhg_gauge_field, qhg_gauge_field, double, int);
void qhg_ape_smear_iter(qhg_gauge_field, qhg_gauge_field, double, int);
void qhg_ape_smear_3d_iter(qhg_gauge_field, qhg_gauge_field, double, int);

#endif /* _QHG_APE_SMEAR */
