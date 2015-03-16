#ifndef _QHG_NN_THRP_H
#define _QHG_NN_THRP_H 1

#include <qhg_types.h>

qhg_correlator qhg_nn_thrp(qhg_spinor_field [], qhg_spinor_field [], qhg_gauge_field, int [], qhg_thrp_nn_sink_params);
void qhg_write_nn_thrp(char [], qhg_correlator);
  
#endif /* _QHG_NN_THRP_H */
