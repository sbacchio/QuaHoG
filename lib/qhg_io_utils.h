#ifndef _QHG_IO_UTILS_H
#define _QHG_IO_UTILS_H 1
#include <stdio.h>

int qhg_is_bigendian(void);
void qhg_byte_swap_float(float *, size_t);
void qhg_byte_swap_double(double *, size_t);
FILE *qhg_fopen(char [], char *);

#endif /* _QHG_IO_UTILS_H */
