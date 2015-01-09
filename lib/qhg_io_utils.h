#ifndef _QHG_IO_UTILS_H
#define _QHG_IO_UTILS_H 1

int qhg_is_bigendian(void);
void qhg_byte_swap_float(float *, size_t);
void qhg_byte_swap_double(double *, size_t);

#endif /* _QHG_IO_UTILS_H */
