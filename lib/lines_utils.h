#ifndef _LINES_UTILS_H
#define _LINES_UTILS_H 1
#include <lines_types.h>

lines lines_new(int);
lines lines_append(lines, line *, int);
lines lines_sorted(lines, int);
void lines_del(lines);
#endif /* _LINES_UTILS_H */
