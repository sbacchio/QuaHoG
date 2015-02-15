#ifndef _LINES_TYPES_H
#define _LINES_TYPES_H 1

typedef struct {
  int n;
  char c[256];
} line;

typedef struct {
  line *l;
  int len;
  int cur;
} lines;

#endif /* _LINES_TYPES_H */
