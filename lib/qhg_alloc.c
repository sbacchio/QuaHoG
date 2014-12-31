#include <stdlib.h>
#include <stdio.h>

/*
 * Malloc, with minimal error handling
 */
void *
qhg_alloc(size_t size)
{
  void *ptr = malloc(size);
  if(ptr == NULL) {
    fprintf(stderr, " malloc() returned NULL. Out of memory?\n");
    exit(-1);
  }
  return ptr;
}
