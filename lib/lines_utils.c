#include <lines_types.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

lines
lines_new(int n)
{
  lines ll;
  ll.len = n;
  ll.cur = 0;
  ll.l = malloc(sizeof(line )*n);
  if(ll.l == NULL) {
    fprintf(stderr, "error allocating memory in %s\n", __func__);
    exit(-1);
  }
  return ll;
}

lines
lines_append(lines ll, line *l, int n)
{
  int nn = ll.cur;
  if(nn + n > ll.len) {
    line *new = malloc(sizeof(line)*(nn+n));
    if(new == NULL) {
      fprintf(stderr, "error allocating memory in %s\n", __func__);
      exit(-1);
    }
    memcpy(new, ll.l, sizeof(line)*nn);
    memcpy(new+nn, l, sizeof(line)*n);
    free(ll.l);
    ll.l = new;
    ll.len = n+nn;
  } else {
    memcpy(ll.l+nn, l, sizeof(line)*n);
  }
  ll.cur = n+nn;
  return ll;
}

lines
lines_sorted(lines ll, int step)
{
  int n = ll.cur;
  line *li = ll.l;
  for(int i=0; i<n; i+=step)
    for(int j=i+step; j<n; j+=step)
      if(li[j].n < li[i].n) {
	line swap[step];
	memcpy(swap, &li[i], sizeof(line)*step);
	memcpy(&li[i], &li[j], sizeof(line)*step);
	memcpy(&li[j], swap, sizeof(line)*step);
      }
  return ll;
}

void
lines_del(lines ll)
{
  free(ll.l);
  ll.len = 0;
  ll.cur = 0;
  return;
}
