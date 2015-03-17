#include <sys/time.h>

double
qhg_stop_watch(double t0)
{
  double time;
  struct timeval t;
  gettimeofday(&t, NULL);
  time = (double) t.tv_sec + (double) t.tv_usec * 1e-6;
  return (time - t0);
}
