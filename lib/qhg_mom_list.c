#include <math.h>

#include <qhg_types.h>
#include <qhg_alloc.h>

qhg_mom_list 
qhg_mom_list_init(int max_mom_sq)
{
  qhg_mom_list mom_list;
  mom_list.max_mom_sq = max_mom_sq;
  int m = (int)(sqrt(max_mom_sq)+0.5);
  /* Find the number of momenta with square <= max_mom_sq */
  int nmom = 0;
  for(int z=-m; z<=m; z++)
    for(int y=-m; y<=m; y++)
      for(int x=-m; x<=m; x++) {
	if(x*x + y*y + z*z <= max_mom_sq)
	  nmom++;
      }
  /* Allocate array to store momenta with square <= max_mom_sq */
  int (*mom)[4] = qhg_alloc(sizeof(int)*4*nmom);
  nmom = 0;
  for(int z=-m; z<=m; z++)
    for(int y=-m; y<=m; y++)
      for(int x=-m; x<=m; x++) {
	int sq = x*x + y*y + z*z;
	if(sq <= max_mom_sq) {
	  mom[nmom][0] = sq;
	  mom[nmom][1] = x;
	  mom[nmom][2] = y;
	  mom[nmom][3] = z;
	  nmom++;
	}
      }
  /* Sort momentum vectors ascending */  
  for(int i=0; i<nmom; i++) {
    for(int j=i; j<nmom; j++)
      if(mom[j][0] <= mom[i][0]) {
	for(int k=0; k<4; k++) {
	  int swap = mom[j][k];
	  mom[j][k] = mom[i][k];
	  mom[i][k] = swap;
	}
      }
  }
  mom_list.mom_vecs = qhg_alloc(sizeof(int)*3*nmom);
  for(int i=0; i<nmom; i++)
    for(int j=0; j<3; j++)
      mom_list.mom_vecs[i][j] = mom[i][j+1];
  free(mom);
  mom_list.n_mom_vecs = nmom;
  return mom_list;
}

void
qhg_mom_list_finalize(qhg_mom_list mom_list)
{
  free(mom_list.mom_vecs);
  return;
}
