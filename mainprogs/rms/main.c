#include <mpi.h>
#include <tmLQCD.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <qhg.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

void
usage(char *argv[])
{
  fprintf(stderr, " Usage: %s CONF_NUMB\n", argv[0]);
  return;
}

int
main(int argc, char *argv[])
{   
  if(argc != 2) {
    usage(argv);
    exit(1);
  }
  char *e;
  int config = (int)strtoul(argv[1], &e, 10);
  if(*e != '\0') {
    usage(argv);
    exit(2);
  }
  
  int dims[ND] = {6, 4, 4, 4}; // t,x,y,z
  int max_mom_sq = 1;  
  int n_ape = 50;
  double alpha_ape = 0.5;
  int n_gauss = 50;
  double alpha_gauss = 4.0;
  int source_coords[ND];

  {
    char *fname_sourcepos;
    asprintf(&fname_sourcepos, "source.%04d.coords", config);  
    FILE *fp = qhg_fopen(fname_sourcepos, "r");
    int *s = &source_coords[0];
    fscanf(fp, "%d %d %d %d", &s[0], &s[1], &s[2], &s[3]);
    fclose(fp);
  }

  qhg_lattice *lat = qhg_lattice_init(dims);
  int am_io_proc = lat->comms->proc_id == 0 ? 1 : 0;
  qhg_gauge_field gf = qhg_gauge_field_init(lat);  
  
  /*
    Use tmLQCD to read config 
  */
  tmLQCD_read_gauge(config);
  _Complex double *gptr;
  tmLQCD_get_gauge_field_pointer((double **)&gptr);

  /*
    import to qhg
  */
  qhg_import_gauge_field(gf, gptr);
  
  /*
    plaquette 
  */
  double p = qhg_plaquette(gf);
  if(am_io_proc)
    printf("Plaquette = %10.8f\n", p);

  /*
    APE smear in 3-dimensions
   */
  double t0 = qhg_stop_watch(0);
  qhg_gauge_field gf_ape = qhg_gauge_field_init(lat);
  qhg_ape_smear_3d(gf_ape, gf, alpha_ape, n_ape);
  if(am_io_proc)
    printf("3D APE smear in %g sec\n", qhg_stop_watch(t0));
  
  /*
    plaquette of smeared gauge-field
  */
  double p_ape = qhg_plaquette(gf_ape);
  if(am_io_proc)
    printf("3D APE plaquette = %10.8f\n", p_ape);
  
  /*
    Spinor fields
   */
  qhg_spinor_field A;
  qhg_spinor_field B;
  A = qhg_spinor_field_init(lat);
  B = qhg_spinor_field_init(lat);
      
  /*
    Smear computing the source r.m.s in each iteration
  */
  t0 = qhg_stop_watch(0);
  if(am_io_proc)
    printf("Smearing the source\n");  

  qhg_point_spinor_field(A, source_coords, 0, 0);

  qhg_spinor_field S[2] = {A, B};
  for(int i=0; i<n_gauss; i++) {
    double *rms;
    qhg_gauss_smear_iter(S[(i+1)%2], S[i%2], gf_ape, alpha_gauss);
    rms = qhg_spinor_field_rms(S[(i+1)%2], source_coords);
    if(am_io_proc)
      printf("Smearing iter = %4d, r.m.s = %+e\n", i, rms[source_coords[0]]);
    free(rms);
  }
    if(am_io_proc)
      printf("Done smearing in %g sec\n", qhg_stop_watch(t0));  
  
  /* 
     Free spinor- and gauge-fields
  */
  qhg_spinor_field_finalize(A);
  qhg_spinor_field_finalize(B);
  qhg_gauge_field_finalize(gf);  
  qhg_gauge_field_finalize(gf_ape);  
  qhg_lattice_finalize(lat);
  return 0;
}
  
