#include <mpi.h>
#include <qhg.h>

void
usage(char *argv[])
{
  fprintf(stderr, " Usage: %s CONF_NAME\n", argv[0]);
  return;
}

int
main(int argc, char *argv[])
{   
  if(argc != 2) {
    usage(argv);
    exit(1);
  }

  char *fname;
  asprintf(&fname, "%s", argv[1]);
  
  int dims[ND] = {12, 6, 8, 10}; // t,x,y,z
  qhg_lattice *lat = qhg_lattice_init(dims);
  int am_io_proc = lat->comms->proc_id == 0 ? 1 : 0;
  qhg_gauge_field gf = qhg_gauge_field_init(lat);  
  
  /*
    Read config assuming ildg
  */
  qhg_read_gauge_field_ildg(gf, fname);
  
  /*
    plaquette 
  */
  double p = qhg_plaquette(gf);
  if(am_io_proc)
    printf("Plaquette = %10.8f\n", p);

  qhg_gauge_field_finalize(gf);  
  qhg_lattice_finalize(lat);
  return 0;
}
  
