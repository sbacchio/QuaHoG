#include <mpi.h>
#include <tmLQCD.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <qhg.h>

static int read_props = 1;
static int write_props = 0;

int
main(int argc, char *argv[])
{   
  int dims[] = {24,12,12,12};
  int n_ape = 0;
  double alpha_ape = 0.5;
  int n_gauss = 0;
  double alpha_gauss = 4.0;
  int source_coords[] = {6,3,10,5}; // t,x,y,z
  int config = 3144;
  int max_mom_sq = 16;
  
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
  qhg_gauge_field gf_ape = qhg_gauge_field_init(lat);
  qhg_ape_smear_3d(gf_ape, gf, alpha_ape, n_ape);
  
  /*
    plaquette of smeared gauge-field
  */
  double p_ape = qhg_plaquette(gf_ape);
  if(am_io_proc)
    printf("Plaquette = %10.8f\n", p_ape);

  
  /*
    Source and solution spinor field
   */
  qhg_spinor_field src[NS*NC];
  qhg_spinor_field sol_u[NS*NC], sol_d[NS*NC];
  for(int sp=0; sp<NS; sp++)
    for(int co=0; co<NC; co++) {
      sol_u[CS(sp,co)] = qhg_spinor_field_init(lat);
      sol_d[CS(sp,co)] = qhg_spinor_field_init(lat);
    }

  if(read_props)
    goto READ;
  
  for(int sp=0; sp<NS; sp++)
    for(int co=0; co<NC; co++) {
      qhg_spinor_field aux = qhg_spinor_field_init(lat);
      qhg_point_spinor_field(aux, source_coords, sp, co);
      src[CS(sp,co)] = qhg_spinor_field_init(lat);
      qhg_gauss_smear(src[CS(sp,co)], aux, gf_ape, alpha_gauss, n_gauss);
      qhg_spinor_field_finalize(aux);
    }

  /*
    Invert for solution, first for up-quark (operator 0)
   */
  int op = 0;
  int write_prop = 0;
  for(int i=0; i<NS*NC; i++)
    tmLQCD_invert((double *) sol_u[i].field,
  		  (double *) src[i].field,
  		  op, write_prop);
  
  /*
    Now for down-quark (operator 1)
   */
  op = 1;
  for(int i=0; i<NS*NC; i++)
    tmLQCD_invert((double *) sol_d[i].field,
  		  (double *) src[i].field,
  		  op, write_prop);

  if(write_props) {
    qhg_write_spinors("prop.up", NC*NS, sol_u);
    qhg_write_spinors("prop.dn", NC*NS, sol_d);  
    
    for(int i=0; i<NS*NC; i++)
      qhg_spinor_field_finalize(src[i]);
    
    MPI_Finalize();
    return 0;
  }
  
READ:
  if(am_io_proc)
    printf("Read\n");
  if(read_props) {
    qhg_read_spinors(sol_u, NC*NS, "prop.up");
    qhg_read_spinors(sol_d, NC*NS, "prop.dn");   
  }
  if(am_io_proc)
    printf("Untwist\n");
  qhg_spinors_untwist_bc(sol_u, NC*NS, 1.0, source_coords[0]);
  qhg_spinors_untwist_bc(sol_d, NC*NS, 1.0, source_coords[0]);

  qhg_correlator mesons = qhg_mesons(sol_u, sol_d, source_coords);
  qhg_correlator nucleons = qhg_nucleons(sol_u, sol_d, source_coords);

  qhg_mom_list mom_list = qhg_mom_list_init(max_mom_sq);

  qhg_correlator mesons_ft = qhg_ft(mesons, &mom_list, "fwd");
  qhg_correlator nucleons_ft = qhg_ft(nucleons, &mom_list, "fwd");
  
  qhg_spinor_field sol_sm_u[NS*NC], sol_sm_d[NS*NC];
  for(int sp=0; sp<NS; sp++)
    for(int co=0; co<NC; co++) {
      sol_sm_u[CS(sp,co)] = qhg_spinor_field_init(lat);
      sol_sm_d[CS(sp,co)] = qhg_spinor_field_init(lat);
      qhg_gauss_smear(sol_sm_u[CS(sp,co)], sol_u[CS(sp,co)], gf_ape, alpha_gauss, n_gauss);
      qhg_gauss_smear(sol_sm_d[CS(sp,co)], sol_d[CS(sp,co)], gf_ape, alpha_gauss, n_gauss);
    }
  
  qhg_correlator mesons_sm = qhg_mesons(sol_sm_u, sol_sm_d, source_coords);
  qhg_correlator mesons_sm_ft = qhg_ft(mesons_sm, &mom_list, "fwd");
  qhg_correlator nucleons_sm = qhg_nucleons(sol_sm_u, sol_sm_d, source_coords);
  qhg_correlator nucleons_sm_ft = qhg_ft(nucleons_sm, &mom_list, "fwd");

  {
    if(am_io_proc)
      printf("Write mesons (local)\n");
    char fname[] = "mesons_sl.dat";
    qhg_write_mesons(fname, mesons_ft);
  }

  {
    if(am_io_proc)
      printf("Write mesons (smeared)\n");
    char fname[] = "mesons_ss.dat";
    qhg_write_mesons(fname, mesons_sm_ft);
  }

  {
    if(am_io_proc)
      printf("Write nucleons (local)\n");
    char fname[] = "nucleons_sl.dat";
    qhg_write_nucleons(fname, nucleons_ft);
  }

  {
    if(am_io_proc)
      printf("Write nucleons (smeared)\n");
    char fname[] = "nucleons_ss.dat";
    qhg_write_nucleons(fname, nucleons_sm_ft);
  }
  
  /*
    Destroy momentum list
   */
  qhg_mom_list_finalize(mom_list);
 
  
  /* 
     Destroy correlators 
  */
  qhg_correlator_finalize(mesons);
  qhg_correlator_finalize(mesons_ft);
  qhg_correlator_finalize(mesons_sm);
  qhg_correlator_finalize(mesons_sm_ft);
  qhg_correlator_finalize(nucleons);
  qhg_correlator_finalize(nucleons_ft);
  qhg_correlator_finalize(nucleons_sm);
  qhg_correlator_finalize(nucleons_sm_ft);
  
  
  /* 
     Destroy spinor- and gauge-fields 
  */
  for(int i=0; i<NS*NC; i++) {
    qhg_spinor_field_finalize(sol_u[i]);
    qhg_spinor_field_finalize(sol_d[i]);
    qhg_spinor_field_finalize(sol_sm_u[i]);
    qhg_spinor_field_finalize(sol_sm_d[i]);
  }
  qhg_gauge_field_finalize(gf);  
  qhg_gauge_field_finalize(gf_ape);  
  qhg_lattice_finalize(lat);
  return 0;
}
