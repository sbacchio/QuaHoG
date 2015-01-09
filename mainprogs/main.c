#include <mpi.h>
#include <tmLQCD.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <qhg.h>

int
main(int argc, char *argv[])
{
  int dims[] = {24,12,12,12};
  int n_ape = 50;
  double alpha_ape = 0.5;
  int n_gauss = 50;
  double alpha_gauss = 4.0;
  int source_coords[] = {6,3,10,5}; // t,x,y,z
  int config = 3144;

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
  
  qhg_write_spinors("prop.up", 12, sol_u);
  qhg_write_spinors("prop.dn", 12, sol_d);  

  for(int i=0; i<NS*NC; i++)
    qhg_spinor_field_finalize(src[i]);

  
  MPI_Finalize();
  return 0;

READ:
  qhg_read_spinors(sol_u, 12, "prop.up");
  qhg_read_spinors(sol_d, 12, "prop.dn");   
  
  qhg_correlator mesons = qhg_mesons(sol_u, sol_d, source_coords);
  int max_mom_sq = 4;
  qhg_mom_list mom_list = qhg_mom_list_init(max_mom_sq);
  qhg_correlator mesons_ft = qhg_ft(mesons, &mom_list, "fwd");

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
   
  { /* Print the two-point function */
    char fname[] = "mesons.dat";
    int am_io_proc = lat->comms->proc_id == 0 ? 1 : 0;
    FILE *fp = NULL;
    if(am_io_proc)
      fp = fopen(fname, "w");
    
    int site_size = mesons_ft.site_size;
    int nm = mesons_ft.mom_list->n_mom_vecs;    
    int (*mv)[3] = mesons_ft.mom_list->mom_vecs;
    int *pd = lat->comms->proc_dims;
    int *ld = lat->ldims;
    int *d = lat->dims;
    for(int tt=0; tt<lat->dims[0]; tt++) {    
      int t = (d[0] + tt + mesons_ft.origin[0]) % d[0];
      int proc_t = t/lat->comms->proc_dims[0];
      for(int n=0; n<nm; n++) {
	int *k = mv[n];
	/* Global coordinates */
	int x[] = {t,
		   (k[0]+d[1]) % d[1],
		   (k[1]+d[2]) % d[2],
		   (k[2]+d[3]) % d[3]};

	/* Process coordinates */
	int px[] = {x[0]/ld[0],
		    x[1]/ld[1],
		    x[2]/ld[2],
		    x[3]/ld[3]};

	/* local coordinates */
	int lx[] = {x[0]%ld[0],
		    x[1]%ld[1],
		    x[2]%ld[2],
		    x[3]%ld[3]};

	int lv = IDX(lx, ld);
	int proc = IDX(px, pd);
	if(am_io_proc)
	  fprintf(fp, " %2d %+d %+d %+d   %+e %+e   %+e %+e LS\n",
		  tt, k[0], k[1], k[2],
		  creal(mesons_ft.C[lv*site_size + 0]),
		  cimag(mesons_ft.C[lv*site_size + 0]),
		  creal(mesons_ft.C[lv*site_size + 1]),
		  cimag(mesons_ft.C[lv*site_size + 1])
		  );
	if(am_io_proc)
	  fprintf(fp, " %2d %+d %+d %+d   %+e %+e   %+e %+e SS\n",
		  tt, k[0], k[1], k[2],
		  creal(mesons_sm_ft.C[lv*site_size + 0]),
		  cimag(mesons_sm_ft.C[lv*site_size + 0]),
		  creal(mesons_sm_ft.C[lv*site_size + 1]),
		  cimag(mesons_sm_ft.C[lv*site_size + 1])
		  );
	MPI_Barrier(lat->comms->comm);
      }
    }
    if(am_io_proc)
      fclose(fp);
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
