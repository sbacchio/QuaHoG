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


static int read_fwd_props = 1;
static int write_fwd_props = 0;

static int read_bwd_props = 1;
static int write_bwd_props = 0;

char *
flt_str(double x)
{
  char *s;
  asprintf(&s, "%g", x);
  unsigned long int n = (unsigned long int)index(s, '.');
  if((void *)n == NULL) {
    n = (unsigned long int)index(s, '\0');
    *((char *)n) = 'p';  
    *((char *)n+1) = '0';
    *((char *)n+2) = '\0';    
  } else {
    *((char *)n) = 'p';
  }
  return s;
}

#define NF 2
#define NSRC 2
#define NSNK (8)

enum flavs {
  up,
  dn
};

char flav_str[NF][3] = {"up","dn"};

static char *
proj_to_str(enum projector proj)
{
  switch(proj) {
  case P0:
    return "P0\0";
  case P3:
    return "P3\0";
  case P4:
    return "P4\0";
  case P5:
    return "P5\0";
  case P6:
    return "P6\0";
  }
  return NULL;
}

static enum projector
str_to_proj(char s[])
{
  if(strcmp(s, "P0\0") == 0)
    return P0;
  if(strcmp(s, "P3\0") == 0)
    return P3;
  if(strcmp(s, "P4\0") == 0)
    return P4;
  if(strcmp(s, "P5\0") == 0)
    return P5;
  if(strcmp(s, "P6\0") == 0)
    return P6;
  return -1;
}

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
  int max_mom_sq = 8;  
  int n_ape = 50;
  double alpha_ape = 0.5;
  int n_gauss = 50;
  double alpha_gauss = 4.0;
  int source_coords[NSRC][ND];
  qhg_thrp_nn_sink_params thrp_snk[NSNK];
  
  {
    char *fname_sinks;
    asprintf(&fname_sinks, "sinks.%04d.list", config);  
    FILE *fp = qhg_fopen(fname_sinks, "r");
    for(int i=0; i<NSNK; i++) {
      int dt;
      char proj[256];
      fscanf(fp, "%d %s", &dt, proj);
      thrp_snk[i].dt = dt;
      thrp_snk[i].proj = str_to_proj(proj);
    }
    fclose(fp);
  }

  {
    char *fname_sourcepos;
    asprintf(&fname_sourcepos, "sources.%04d.list", config);  
    FILE *fp = qhg_fopen(fname_sourcepos, "r");
    for(int i=0; i<NSRC; i++) {
      int *s = &source_coords[i][0];
      fscanf(fp, "%d %d %d %d", &s[0], &s[1], &s[2], &s[3]);
    }
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
  qhg_spinor_field sol_sm_u[NS*NC], sol_sm_d[NS*NC];
  for(int sp=0; sp<NS; sp++)
    for(int co=0; co<NC; co++) {
      sol_u[CS(sp,co)] = qhg_spinor_field_init(lat);
      sol_d[CS(sp,co)] = qhg_spinor_field_init(lat);
      sol_sm_u[CS(sp,co)] = qhg_spinor_field_init(lat);
      sol_sm_d[CS(sp,co)] = qhg_spinor_field_init(lat);
    }

  qhg_mom_list mom_list = qhg_mom_list_init(max_mom_sq);  

  /*
    APE string used in filenames
   */
  char *apestr;
  asprintf(&apestr, "aN%da%s", n_ape, flt_str(alpha_ape));

  /*
    Source smearing string used in filenames
  */
  char *smrstr;
  asprintf(&smrstr, "gN%da%s", n_gauss, flt_str(alpha_gauss));
  
  /*
    Loop over source positions
  */  
  for(int isrc=0; isrc<NSRC; isrc++) {
    /* 
       This source position coordinates 
    */
    int *sco = &source_coords[isrc][0];
    if(am_io_proc)
      printf("Source coords (t, x, y, z) = (%d, %d, %d, %d)\n", sco[0], sco[1], sco[2], sco[3]);

    /*
      Source position string to be used in filenames
    */
    char *srcstr;
    asprintf(&srcstr, "sx%02dsy%02dsz%02dst%02d", sco[1], sco[2], sco[3], sco[0]);

    /*
      Jump to reading forward props if available
     */
    if(read_fwd_props)
      goto FREAD;
    
    /*
      Smear the source for this source position
    */
    if(am_io_proc)
      printf("Smearing the source\n");  
    for(int sp=0; sp<NS; sp++)
      for(int co=0; co<NC; co++) {
	if(am_io_proc)
	  printf("\tcol=%d, spin=%d\n", co, sp);  
	qhg_spinor_field aux = qhg_spinor_field_init(lat);
	qhg_point_spinor_field(aux, sco, sp, co);
	src[CS(sp,co)] = qhg_spinor_field_init(lat);
	qhg_gauss_smear(src[CS(sp,co)], aux, gf_ape, alpha_gauss, n_gauss);
	qhg_spinor_field_finalize(aux);
      }
    if(am_io_proc)
      printf("Done\n");  
    
    /*
      Invert for up- and down-quark in the same color-spin loop
    */
    int op[] = {0, 1};
    int write_prop = 0;
    qhg_spinor_field *sol_f[] = {sol_u, sol_d};
    for(int i=0; i<NS*NC; i++)
      for(int flav=0; flav<NF; flav++)
	tmLQCD_invert((double *) sol_f[flav][i].field,
		      (double *) src[i].field,
		      op[flav], write_prop);
    /*
      Write the propagators if selected
    */
    if(write_fwd_props) {
      char *propname;
      
      asprintf(&propname, "prop_%s.up", srcstr);      
      qhg_write_spinors(propname, NC*NS, sol_u);
      free(propname);
      
      asprintf(&propname, "prop_%s.dn", srcstr);      
      qhg_write_spinors(propname, NC*NS, sol_d);
      free(propname);      
    }
        
    /*
      We don't need the source any more
    */
    for(int i=0; i<NS*NC; i++)
      qhg_spinor_field_finalize(src[i]);

  FREAD: if(read_fwd_props) {
      char *propname;
      
      asprintf(&propname, "prop_%s.up", srcstr);      
      qhg_read_spinors(sol_u, NC*NS, propname);
      free(propname);
      
      asprintf(&propname, "prop_%s.dn", srcstr);      
      qhg_read_spinors(sol_d, NC*NS, propname);
      free(propname);      
    }
    
    /*
      Twist t-phase to anti-periodic boundary conditions
    */
    qhg_spinors_untwist_bc(sol_u, NC*NS, 1.0, sco[0]);
    qhg_spinors_untwist_bc(sol_d, NC*NS, 1.0, sco[0]);
    
    /*
      Smear the propagators sink-side
    */
    if(am_io_proc)
      printf("Smearing the propagator sink\n");  
    for(int sp=0; sp<NS; sp++)
      for(int co=0; co<NC; co++) {
	if(am_io_proc)
	  printf("\tcol=%d, spin=%d\n", co, sp);  
	qhg_gauss_smear(sol_sm_u[CS(sp,co)], sol_u[CS(sp,co)], gf_ape, alpha_gauss, n_gauss);
	qhg_gauss_smear(sol_sm_d[CS(sp,co)], sol_d[CS(sp,co)], gf_ape, alpha_gauss, n_gauss);
      }
    if(am_io_proc)
      printf("Done\n");  
	
    /*
      Smeared nucleon and meson correlators and fourier transform
    */
    if(am_io_proc)
      printf("Meson correlator (smeared)\n"); 
    qhg_correlator mesons = qhg_mesons(sol_sm_u, sol_sm_d, sco);
    qhg_correlator mesons_ft = qhg_ft(mesons, &mom_list, "fwd");
    qhg_correlator_finalize(mesons);
    
    if(am_io_proc)
      printf("Nucleon correlator (smeared)\n"); 
    qhg_correlator nucleons = qhg_nucleons(sol_sm_u, sol_sm_d, sco);
    qhg_correlator nucleons_ft = qhg_ft(nucleons, &mom_list, "fwd");
    qhg_correlator_finalize(nucleons);
    
    /*
      Write the correlators
    */
    {
      char *fname;
      asprintf(&fname, "mesons_%s_%s_%s.dat", srcstr, smrstr, apestr);
      qhg_write_mesons(fname, mesons_ft);
      if(am_io_proc)
	printf("Wrote %s\n", fname); 
      free(fname);
    }
    qhg_correlator_finalize(mesons_ft);
    
    {
      char *fname;
      asprintf(&fname, "nucleons_%s_%s_%s.dat", srcstr, smrstr, apestr);
      qhg_write_nucleons(fname, nucleons_ft);
      if(am_io_proc)
	printf("Wrote %s\n", fname); 
      free(fname);
    }
    qhg_correlator_finalize(nucleons_ft);

    /*
      Twist t-phase to back to twisted boundary conditions
    */
    qhg_spinors_untwist_bc(sol_u, NC*NS, -1.0, sco[0]);
    qhg_spinors_untwist_bc(sol_d, NC*NS, -1.0, sco[0]);
    qhg_spinors_untwist_bc(sol_sm_u, NC*NS, -1.0, sco[0]);
    qhg_spinors_untwist_bc(sol_sm_d, NC*NS, -1.0, sco[0]);

    /*
      Allocate sequential source and solution
     */
    qhg_spinor_field seq_src[NC*NS], seq_sol[NC*NS];
    for(int i=0; i<NC*NS; i++) {
      seq_src[i] = qhg_spinor_field_init(lat);
      seq_sol[i] = qhg_spinor_field_init(lat);
    }
    
    /*
      Loop over sequential sink parameters
    */    
    for(int isnk=0; isnk<NSNK; isnk++) {
      for(int flav=0; flav<NF; flav++) {
	/* 
	   Jump to reading backward props if available 
	*/
	if(read_bwd_props)
	  goto BREAD;
	
	switch(flav) {
	case up:
	  qhg_nn_sequential_sink_u(seq_src, sol_sm_u, sol_sm_d, sco[0], thrp_snk[isnk]);
	  break;
	case dn:
	  qhg_nn_sequential_sink_d(seq_src, sol_sm_u, sco[0], thrp_snk[isnk]);
	  break;
	}
	
	/*
	  Smear the sequential source
	*/
	if(am_io_proc)
	  printf("Smearing the sequential sink, isnk = %2d, Proj = %s, sink-source = %2d, flav = %s\n",
		 isnk, proj_to_str(thrp_snk[isnk].proj), thrp_snk[isnk].dt, flav_str[flav]);  
	for(int sp=0; sp<NS; sp++)
	  for(int co=0; co<NC; co++) {
	    if(am_io_proc)
	      printf("\tcol=%d, spin=%d\n", co, sp);  
	    qhg_gauss_smear(seq_src[CS(sp,co)], seq_src[CS(sp,co)], gf_ape, alpha_gauss, n_gauss);
	  }

	int op_id[2] = {1, 0};
	for(int i=0; i<NS*NC; i++)
	  tmLQCD_invert((double *) seq_sol[i].field,
			(double *) seq_src[i].field,
			op_id[flav], write_prop);
	
	/*
	  Write the propagators if selected
	*/
	if(write_bwd_props) {
	  char *propname;	
	  asprintf(&propname, "backprop_%s_%s_dt%02d.%s", srcstr, proj_to_str(thrp_snk[isnk].proj),
		   thrp_snk[isnk].dt, flav_str[flav]);
	  qhg_write_spinors(propname, NC*NS, seq_sol);
	  free(propname);
	}

      BREAD: if(read_bwd_props) {
	  char *propname;
      
	  asprintf(&propname, "backprop_%s_%s_dt%02d.%s", srcstr, proj_to_str(thrp_snk[isnk].proj),
		   thrp_snk[isnk].dt, flav_str[flav]);
	  qhg_read_spinors(seq_sol, NC*NS, propname);
	  free(propname);
	}
	
	/*
	  Untwist boundary condition for both forward and backward prop
	*/
	qhg_spinor_field *fwd[2] = {sol_u, sol_d};
	qhg_spinors_untwist_bc(fwd[flav], NC*NS, 1.0, sco[0]);
	qhg_spinors_untwist_bc(seq_sol, NC*NS, 1.0, sco[0]);

	/*
	  Backward prop takes a 3x contribution of ts-t0 from three combined props
	*/
	double angle = (3.0*thrp_snk[isnk].dt)/dims[0];
	_Complex double phase = cos(angle*M_PI) - _Complex_I*sin(angle*M_PI);
	qhg_spinors_ax(phase, seq_sol, NC*NS);
	/*
	  backprop = (\gamma_5 backprop)^\dagger
	 */
	qhg_prop_field_g5_G(seq_sol);
	qhg_prop_field_Gdag(seq_sol);

	/*
	  Three-point function. Needs gauge-field for derivative
	  operators. 
	 */
	qhg_correlator thrp = qhg_nn_thrp(fwd[flav], seq_sol, gf, sco, thrp_snk[isnk]);

	/* 
	   Fourier transform over momuntum insertion
	*/
	qhg_correlator thrp_ft = qhg_ft(thrp, &mom_list, "bwd");

	/*
	  Write three-point function
	 */
	{
	  char *fname;
	  asprintf(&fname, "thrp_%s_%s_%s_%s_dt%02d.%s",
		   srcstr, smrstr, apestr, proj_to_str(thrp_snk[isnk].proj),
		   thrp_snk[isnk].dt, flav_str[flav]);
	  qhg_write_nn_thrp(fname, thrp_ft);
	  if(am_io_proc)
	    printf("Wrote %s\n", fname); 
	  free(fname);
	}

	/*
	  Clean-up
	 */
	qhg_correlator_finalize(thrp_ft);
	qhg_correlator_finalize(thrp);

	/*
	  Need to twist forward prop, to be used in next sequential
	  sink
	 */
	qhg_spinors_untwist_bc(fwd[flav], NC*NS, -1.0, sco[0]);
      }
    }
      
    /*
      Free sequential source and solution
    */
    for(int i=0; i<NC*NS; i++) {
      qhg_spinor_field_finalize(seq_src[i]);
      qhg_spinor_field_finalize(seq_sol[i]);
    }        
    free(srcstr);
  }
  free(smrstr);
  
  /*
    Destroy momentum list
  */
  qhg_mom_list_finalize(mom_list);
  
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
  