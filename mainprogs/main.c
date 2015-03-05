#include <mpi.h>
#include <tmLQCD.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <qhg.h>

static int read_props = 0;
static int write_props = 1;

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

#define NSRC 5
#define NSMR 2

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
  int config = atoi(argv[1]);
  int dims[ND] = {6, 4, 4, 4}; // t,x,y,z
  int n_ape = 50;
  double alpha_ape = 0.5;
  int n_gausses[NSMR];
  double alpha_gausses[NSMR];
  int source_coords[NSRC][ND];
  {
    char *fname_sourcepos;
    asprintf(&fname_sourcepos, "sources.%04d.list", config);  
    FILE *fp = fopen(fname_sourcepos, "r");
    for(int i=0; i<NSRC; i++) {
      int *s = &source_coords[i][0];
      fscanf(fp, "%d %d %d %d", &s[0], &s[1], &s[2], &s[3]);
      /* printf("%d %d %d %d\n", */
      /* 	     source_coords[i][0], */
      /* 	     source_coords[i][1], */
      /* 	     source_coords[i][2], */
      /* 	     source_coords[i][3]); */
    }
    fclose(fp);
  }

  {
    char *fname_smearings;
    asprintf(&fname_smearings, "smearings.%04d.list", config);
    FILE *fp = fopen(fname_smearings, "r");
    for(int i=0; i<NSMR; i++) {
      int *n = &n_gausses[i];
      double *a = &alpha_gausses[i];
      fscanf(fp, "%d %lf", n, a);
      //printf("%d %g\n", n_gausses[i], alpha_gausses[i]);
    }
    fclose(fp);
  }
  int max_mom_sq = 8;
  
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
    Loop over smearings
   */
  for(int ismr=0; ismr<NSMR; ismr++) {
    double alpha_gauss = alpha_gausses[ismr];
    int n_gauss = n_gausses[ismr];
    /*
      Source smearing string used in filenames
    */
    char *smr_0_str;
    asprintf(&smr_0_str, "g0N%da%s", n_gauss, flt_str(alpha_gauss));
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
	for(int flav=0; flav<2; flav++)
	  tmLQCD_invert((double *) sol_f[flav][i].field,
			(double *) src[i].field,
			op[flav], write_prop);
      /*
	Write the propagators if selected
      */
      if(write_props) {
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

      /*
	Twist t-phase to anti-periodic boundary conditions
      */
      qhg_spinors_untwist_bc(sol_u, NC*NS, 1.0, sco[0]);
      qhg_spinors_untwist_bc(sol_d, NC*NS, 1.0, sco[0]);

      /* 
	 Loop over sink smearings
      */
      for(int ismr_snk=0; ismr_snk<NSMR; ismr_snk++) {
	double alpha_gauss_snk = alpha_gausses[ismr_snk];
	int n_gauss_snk = n_gausses[ismr_snk];    
	/*
	  Sink smearing string used in filenames
	*/
	char *smr_x_str;
	asprintf(&smr_x_str, "gxN%da%s", n_gauss_snk, flt_str(alpha_gauss_snk));
	
	/*
	  Smear the propagators sink-side
	*/
	if(am_io_proc)
	  printf("Smearing the propagator sink\n");  
	for(int sp=0; sp<NS; sp++)
	  for(int co=0; co<NC; co++) {
	    if(am_io_proc)
	      printf("\tcol=%d, spin=%d\n", co, sp);  
	    qhg_gauss_smear(sol_sm_u[CS(sp,co)], sol_u[CS(sp,co)], gf_ape, alpha_gauss_snk, n_gauss_snk);
	    qhg_gauss_smear(sol_sm_d[CS(sp,co)], sol_d[CS(sp,co)], gf_ape, alpha_gauss_snk, n_gauss_snk);
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
	  if(am_io_proc)
	    printf("Writing the meson correlator\n"); 
	  char *fname;
	  asprintf(&fname, "mesons_%s_%s_%s_%s.dat", srcstr, smr_x_str, smr_0_str, apestr);
	  qhg_write_mesons(fname, mesons_ft);
	  if(am_io_proc)
	    printf("Done\n");
	  free(fname);
	}
	qhg_correlator_finalize(mesons_ft);
	
	{
	  if(am_io_proc)
	    printf("Writing the nucleon correlator\n"); 
	  char *fname;
	  asprintf(&fname, "nucleons_%s_%s_%s_%s.dat", srcstr, smr_x_str, smr_0_str, apestr);
	  qhg_write_nucleons(fname, nucleons_ft);
	  if(am_io_proc)
	    printf("Done\n");
	  free(fname);
	}
	qhg_correlator_finalize(nucleons_ft);
	free(smr_x_str);
      }
      free(srcstr);
    }
    free(smr_0_str);    
  }
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
  
