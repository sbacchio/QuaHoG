#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>

#include <complex.h>
#include <qhg_io_utils.h>
#include <qhg_idx.h>
#include <qhg_defs.h>
#include <qhg_types.h>
#include <qhg_alloc.h>
#include <qhg_nucleon_defs.h>
#include <qhg_correlator.h>
#include <qhg_correlator_shift.h>

#if 0 /* Deprecated function for writing to text files */
#include <lines_types.h>
#include <lines_utils.h>

void
qhg_write_mom_nucleons(char fname[], qhg_correlator corr)
{
  qhg_lattice *lat = corr.lat;
  int proc_id = lat->comms->proc_id;
  int am_io_proc = proc_id == 0 ? 1 : 0;
  int site_size = corr.site_size;
  int nm = corr.mom_list->n_mom_vecs;
  int lvol = lat->lvol;
  int (*mv)[3] = corr.mom_list->mom_vecs;
  int *pd = lat->comms->proc_dims;
  int np = lat->comms->nprocs;
  int *ld = lat->ldims;
  int *d = lat->dims;
  int nlines[np];
  for(int i=0; i<np; i++)
    nlines[i] = 0;
  
  lines lines_loc = lines_new(d[0]*nm*NFLAV*NCHAN*NCOMP*NCOMP);
  for(int t=0; t<d[0]; t++) {    
    for(int n=0; n<nm; n++) {
      /* Momentum vector */
      int *k = mv[n];
      /* Global coordinates of momentum vector */
      int x[] = {t, (k[0]+d[1]) % d[1], (k[1]+d[2]) % d[2], (k[2]+d[3]) % d[3]};
      /* Coordinates of process on which x[] resides */
      int pc[] = {x[0]/ld[0], x[1]/ld[1], x[2]/ld[2], x[3]/ld[3]};
      /* Lexico index of pc[] */
      int ip = IDX(pc, pd);
      /* Local coordinates of x[] within process ip */
      int lc[] = {x[0]%ld[0], x[1]%ld[1], x[2]%ld[2], x[3]%ld[3]};
      /* Lexico index of lc[] */
      int ix = IDX(lc,ld);
      _Complex double *c = corr.C;
      line li[NFLAV*NCHAN*NCOMP];
      /* Shift to time relative to source */
      int tt = (d[0] + t - corr.origin[0]) % d[0];      
      for(int ifl=0; ifl<NFLAV; ifl++)
	for(int ich=0; ich<NCHAN; ich++)
	  for(int si0=0; si0<NCOMP; si0++) {
	    int j = si0 + NCOMP*(ich + NCHAN*ifl);
	    li[j].n = j + NFLAV*NCHAN*NCOMP*(n + nm*tt);
	    sprintf(li[j].c, "%4d %+d %+d %+d  %+e %+e  %+e %+e  %+e %+e  %+e %+e  %s %s\n",
		    tt, k[0], k[1], k[2],
		    creal(c[NIDX(ix, ifl, ich, si0, 0)]), cimag(c[NIDX(ix, ifl, ich, si0, 0)]),
		    creal(c[NIDX(ix, ifl, ich, si0, 1)]), cimag(c[NIDX(ix, ifl, ich, si0, 1)]),
		    creal(c[NIDX(ix, ifl, ich, si0, 2)]), cimag(c[NIDX(ix, ifl, ich, si0, 2)]),
		    creal(c[NIDX(ix, ifl, ich, si0, 3)]), cimag(c[NIDX(ix, ifl, ich, si0, 3)]),
		    flav_tags[ifl],
		    chan_tags[ich]);
	  }
      nlines[ip] += NFLAV*NCHAN*NCOMP;
      if(proc_id == ip)
	lines_loc = lines_append(lines_loc, li, NFLAV*NCHAN*NCOMP);
      
    }
  }
  
  int nl = 0;
  int displs[np];
  for(int i=0; i<np; i++)
    displs[i] = 0;
  
  for(int i=0; i<np; i++) {
    nl += nlines[i];
    nlines[i] *= sizeof(line);
    for(int j=0; j<i; j++)
      displs[i] += nlines[j];
  }
  
  lines lines_glob = lines_new(nl);
  MPI_Gatherv(lines_loc.l, nlines[proc_id], MPI_BYTE,
	      lines_glob.l, nlines, displs, MPI_BYTE, 0, lat->comms->comm);
  lines_glob.cur = lines_glob.len;
  lines_glob = lines_sorted(lines_glob, NFLAV*NCHAN*NCOMP);

  if(am_io_proc) {
    FILE *fp = qhg_fopen(fname, "w");
    for(int i=0; i<nl; i++)
      fprintf(fp, "%s", lines_glob.l[i].c);
    fclose(fp);
  }
  
  lines_del(lines_loc);
  lines_del(lines_glob);  
  return;
}
#endif

void
qhg_write_mom_nucleons(char fname[], qhg_correlator corr_in)
{
  qhg_lattice *lat = corr_in.lat;
  int proc_id = lat->comms->proc_id;
  int am_io_proc = proc_id == 0 ? 1 : 0;
  int site_size = corr_in.site_size;
  int (*mv)[3] = corr_in.mom_list->mom_vecs;
  int nm = corr_in.mom_list->n_mom_vecs;
  int lvol = lat->lvol;
  int *pd = lat->comms->proc_dims;
  int *pc = lat->comms->proc_coords;
  int np = lat->comms->nprocs;
  int *ld = lat->ldims;
  int *d = lat->dims;
  int t_rank = pc[0];
  /*
    This function requires that the correlator is shifted to the
    origin. Copy into new correlator struct and shift to the origin.
   */
  qhg_correlator corr = qhg_correlator_copy(corr_in);
  qhg_correlator_shift(corr, corr.origin);
  
  hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
  H5Pclose(fapl_id);

  /*
    Attributes (metadata) are: 1) the origin (source position) and 2)
    the index order in the file
   */
  hsize_t n = 4;
  hid_t attrdat_id = H5Screate_simple(1, &n, NULL);
  hid_t attr_id = H5Acreate2(file_id, "Origin", H5T_NATIVE_INT, attrdat_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr_id, H5T_NATIVE_INT, corr.origin);
  H5Aclose(attr_id);
  H5Sclose(attrdat_id);

  char order[] = "C-order: [t,spin0,spin1,real/imag]\0";
  attrdat_id = H5Screate(H5S_SCALAR);
  hid_t type_id = H5Tcopy(H5T_C_S1);
  H5Tset_size(type_id, strlen(order));
  attr_id = H5Acreate1(file_id, "IndexOrder", type_id, attrdat_id, H5P_DEFAULT);
  H5Awrite(attr_id, type_id, &order);

  H5Aclose(attr_id);
  H5Tclose(type_id);
  H5Sclose(attrdat_id);
  
  double *buf = qhg_alloc(sizeof(double)*NS*NS*2*ld[0]);
  _Complex double *c = corr.C;
  /* 
     Hierarchy is /channel/interpolator/ with channel one of ppm or
     pmm (proton or neutron), and with interpolator one of 1-1, 1-2,
     2-1, 2-2 (chi_1, chi_2 combinations)
  */
  for(int ifl=0; ifl<NFLAV; ifl++) {
    char *group_tag;
    /* 
       Group /flav 
    */
    asprintf(&group_tag, "%s", flav_tags[ifl]);
    hid_t group1_id = H5Gcreate(file_id, group_tag, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);    
    for(int ich=0; ich<NCHAN; ich++) {
      /* 
	 Group /flav/chan 
      */
      asprintf(&group_tag, "%s", chan_tags[ich]);      
      hid_t group2_id = H5Gcreate(group1_id, group_tag, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);      
      for(int n=0; n<nm; n++) {
	/* Momentum vector */
	int *k = mv[n];	
	/* 
	   Group /flav/chan/mom
	*/
	asprintf(&group_tag, "px%+dpy%+dpz%+d", k[0], k[1], k[2]);      
	hid_t group3_id = H5Gcreate(group2_id, group_tag, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	/*
	  All ranks must call the dataset create function
	 */
	hsize_t dims[] = {d[0], NS, NS, 2};	  
	hid_t filespace = H5Screate_simple(4, dims, NULL); 
	hid_t dataset_id = H5Dcreate(group3_id, "corr_p", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	/*
	  Subspace of this rank 
	 */
	hsize_t ldims[] = {ld[0], NS, NS, 2};
	hsize_t starts[] = {pc[0]*ld[0], 0, 0, 0};	
	hid_t subspace = H5Screate_simple(4, ldims, NULL);
	filespace = H5Dget_space(dataset_id);
	H5Sselect_hyperslab(filespace, H5S_SELECT_SET, starts, NULL, ldims, NULL);
	hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);      

	/* 
	   The "3d volume rank" of this process
	 */
	int p3d[ND-1] = {pd[1], pd[2], pd[3]};
	int p3c[ND-1] = {pc[1], pc[2], pc[3]};
	int v3_rank = IDX3(p3c, p3d);

	/* 
	   The "3d volume rank" of the processes which this momentum
	   resides on
	 */
	int q[] = {(k[0]+d[1]) % d[1], (k[1]+d[2]) % d[2], (k[2]+d[3]) % d[3]};
	/* process x, y, z coords on which this momentum resides */
	int q3c[] = {q[0]/ld[1], q[1]/ld[2], q[2]/ld[3]};
	int q3_rank = IDX3(q3c, p3d);
	if(q3_rank == v3_rank) {
	  for(int lt=0; lt<ld[0]; lt++) {
	    /* Local coordinates of q[] within process */
	    int lq[] = {lt, q[0]%ld[1], q[1]%ld[2], q[2]%ld[3]};
	    /* Lexico index of lq[] */
	    int iq = IDX(lq,ld);
	    _Complex double *c = corr.C;
	    for(int si0=0; si0<NCOMP; si0++) 
	      for(int si1=0; si1<NCOMP; si1++) {
		buf[0 + 2*(si1 + NS*(si0 + lt*NS))] = creal(corr.C[NIDX(iq, ifl, ich, si0, si1)]);
		buf[1 + 2*(si1 + NS*(si0 + lt*NS))] = cimag(corr.C[NIDX(iq, ifl, ich, si0, si1)]);
	      }
	  }
	  herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, subspace, filespace, plist_id, buf);
	}
	H5Dclose(dataset_id);
	H5Sclose(filespace);
	H5Sclose(subspace);
	H5Pclose(plist_id);	
	H5Gclose(group3_id);
      }
      H5Gclose(group2_id);
    }
    H5Gclose(group1_id);
  }
  H5Fclose(file_id);
  free(buf);

  qhg_correlator_finalize(corr);
  return;
}
