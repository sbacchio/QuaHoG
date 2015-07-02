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
#include <qhg_nn_thrp_defs.h>
#include <qhg_correlator.h>
#include <qhg_correlator_shift.h>

void
qhg_write_mom_nn_thrp(char fname[], qhg_thrp_correlator corr_thrp_in)
{
  qhg_lattice *lat = corr_thrp_in.corr.lat;
  int proc_id = lat->comms->proc_id;
  int am_io_proc = proc_id == 0 ? 1 : 0;
  int site_size = corr_thrp_in.corr.site_size;
  int (*mv)[3] = corr_thrp_in.corr.mom_list->mom_vecs;
  int nm = corr_thrp_in.corr.mom_list->n_mom_vecs;
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
  qhg_correlator corr = qhg_correlator_copy(corr_thrp_in.corr);
  qhg_correlator_shift(corr, corr.origin);

  hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
  H5Pclose(fapl_id);
  
  /*
    Attributes (metadata) are:
    1) the origin (source position) [1-dimensional integer array of 4 elements]
    2) the index order in the file [string]
    3) source-sink separation [single integer]
    4) flavor [string]
    5) projector [string]
   */

  {
    hsize_t n = 4;
    hid_t attrdat_id = H5Screate_simple(1, &n, NULL);
    hid_t attr_id = H5Acreate2(file_id, "Origin", H5T_NATIVE_INT, attrdat_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_INT, corr.origin);
    H5Aclose(attr_id);
    H5Sclose(attrdat_id);
  }

  {
    char order[] = "C-order: [t,real/imag]\0";
    hid_t attrdat_id = H5Screate(H5S_SCALAR);
    hid_t type_id = H5Tcopy(H5T_C_S1);
    H5Tset_size(type_id, strlen(order));
    hid_t attr_id = H5Acreate1(file_id, "IndexOrder", type_id, attrdat_id, H5P_DEFAULT);
    H5Awrite(attr_id, type_id, &order);
    H5Aclose(attr_id);
    H5Tclose(type_id);
    H5Sclose(attrdat_id);
  }

  {
    hid_t attrdat_id = H5Screate(H5S_SCALAR);
    hid_t attr_id = H5Acreate1(file_id, "Sink-source separation", H5T_NATIVE_INT, attrdat_id, H5P_DEFAULT);
    H5Awrite(attr_id, H5T_NATIVE_INT, &corr_thrp_in.dt);
    H5Aclose(attr_id);
    H5Sclose(attrdat_id);
  }
  
  {
    char *flav;
    switch(corr_thrp_in.flav) {
    case up:
      flav = strdup("up");
      break;
    case dn:
      flav = strdup("dn");
      break;
    }      
    hid_t attrdat_id = H5Screate(H5S_SCALAR);
    hid_t type_id = H5Tcopy(H5T_C_S1);
    H5Tset_size(type_id, strlen(flav));
    hid_t attr_id = H5Acreate1(file_id, "Flavor", type_id, attrdat_id, H5P_DEFAULT);
    H5Awrite(attr_id, type_id, flav);
    H5Aclose(attr_id);
    H5Tclose(type_id);
    H5Sclose(attrdat_id);
  }
  
  {
    char *proj = strdup(proj_to_str(corr_thrp_in.proj));
    hid_t attrdat_id = H5Screate(H5S_SCALAR);
    hid_t type_id = H5Tcopy(H5T_C_S1);
    H5Tset_size(type_id, strlen(proj));
    hid_t attr_id = H5Acreate1(file_id, "Projector", type_id, attrdat_id, H5P_DEFAULT);
    H5Awrite(attr_id, type_id, proj);
    H5Aclose(attr_id);
    H5Tclose(type_id);
    H5Sclose(attrdat_id);
  }
    
  /*
    Hierarchy is /insertion/ with insertions listed in the
    qhg_nn_thrp_defs.h header file
  */

  double *buf = qhg_alloc(sizeof(double)*2*ld[0]);
  _Complex double *c = corr.C;
  for(int ichan=0; ichan<NCHAN; ichan++) {
    char *group_tag;
    /* 
       Group /insertion 
    */
    asprintf(&group_tag, "%s", chan_tags[ichan]);
    hid_t group_id = H5Gcreate(file_id, group_tag, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);        
    for(int n=0; n<nm; n++) {
      /* Momentum vector */
      int *k = mv[n];	
      /* 
	 Group /insertion/momentum
      */
      asprintf(&group_tag, "px%+dpy%+dpz%+d", k[0], k[1], k[2]);      
      hid_t group1_id = H5Gcreate(group_id, group_tag, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      /*
	All ranks must call the dataset create function
      */
      int nt = corr_thrp_in.dt+1;
      hsize_t dims[] = {nt, 2};	  
      hid_t filespace = H5Screate_simple(2, dims, NULL); 
      hid_t dataset_id = H5Dcreate(group1_id, "corr_p", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      /*
	Subspace of this rank 
      */
      hsize_t ldims[] = {ld[0], 2};
      hsize_t starts[] = {pc[0]*ld[0], 0};
      /* Ranks beyond the tsink should write nothing */
      if(starts[0] > nt)
	ldims[0] = 0;
      /* The ranks sitting on tsink should only write the piece up to tsink */
      if((starts[0] < nt) && (starts[0]+ldims[0] > nt))
	ldims[0] = nt % ld[0];
	
      hid_t subspace = H5Screate_simple(2, ldims, NULL);
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
	  buf[0 + 2*lt] = creal(c[TIDX(iq, ichan)]);
	  buf[1 + 2*lt] = cimag(c[TIDX(iq, ichan)]);
	}
	herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, subspace, filespace, plist_id, buf);
      }
      H5Dclose(dataset_id);
      H5Sclose(filespace);
      H5Sclose(subspace);
      H5Pclose(plist_id);	
      H5Gclose(group1_id);
    }
    H5Gclose(group_id);
  }
  H5Fclose(file_id);
  free(buf);
  
  qhg_correlator_finalize(corr);
  return;
}
