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
#include <qhg_meson_defs.h>
#include <qhg_correlator.h>
#include <qhg_correlator_shift.h>

void
qhg_write_mom_mesons(char fname[], qhg_correlator corr_in)
{
  qhg_lattice *lat = corr_in.lat;
  int proc_id = lat->comms->proc_id;
  int am_io_proc = proc_id == 0 ? 1 : 0;
  int site_size = corr_in.site_size;
  int (*mv)[3] = corr_in.mom_list->mom_vecs;
  int nm = corr_in.mom_list->n_mom_vecs;
  int lvol = lat->lvol;
  int *pd = lat->comms->proc_dims;
  int pd3[] = {pd[1], pd[2], pd[3]};
  int *pc = lat->comms->proc_coords;
  int np = lat->comms->nprocs;
  int *ld = lat->ldims;
  int ld3[] = {ld[1], ld[2], ld[3]};
  int lv3 = ld3[0]*ld3[1]*ld3[2];
  int *d = lat->dims;
  int t_rank = pc[0];
  int v3_rank = pc[3] + pd[3]*(pc[2] + pd[2]*pc[1]);
  int np3 = pd[3]*pd[2]*pd[1];

  /*
    This function requires that the correlator is shifted to the
    origin. Copy into new correlator struct and shift to the origin.
   */
  qhg_correlator corr = qhg_correlator_copy(corr_in);
  qhg_correlator_shift(corr, corr.origin);

  /* 
     Loop over momenta 
       Get rank of momentum
       Append to datatype for sendtype of this rank with displacement of momentum in corr.C
       Append to datatype of recvtype array: recvtype[np3], with displacement within the nm momenta
       
     Loop over ranks_3d: rank_3 = 0; rank_3<np3; rank_3++
       if rank == rank_3, send corr.C to root with sendtype
       if rank == root, recv with recvtype[rank_3]
  */
  
  /* Count the number of momenta on each rank */
  int send_displ[nm];
  int recv_displ[np3][nm];
  int send_count = 0;
  int recv_count[np3];
  for(int i=0; i<np3; i++)
    recv_count[i] = 0;
  for(int im=0; im<nm; im++) {
    /* Momentum vector */
    int *k = mv[im];	
    /* 
       The "3d volume rank" of the processes which this momentum
       resides on
    */
    int q[] = {(k[0]+d[1]) % d[1], (k[1]+d[2]) % d[2], (k[2]+d[3]) % d[3]};
    /* process x, y, z coords on which this momentum resides */
    int q3c[] = {q[0]/ld[1], q[1]/ld[2], q[2]/ld[3]};
    int q3_rank = IDX3(q3c, pd3);

    /* Local coordinates of q[] within process */
    int lq[] = {q[0]%ld[1], q[1]%ld[2], q[2]%ld[3]};
    /* Lexico index of lq[] */
    int iq = IDX3(lq,ld3);
    if(v3_rank == q3_rank) {
      send_displ[send_count] = iq;
      send_count++;
    }
    recv_displ[q3_rank][recv_count[q3_rank]] = im;
    recv_count[q3_rank]++;
  }

  /* The send datatype */
  MPI_Datatype sendtype;
  MPI_Datatype temptype;
  int blocklens[send_count];
  for(int i=0; i<send_count; i++) {
    blocklens[i] = site_size*sizeof(_Complex double);
    send_displ[i] *= site_size*sizeof(_Complex double);
  }
  MPI_Type_indexed(send_count, blocklens, send_displ, MPI_BYTE, &temptype);
  MPI_Type_commit(&temptype);
  MPI_Type_create_hvector(ld[0], 1, sizeof(_Complex double)*site_size*lv3, temptype, &sendtype);
  MPI_Type_commit(&sendtype);
  MPI_Type_free(&temptype);
  
  /* The recv datatype */
  MPI_Datatype recvtype[np3];
  for(int p=0; p<np3; p++) {
    int blocklens[recv_count[p]];
    for(int i=0; i<recv_count[p]; i++) {
      blocklens[i] = site_size*sizeof(_Complex double);
      recv_displ[p][i] *= site_size*sizeof(_Complex double);
    }
    MPI_Type_indexed(recv_count[p], blocklens, recv_displ[p], MPI_BYTE, &temptype);
    MPI_Type_commit(&temptype);
    MPI_Type_create_hvector(ld[0], 1, sizeof(_Complex double)*site_size*nm, temptype, &recvtype[p]);
    MPI_Type_commit(&recvtype[p]);
    MPI_Type_free(&temptype);
  }
  
  /*
    Split the communicator into nproc_t communicators, one
    communicator for each time-slab
   */
  MPI_Comm t_comm;
  MPI_Comm_split(MPI_COMM_WORLD, t_rank, proc_id, &t_comm);
  _Complex double *buf;
  if(v3_rank == 0)
    buf = qhg_alloc(nm*site_size*ld[0]*sizeof(_Complex double));    

  MPI_Request req;
  MPI_Isend(corr.C, 1, sendtype, 0, v3_rank, t_comm, &req);

  if(v3_rank == 0)
    for(int i=0; i<np3; i++)
      MPI_Recv(buf, 1, recvtype[i], i, i, t_comm, MPI_STATUS_IGNORE);

  MPI_Wait(&req, MPI_STATUS_IGNORE);

  MPI_Type_free(&sendtype);
  for(int i=0; i<np3; i++)
    MPI_Type_free(&recvtype[i]);
  
  MPI_Group world_group;
  MPI_Comm_group(MPI_COMM_WORLD, &world_group);
  /*
    Group including all ranks with v3_rank == 0, which will handle I/O
   */
  MPI_Group io_group;
  int v3_root_ranks[pd[0]];
  for(int i=0; i<pd[0]; i++)
    v3_root_ranks[i] = i*np3;

  MPI_Group_incl(world_group, pd[0], v3_root_ranks, &io_group);
  MPI_Comm io_comm;
  MPI_Comm_create(MPI_COMM_WORLD, io_group, &io_comm);

  if(v3_rank == 0) {
    hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl_id, io_comm, MPI_INFO_NULL);
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

    char order[] = "C-order: [t, real/imag]\0";
    attrdat_id = H5Screate(H5S_SCALAR);
    hid_t type_id = H5Tcopy(H5T_C_S1);
    H5Tset_size(type_id, strlen(order));
    attr_id = H5Acreate1(file_id, "IndexOrder", type_id, attrdat_id, H5P_DEFAULT);
    H5Awrite(attr_id, type_id, &order);

    H5Aclose(attr_id);
    H5Tclose(type_id);
    H5Sclose(attrdat_id);
  
    double *subdat = qhg_alloc(sizeof(double)*2*ld[0]);
    _Complex double *c = corr.C;
    /* 
       Hierarchy is /flavor/interpolator/momentum with flavor one of up or dn,
       and with interpolator one of the gamma combinations in qhg_meson_defs.h
    */
    for(int ifl=0; ifl<NFLAV; ifl++) {
      char *group_tag;
      /* 
	 Group /flav 
      */
      asprintf(&group_tag, "%s", flav_tags[ifl]);
      hid_t group_id = H5Gcreate(file_id, group_tag, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);    
      for(int interp=0; interp<NGAMMAS; interp++) {
	/* 
	   Group /flav/interpolator
	*/
	asprintf(&group_tag, "%s", gamma_tags[interp]);      
	hid_t group1_id = H5Gcreate(group_id, group_tag, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);      
	for(int im=0; im<nm; im++) {
	  /* Momentum vector */
	  int *k = mv[im];
	  /* 
	     Group /flav/interpolator/momentum
	  */
	  asprintf(&group_tag, "px%+dpy%+dpz%+d", k[0], k[1], k[2]);      
	  hid_t group2_id = H5Gcreate(group1_id, group_tag, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	  
	  /*
	    Create the (global) dataset
	  */
	  hsize_t dims[] = {d[0], 2};	  
	  hid_t filespace = H5Screate_simple(2, dims, NULL); 
	  hid_t dataset_id = H5Dcreate(group2_id, "corr_p", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	  /*
	    Subspace of this rank 
	  */
	  hsize_t ldims[] = {ld[0], 2};
	  hsize_t starts[] = {pc[0]*ld[0], 0};	
	  hid_t subspace = H5Screate_simple(2, ldims, NULL);
	  filespace = H5Dget_space(dataset_id);
	  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, starts, NULL, ldims, NULL);
	  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
	  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);      

	  for(int t=0; t<ld[0]; t++)	    
	    memcpy(&subdat[t*2], &buf[VGF((nm*t + im), interp, ifl)], sizeof(_Complex double));
	  
	  herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, subspace, filespace, plist_id, subdat);
	  H5Dclose(dataset_id);
	  H5Sclose(filespace);
	  H5Sclose(subspace);
	  H5Pclose(plist_id);	
	  H5Gclose(group2_id);
	}
	H5Gclose(group1_id);
      }
      H5Gclose(group_id);
    }
    H5Fclose(file_id);
    free(subdat);
  }
  free(buf);

  qhg_correlator_finalize(corr);
  return;
}
