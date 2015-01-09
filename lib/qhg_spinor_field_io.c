#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <mpi.h>

#include <qhg_types.h>
#include <qhg_defs.h>
#include <qhg_alloc.h>
#include <qhg_io_utils.h>

static void
get_file_types(MPI_Datatype *etype, MPI_Datatype *ftype, qhg_lattice *lat)
{
  MPI_Type_contiguous(2*NC*NS, MPI_FLOAT, etype);
  MPI_Type_commit(etype);
  
  int *ldims = lat->ldims;
  int *dims = lat->dims;  
  int *proc_coords = lat->comms->proc_coords;
  int starts[ND];
  for(int i=0; i<ND; i++)
    starts[i] = proc_coords[i]*ldims[i];

  int ierr = MPI_Type_create_subarray(ND, dims, ldims, starts, MPI_ORDER_C, *etype, ftype);
  MPI_Type_commit(ftype);
  return;
}
 
void
qhg_read_spinors(qhg_spinor_field psi[], int n_spinors, char fname[])
{
  qhg_lattice *lat = psi[0].lat;
  int am_io_proc = lat->comms->proc_id == 0 ? 1 : 0;

  MPI_Datatype elemtype, filetype;
  get_file_types(&elemtype, &filetype, lat);
  
  MPI_File fhandle;
  MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fhandle);
  MPI_File_set_size(fhandle, 0); 
  MPI_File_set_view(fhandle, 0, elemtype, filetype, "native", MPI_INFO_NULL);
  
  int lvol = lat->lvol;
  float *buffer = qhg_alloc(lvol*2*sizeof(float)*NC*NS);
  for(int i_sp=0; i_sp<n_spinors; i_sp++) {
    MPI_File_read_all(fhandle, buffer, lvol, elemtype, MPI_STATUS_IGNORE);
    if(!qhg_is_bigendian())
      qhg_byte_swap_float(buffer, lvol*NC*NS*2);
    
    for(int v=0; v<lvol; v++)
      for(int sp=0; sp<NC*NS; sp++) {
	int x = (v*NC*NS + sp);
	double re = buffer[x*2 + 0];
	double im = buffer[x*2 + 1];
	psi[i_sp].field[x] = re + _Complex_I*im;
      }    
  }
  free(buffer);    

  MPI_File_close(&fhandle);
  MPI_Type_free(&elemtype);
  MPI_Type_free(&filetype);
  
  return;
}

void
qhg_write_spinors(char fname[], int n_spinors, qhg_spinor_field psi[])
{
  qhg_lattice *lat = psi[0].lat;
  int am_io_proc = lat->comms->proc_id == 0 ? 1 : 0;

  MPI_Datatype elemtype, filetype;
  get_file_types(&elemtype, &filetype, lat);
  
  MPI_File fhandle;
  MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fhandle);
  MPI_File_set_size(fhandle, 0); 
  MPI_File_set_view(fhandle, 0, elemtype, filetype, "native", MPI_INFO_NULL);

  int lvol = lat->lvol;
  float *buffer = qhg_alloc(lvol*2*sizeof(float)*NC*NS);
  for(int i_sp=0; i_sp<n_spinors; i_sp++) {
    for(int v=0; v<lvol; v++)
      for(int sp=0; sp<NC*NS; sp++) {
	int x = (v*NC*NS + sp);
	buffer[x*2 + 0] = creal(psi[i_sp].field[x]);
	buffer[x*2 + 1] = cimag(psi[i_sp].field[x]);
      }

    if(!qhg_is_bigendian())
      qhg_byte_swap_float(buffer, lvol*NC*NS*2);
    
    MPI_File_write_all(fhandle, buffer, lvol, elemtype, MPI_STATUS_IGNORE);
  }
  free(buffer);    

  MPI_File_close(&fhandle);
  MPI_Type_free(&elemtype);
  MPI_Type_free(&filetype);
  
  return;
}
