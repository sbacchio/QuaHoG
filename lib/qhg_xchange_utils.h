#ifndef _QHG_XCHANGE_UTILS_H
#define _QHG_XCHANGE_UTILS_H 1

/*
 * For a lattice where each site is site_size bytes, with ndims number
 * of dimensions, and with dimension lengths in dims[], create a slice
 * in direction dir. That is, create the slice whith dims[dir] = 1;
 * and return the MPI_Datatype that corrsponds to it. slice is the
 * coordinate in direction dir of the slice. That is, the MPI_Datatype
 * will return the slice: (0:dims[0]-1, 0:dims[1]-1, ..., slice,
 *                         0:dims[dir+1]-1, ..., 0:dims[ndims-1]-1)
 */
static MPI_Datatype
get_slice(int ndims, int dims[], int dir, int slice, size_t site_size)
{
  MPI_Datatype elem;
  MPI_Type_contiguous(site_size, MPI_BYTE, &elem);
  MPI_Type_commit(&elem);

  unsigned long int vol = 1;
  for(int i=0; i<ND; i++)
    vol *= (unsigned long int)dims[i];
  
  unsigned long int nelems = 1;
  for(int i=ndims-1; i>dir; i--)
    nelems *= (unsigned long int)dims[i];

  int count = (int)(vol/nelems/(unsigned long int)dims[dir]);
  
  unsigned long int stride = (unsigned long int)dims[dir]*nelems;
  unsigned long int offset = (unsigned long int)slice*nelems;
  int *nelem_arr = qhg_alloc(sizeof(int)*count);
  int *displ_arr = qhg_alloc(sizeof(int)*count);
  for(int i=0; i<count; i++) {
    nelem_arr[i] = (int)nelems;
    displ_arr[i] = (int)(offset + i*stride);
  }

  MPI_Datatype slice_type;
  MPI_Type_indexed(count, nelem_arr, displ_arr, elem, &slice_type);  
  MPI_Type_commit(&slice_type);
  MPI_Type_free(&elem);
  
  free(nelem_arr);
  free(displ_arr);
  return slice_type;
}

#endif /* _QHG_XCHANGE_UTILS_H */
