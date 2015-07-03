#include <stdlib.h>

#include <complex.h>
#include <qhg_defs.h>
#include <qhg_idx.h>
#include <qhg_alloc.h>
#include <qhg_spinor_field.h>

#define G(v,mu) (((v)*ND + mu)*NC*NC)
#define S(v) ((v)*NS*NC)

double *
qhg_spinor_field_rms(qhg_spinor_field in, int origin[ND])
{
  int vol = in.lat->vol;
  int lvol = in.lat->lvol;
  int **nn = in.lat->nn;
  int *dims = in.lat->dims;
  int *ldims = in.lat->ldims;
  int *pc = in.lat->comms->proc_coords;
  MPI_Comm comm = in.lat->comms->comm;
  
  _Complex double *phi = in.field;  
  double *rms = qhg_alloc(dims[0]*sizeof(double));
  
  double numer[dims[0]];
  double denom[dims[0]];
  for(int i=0; i<dims[0]; i++) {
    numer[i] = 0;
    denom[i] = 0;
  }


#ifdef QHG_OMP
#pragma omp parallel shared(denom, numer)
  {
#endif

    double numer_priv[dims[0]];
    double denom_priv[dims[0]];
    for(int i=0; i<dims[0]; i++) {
      numer_priv[i] = 0;
      denom_priv[i] = 0;
    }

#ifdef QHG_OMP
#pragma omp for 
#endif
    for(int v=0; v<lvol; v++) {
      _Complex double *s = &phi[S(v)];
      int lx[ND] = CO(v, ldims);
      int x[ND];
      x[0] = lx[0] + pc[0]*ldims[0];
      for(int i=1; i<ND; i++)
	x[i] = (dims[i] + lx[i] + pc[i]*ldims[i] - origin[i]) % dims[i];

      for(int i=1; i<ND; i++)
	x[i] = x[i] > dims[i]/2 ? dims[i] - x[i] : x[i];
      
      

      int rsq = 0;
      for(int i=1; i<ND; i++)
	rsq += x[i]*x[i];

      double n = 0;
      for(int sp=0; sp<NC*NS; sp++)
	n += s[sp]*conj(s[sp]);

      numer_priv[x[0]] += n*rsq;
      denom_priv[x[0]] += n;
    }

#ifdef QHG_OMP
#pragma omp critical
    {
#endif  
      for(int i=0; i<dims[0]; i++) {
	numer[i] += numer_priv[i];
	denom[i] += denom_priv[i];
      }      
#ifdef QHG_OMP
    }
#endif  
    
#ifdef QHG_OMP
  }
#endif

  MPI_Allreduce(MPI_IN_PLACE, numer, dims[0], MPI_DOUBLE, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, denom, dims[0], MPI_DOUBLE, MPI_SUM, comm);  

  for(int i=0; i<dims[0]; i++)
    rms[i] = numer[i]/denom[i];
  
  return rms;
}
