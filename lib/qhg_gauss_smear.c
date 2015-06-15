#include <stdlib.h>

#include <complex.h>
#include <qhg_defs.h>
#include <qhg_spinor_linalg.h>
#include <qhg_spinor_field.h>
#include <qhg_gauge_field.h>
#include <qhg_xchange_gauge.h>
#include <qhg_xchange_spinor.h>

#define G(v,mu) (((v)*ND + mu)*NC*NC)
#define S(v) ((v)*NS*NC)
void
qhg_gauss_smear_iter(qhg_spinor_field out, qhg_spinor_field in, qhg_gauge_field gf, double alpha)
{
#ifdef QHG_OMP
#pragma omp parallel
  {
#endif

  double norm = 1.0/(1.0 + 6.0*alpha);
  int vol = in.lat->vol;
  int lvol = in.lat->lvol;
  int **nn = in.lat->nn;
    
#ifdef QHG_OMP
#pragma omp single
  {
#endif  

  qhg_xchange_spinor(in);

#ifdef QHG_OMP
  }
#endif  

  _Complex double *U = gf.field;
  _Complex double *phi = in.field;  
  _Complex double *psi = out.field;  
#ifdef QHG_OMP
#pragma omp for
#endif
  for(int v0=0; v0<lvol; v0++) {
    _Complex double *s = &psi[S(v0)];
    spinor_linalg_zero(s);
    for(int mu=1; mu<ND; mu++) {
      int vp = nn[mu][v0];
      int vm = nn[mu+ND][v0];
      _Complex double *p;
      _Complex double *u;

      p = &phi[S(vp)];
      u = &U[G(v0, mu)];
      spinor_linalg_ypequx(s, u, p);

      p = &phi[S(vm)];
      u = &U[G(vm, mu)];
      spinor_linalg_ypeqdx(s, u, p);
    }

    _Complex double *p = &phi[S(v0)];
    spinor_linalg_ax(alpha, s);
    spinor_linalg_ypeqx(s, p);
    spinor_linalg_ax(norm, s);
  }

#ifdef QHG_OMP
  }
#endif

  return;
}

void
qhg_gauss_smear(qhg_spinor_field out, qhg_spinor_field in, qhg_gauge_field gf,
		double alpha, int niter)
{
  qhg_spinor_field aux[2];
  aux[0] = qhg_spinor_field_init(in.lat);
  aux[1] = qhg_spinor_field_init(in.lat);

  qhg_xchange_gauge(gf);

  qhg_spinor_field_copy(aux[0], in);
  for(int i=0; i<niter; i++)
    qhg_gauss_smear_iter(aux[(i+1) % 2], aux[i % 2], gf, alpha);
  qhg_spinor_field_copy(out, aux[niter % 2]);
  
  qhg_spinor_field_finalize(aux[0]);
  qhg_spinor_field_finalize(aux[1]);  
  return;
}
