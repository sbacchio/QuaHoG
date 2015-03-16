#include <stdlib.h>

#include <complex.h>
#include <qhg_defs.h>
#include <qhg_su3_ops.h>
#include <qhg_su3_project.h>
#include <qhg_gauge_field.h>
#include <qhg_xchange_gauge.h>

#define G(v,mu) (((v)*ND + mu)*NC*NC)

void
qhg_ape_smear_3d_iter(qhg_gauge_field out, qhg_gauge_field in, double alpha_in)
{
  int vol = in.lat->vol;
  int lvol = in.lat->lvol;
  int **nn = in.lat->nn;
  double alpha = 4.0*alpha_in / (1.0 + 4.0*alpha_in);

  qhg_xchange_gauge(in);

  _Complex double *U = in.field;
  _Complex double *V = out.field;  
  _Complex double *u0, *u1, *u2;
  _Complex double staple[NC*NC], w[NC*NC], u[NC*NC];
  _Complex double one_minus_alpha[NC*NC];
  _Complex double z[NC] = {1-alpha, 1-alpha, 1-alpha};

  su3_linalg_diag(one_minus_alpha, z);
  
  for(int mu=1; mu<ND; mu++)
    for(int v00=0; v00<lvol; v00++) {
      su3_linalg_zero(staple);

      /* create the staple term for this mu direction */
      for(int nu=(mu==1?2:1); nu<ND; nu==mu-1?nu+=2:nu++)	{
      	int vp0 = nn[mu][v00];
      	int v0p = nn[nu][v00];
      	int v0m = nn[nu+ND][v00];	
      	int vpm = nn[nu+ND][nn[mu][v00]];
	
      	/* Fwd staple */
      	u0 = &U[G(v00, nu)];
      	u1 = &U[G(v0p, mu)];
      	u2 = &U[G(vp0, nu)];
      	su3_mul_uu(u, u0, u1);
      	su3_mul_ud(w, u, u2);
	
      	su3_linalg_upeqv(staple, w);

      	/* Bwd staple */
      	u0 = &U[G(v0m, nu)];
      	u1 = &U[G(v0m, mu)];
      	u2 = &U[G(vpm, nu)];
      	su3_mul_du(u, u0, u1);
      	su3_mul_uu(w, u, u2);
	
      	su3_linalg_upeqv(staple, w);
      }
      u0 = &U[G(v00, mu)];
      u1 = &V[G(v00, mu)];
      /* multiply staple with u^+*alpha/4 and add to 1-alpha */

      su3_mul_ud(w, staple, u0);
      su3_linalg_au(alpha*0.25, w);
      su3_linalg_upeqv(w, one_minus_alpha);

      qhg_su3_project(w, 1);
      
      su3_mul_uu(u1, w, u0);
    }
  
  return;
}

void
qhg_ape_smear_3d(qhg_gauge_field out, qhg_gauge_field in, double alpha, int niter)
{
  qhg_gauge_field aux[2] = {qhg_gauge_field_init(in.lat),
			    qhg_gauge_field_init(in.lat)};

  qhg_gauge_field_copy(aux[0], in);
  qhg_gauge_field_copy(aux[1], in); /* This will copy the temporal links which will be left untouched */
  for(int i=0; i<niter; i++)
    qhg_ape_smear_3d_iter(aux[(i+1) % 2], aux[i % 2], alpha);
  qhg_gauge_field_copy(out, aux[niter % 2]);
  
  qhg_gauge_field_finalize(aux[0]);
  qhg_gauge_field_finalize(aux[1]);  
  return;
}
