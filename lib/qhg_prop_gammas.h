#ifndef _QHG_PROP_GAMMAS_H
#define _QHG_PROP_GAMMAS_H 1

#include <qhg_prop_gammas_decl.h>

/* multiply prop by 1 from the left */
static inline void
prop_1_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_1_G_nva; k++)
           out[i][j] += (prop_1_G_val[i][j][k])*in[prop_1_G_idx[i][j][k][0]][prop_1_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{1} from the right */
static inline void
prop_G_1(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_1_nva; k++)
           out[i][j] += (prop_G_1_val[i][j][k])*in[prop_G_1_idx[i][j][k][0]][prop_G_1_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by \gamma_0 from the left */
static inline void
prop_g0_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_g0_G_nva; k++)
           out[i][j] += (prop_g0_G_val[i][j][k])*in[prop_g0_G_idx[i][j][k][0]][prop_g0_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{\gamma_0} from the right */
static inline void
prop_G_g0(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_g0_nva; k++)
           out[i][j] += (prop_G_g0_val[i][j][k])*in[prop_G_g0_idx[i][j][k][0]][prop_G_g0_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by \gamma_x from the left */
static inline void
prop_gx_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_gx_G_nva; k++)
           out[i][j] += (prop_gx_G_val[i][j][k])*in[prop_gx_G_idx[i][j][k][0]][prop_gx_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{\gamma_x} from the right */
static inline void
prop_G_gx(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_gx_nva; k++)
           out[i][j] += (prop_G_gx_val[i][j][k])*in[prop_G_gx_idx[i][j][k][0]][prop_G_gx_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by \gamma_y from the left */
static inline void
prop_gy_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_gy_G_nva; k++)
           out[i][j] += (prop_gy_G_val[i][j][k])*in[prop_gy_G_idx[i][j][k][0]][prop_gy_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{\gamma_y} from the right */
static inline void
prop_G_gy(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_gy_nva; k++)
           out[i][j] += (prop_G_gy_val[i][j][k])*in[prop_G_gy_idx[i][j][k][0]][prop_G_gy_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by \gamma_z from the left */
static inline void
prop_gz_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_gz_G_nva; k++)
           out[i][j] += (prop_gz_G_val[i][j][k])*in[prop_gz_G_idx[i][j][k][0]][prop_gz_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{\gamma_z} from the right */
static inline void
prop_G_gz(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_gz_nva; k++)
           out[i][j] += (prop_G_gz_val[i][j][k])*in[prop_G_gz_idx[i][j][k][0]][prop_G_gz_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by \gamma_t from the left */
static inline void
prop_gt_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_gt_G_nva; k++)
           out[i][j] += (prop_gt_G_val[i][j][k])*in[prop_gt_G_idx[i][j][k][0]][prop_gt_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{\gamma_t} from the right */
static inline void
prop_G_gt(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_gt_nva; k++)
           out[i][j] += (prop_G_gt_val[i][j][k])*in[prop_G_gt_idx[i][j][k][0]][prop_G_gt_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by \gamma_5 from the left */
static inline void
prop_g5_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_g5_G_nva; k++)
           out[i][j] += (prop_g5_G_val[i][j][k])*in[prop_g5_G_idx[i][j][k][0]][prop_g5_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{\gamma_5} from the right */
static inline void
prop_G_g5(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_g5_nva; k++)
           out[i][j] += (prop_G_g5_val[i][j][k])*in[prop_G_g5_idx[i][j][k][0]][prop_G_g5_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by \gamma_5\gamma_0 from the left */
static inline void
prop_g5g0_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_g5g0_G_nva; k++)
           out[i][j] += (prop_g5g0_G_val[i][j][k])*in[prop_g5g0_G_idx[i][j][k][0]][prop_g5g0_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{\gamma_5\gamma_0} from the right */
static inline void
prop_G_g5g0(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_g5g0_nva; k++)
           out[i][j] += (prop_G_g5g0_val[i][j][k])*in[prop_G_g5g0_idx[i][j][k][0]][prop_G_g5g0_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by \gamma_5\gamma_x from the left */
static inline void
prop_g5gx_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_g5gx_G_nva; k++)
           out[i][j] += (prop_g5gx_G_val[i][j][k])*in[prop_g5gx_G_idx[i][j][k][0]][prop_g5gx_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{\gamma_5\gamma_x} from the right */
static inline void
prop_G_g5gx(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_g5gx_nva; k++)
           out[i][j] += (prop_G_g5gx_val[i][j][k])*in[prop_G_g5gx_idx[i][j][k][0]][prop_G_g5gx_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by \gamma_5\gamma_y from the left */
static inline void
prop_g5gy_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_g5gy_G_nva; k++)
           out[i][j] += (prop_g5gy_G_val[i][j][k])*in[prop_g5gy_G_idx[i][j][k][0]][prop_g5gy_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{\gamma_5\gamma_y} from the right */
static inline void
prop_G_g5gy(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_g5gy_nva; k++)
           out[i][j] += (prop_G_g5gy_val[i][j][k])*in[prop_G_g5gy_idx[i][j][k][0]][prop_G_g5gy_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by \gamma_5\gamma_z from the left */
static inline void
prop_g5gz_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_g5gz_G_nva; k++)
           out[i][j] += (prop_g5gz_G_val[i][j][k])*in[prop_g5gz_G_idx[i][j][k][0]][prop_g5gz_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{\gamma_5\gamma_z} from the right */
static inline void
prop_G_g5gz(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_g5gz_nva; k++)
           out[i][j] += (prop_G_g5gz_val[i][j][k])*in[prop_G_g5gz_idx[i][j][k][0]][prop_G_g5gz_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by \gamma_5\gamma_t from the left */
static inline void
prop_g5gt_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_g5gt_G_nva; k++)
           out[i][j] += (prop_g5gt_G_val[i][j][k])*in[prop_g5gt_G_idx[i][j][k][0]][prop_g5gt_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{\gamma_5\gamma_t} from the right */
static inline void
prop_G_g5gt(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_g5gt_nva; k++)
           out[i][j] += (prop_G_g5gt_val[i][j][k])*in[prop_G_g5gt_idx[i][j][k][0]][prop_G_g5gt_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by C from the left */
static inline void
prop_C_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_C_G_nva; k++)
           out[i][j] += (prop_C_G_val[i][j][k])*in[prop_C_G_idx[i][j][k][0]][prop_C_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{C} from the right */
static inline void
prop_G_C(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_C_nva; k++)
           out[i][j] += (prop_G_C_val[i][j][k])*in[prop_G_C_idx[i][j][k][0]][prop_G_C_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by C\gamma_5 from the left */
static inline void
prop_Cg5_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_Cg5_G_nva; k++)
           out[i][j] += (prop_Cg5_G_val[i][j][k])*in[prop_Cg5_G_idx[i][j][k][0]][prop_Cg5_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{C\gamma_5} from the right */
static inline void
prop_G_Cg5(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_Cg5_nva; k++)
           out[i][j] += (prop_G_Cg5_val[i][j][k])*in[prop_G_Cg5_idx[i][j][k][0]][prop_G_Cg5_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by C\gamma_x from the left */
static inline void
prop_Cgx_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_Cgx_G_nva; k++)
           out[i][j] += (prop_Cgx_G_val[i][j][k])*in[prop_Cgx_G_idx[i][j][k][0]][prop_Cgx_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{C\gamma_x} from the right */
static inline void
prop_G_Cgx(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_Cgx_nva; k++)
           out[i][j] += (prop_G_Cgx_val[i][j][k])*in[prop_G_Cgx_idx[i][j][k][0]][prop_G_Cgx_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by C\gamma_y from the left */
static inline void
prop_Cgy_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_Cgy_G_nva; k++)
           out[i][j] += (prop_Cgy_G_val[i][j][k])*in[prop_Cgy_G_idx[i][j][k][0]][prop_Cgy_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{C\gamma_y} from the right */
static inline void
prop_G_Cgy(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_Cgy_nva; k++)
           out[i][j] += (prop_G_Cgy_val[i][j][k])*in[prop_G_Cgy_idx[i][j][k][0]][prop_G_Cgy_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by C\gamma_z from the left */
static inline void
prop_Cgz_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_Cgz_G_nva; k++)
           out[i][j] += (prop_Cgz_G_val[i][j][k])*in[prop_Cgz_G_idx[i][j][k][0]][prop_Cgz_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{C\gamma_z} from the right */
static inline void
prop_G_Cgz(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_Cgz_nva; k++)
           out[i][j] += (prop_G_Cgz_val[i][j][k])*in[prop_G_Cgz_idx[i][j][k][0]][prop_G_Cgz_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by C\gamma_t from the left */
static inline void
prop_Cgt_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_Cgt_G_nva; k++)
           out[i][j] += (prop_Cgt_G_val[i][j][k])*in[prop_Cgt_G_idx[i][j][k][0]][prop_Cgt_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{C\gamma_t} from the right */
static inline void
prop_G_Cgt(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_Cgt_nva; k++)
           out[i][j] += (prop_G_Cgt_val[i][j][k])*in[prop_G_Cgt_idx[i][j][k][0]][prop_G_Cgt_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by \gamma_5\sigma_{t,x} = -0.5\gamma_5[\gamma_t, \gamma_x] from the left */
static inline void
prop_g5sitx_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_g5sitx_G_nva; k++)
           out[i][j] += (prop_g5sitx_G_val[i][j][k])*in[prop_g5sitx_G_idx[i][j][k][0]][prop_g5sitx_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{\gamma_5\sigma_{t,x} = -0.5\gamma_5[\gamma_t, \gamma_x]} from the right */
static inline void
prop_G_g5sitx(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_g5sitx_nva; k++)
           out[i][j] += (prop_G_g5sitx_val[i][j][k])*in[prop_G_g5sitx_idx[i][j][k][0]][prop_G_g5sitx_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by \gamma_5\sigma_{t,y} = -0.5\gamma_5[\gamma_t, \gamma_y] from the left */
static inline void
prop_g5sity_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_g5sity_G_nva; k++)
           out[i][j] += (prop_g5sity_G_val[i][j][k])*in[prop_g5sity_G_idx[i][j][k][0]][prop_g5sity_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{\gamma_5\sigma_{t,y} = -0.5\gamma_5[\gamma_t, \gamma_y]} from the right */
static inline void
prop_G_g5sity(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_g5sity_nva; k++)
           out[i][j] += (prop_G_g5sity_val[i][j][k])*in[prop_G_g5sity_idx[i][j][k][0]][prop_G_g5sity_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by \gamma_5\sigma_{t,z} = -0.5\gamma_5[\gamma_t, \gamma_z] from the left */
static inline void
prop_g5sitz_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_g5sitz_G_nva; k++)
           out[i][j] += (prop_g5sitz_G_val[i][j][k])*in[prop_g5sitz_G_idx[i][j][k][0]][prop_g5sitz_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{\gamma_5\sigma_{t,z} = -0.5\gamma_5[\gamma_t, \gamma_z]} from the right */
static inline void
prop_G_g5sitz(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_g5sitz_nva; k++)
           out[i][j] += (prop_G_g5sitz_val[i][j][k])*in[prop_G_g5sitz_idx[i][j][k][0]][prop_G_g5sitz_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by \gamma_5\sigma_{0,x} = -0.5\gamma_5[\gamma_0, \gamma_x] from the left */
static inline void
prop_g5si0x_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_g5si0x_G_nva; k++)
           out[i][j] += (prop_g5si0x_G_val[i][j][k])*in[prop_g5si0x_G_idx[i][j][k][0]][prop_g5si0x_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{\gamma_5\sigma_{0,x} = -0.5\gamma_5[\gamma_0, \gamma_x]} from the right */
static inline void
prop_G_g5si0x(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_g5si0x_nva; k++)
           out[i][j] += (prop_G_g5si0x_val[i][j][k])*in[prop_G_g5si0x_idx[i][j][k][0]][prop_G_g5si0x_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by \gamma_5\sigma_{0,y} = -0.5\gamma_5[\gamma_0, \gamma_y] from the left */
static inline void
prop_g5si0y_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_g5si0y_G_nva; k++)
           out[i][j] += (prop_g5si0y_G_val[i][j][k])*in[prop_g5si0y_G_idx[i][j][k][0]][prop_g5si0y_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{\gamma_5\sigma_{0,y} = -0.5\gamma_5[\gamma_0, \gamma_y]} from the right */
static inline void
prop_G_g5si0y(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_g5si0y_nva; k++)
           out[i][j] += (prop_G_g5si0y_val[i][j][k])*in[prop_G_g5si0y_idx[i][j][k][0]][prop_G_g5si0y_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by \gamma_5\sigma_{0,z} = -0.5\gamma_5[\gamma_0, \gamma_z] from the left */
static inline void
prop_g5si0z_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_g5si0z_G_nva; k++)
           out[i][j] += (prop_g5si0z_G_val[i][j][k])*in[prop_g5si0z_G_idx[i][j][k][0]][prop_g5si0z_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{\gamma_5\sigma_{0,z} = -0.5\gamma_5[\gamma_0, \gamma_z]} from the right */
static inline void
prop_G_g5si0z(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_g5si0z_nva; k++)
           out[i][j] += (prop_G_g5si0z_val[i][j][k])*in[prop_G_g5si0z_idx[i][j][k][0]][prop_G_g5si0z_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by \gamma_5\sigma_{x,y} = -0.5\gamma_5[\gamma_x, \gamma_y] from the left */
static inline void
prop_g5sixy_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_g5sixy_G_nva; k++)
           out[i][j] += (prop_g5sixy_G_val[i][j][k])*in[prop_g5sixy_G_idx[i][j][k][0]][prop_g5sixy_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{\gamma_5\sigma_{x,y} = -0.5\gamma_5[\gamma_x, \gamma_y]} from the right */
static inline void
prop_G_g5sixy(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_g5sixy_nva; k++)
           out[i][j] += (prop_G_g5sixy_val[i][j][k])*in[prop_G_g5sixy_idx[i][j][k][0]][prop_G_g5sixy_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by \gamma_5\sigma_{x,z} = -0.5\gamma_5[\gamma_x, \gamma_z] from the left */
static inline void
prop_g5sixz_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_g5sixz_G_nva; k++)
           out[i][j] += (prop_g5sixz_G_val[i][j][k])*in[prop_g5sixz_G_idx[i][j][k][0]][prop_g5sixz_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{\gamma_5\sigma_{x,z} = -0.5\gamma_5[\gamma_x, \gamma_z]} from the right */
static inline void
prop_G_g5sixz(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_g5sixz_nva; k++)
           out[i][j] += (prop_G_g5sixz_val[i][j][k])*in[prop_G_g5sixz_idx[i][j][k][0]][prop_G_g5sixz_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by \gamma_5\sigma_{y,z} = -0.5\gamma_5[\gamma_y, \gamma_z] from the left */
static inline void
prop_g5siyz_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_g5siyz_G_nva; k++)
           out[i][j] += (prop_g5siyz_G_val[i][j][k])*in[prop_g5siyz_G_idx[i][j][k][0]][prop_g5siyz_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{\gamma_5\sigma_{y,z} = -0.5\gamma_5[\gamma_y, \gamma_z]} from the right */
static inline void
prop_G_g5siyz(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_g5siyz_nva; k++)
           out[i][j] += (prop_G_g5siyz_val[i][j][k])*in[prop_G_g5siyz_idx[i][j][k][0]][prop_G_g5siyz_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by \gamma_5\sigma_{x,t} = -0.5\gamma_5[\gamma_x, \gamma_t] from the left */
static inline void
prop_g5sixt_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_g5sixt_G_nva; k++)
           out[i][j] += (prop_g5sixt_G_val[i][j][k])*in[prop_g5sixt_G_idx[i][j][k][0]][prop_g5sixt_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{\gamma_5\sigma_{x,t} = -0.5\gamma_5[\gamma_x, \gamma_t]} from the right */
static inline void
prop_G_g5sixt(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_g5sixt_nva; k++)
           out[i][j] += (prop_G_g5sixt_val[i][j][k])*in[prop_G_g5sixt_idx[i][j][k][0]][prop_G_g5sixt_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by \gamma_5\sigma_{y,t} = -0.5\gamma_5[\gamma_y, \gamma_t] from the left */
static inline void
prop_g5siyt_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_g5siyt_G_nva; k++)
           out[i][j] += (prop_g5siyt_G_val[i][j][k])*in[prop_g5siyt_G_idx[i][j][k][0]][prop_g5siyt_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{\gamma_5\sigma_{y,t} = -0.5\gamma_5[\gamma_y, \gamma_t]} from the right */
static inline void
prop_G_g5siyt(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_g5siyt_nva; k++)
           out[i][j] += (prop_G_g5siyt_val[i][j][k])*in[prop_G_g5siyt_idx[i][j][k][0]][prop_G_g5siyt_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by \gamma_5\sigma_{z,t} = -0.5\gamma_5[\gamma_z, \gamma_t] from the left */
static inline void
prop_g5sizt_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_g5sizt_G_nva; k++)
           out[i][j] += (prop_g5sizt_G_val[i][j][k])*in[prop_g5sizt_G_idx[i][j][k][0]][prop_g5sizt_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{\gamma_5\sigma_{z,t} = -0.5\gamma_5[\gamma_z, \gamma_t]} from the right */
static inline void
prop_G_g5sizt(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_g5sizt_nva; k++)
           out[i][j] += (prop_G_g5sizt_val[i][j][k])*in[prop_G_g5sizt_idx[i][j][k][0]][prop_G_g5sizt_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by \gamma_5\sigma_{x,0} = -0.5\gamma_5[\gamma_x, \gamma_0] from the left */
static inline void
prop_g5six0_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_g5six0_G_nva; k++)
           out[i][j] += (prop_g5six0_G_val[i][j][k])*in[prop_g5six0_G_idx[i][j][k][0]][prop_g5six0_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{\gamma_5\sigma_{x,0} = -0.5\gamma_5[\gamma_x, \gamma_0]} from the right */
static inline void
prop_G_g5six0(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_g5six0_nva; k++)
           out[i][j] += (prop_G_g5six0_val[i][j][k])*in[prop_G_g5six0_idx[i][j][k][0]][prop_G_g5six0_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by \gamma_5\sigma_{y,0} = -0.5\gamma_5[\gamma_y, \gamma_0] from the left */
static inline void
prop_g5siy0_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_g5siy0_G_nva; k++)
           out[i][j] += (prop_g5siy0_G_val[i][j][k])*in[prop_g5siy0_G_idx[i][j][k][0]][prop_g5siy0_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{\gamma_5\sigma_{y,0} = -0.5\gamma_5[\gamma_y, \gamma_0]} from the right */
static inline void
prop_G_g5siy0(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_g5siy0_nva; k++)
           out[i][j] += (prop_G_g5siy0_val[i][j][k])*in[prop_G_g5siy0_idx[i][j][k][0]][prop_G_g5siy0_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by \gamma_5\sigma_{z,0} = -0.5\gamma_5[\gamma_z, \gamma_0] from the left */
static inline void
prop_g5siz0_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_g5siz0_G_nva; k++)
           out[i][j] += (prop_g5siz0_G_val[i][j][k])*in[prop_g5siz0_G_idx[i][j][k][0]][prop_g5siz0_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{\gamma_5\sigma_{z,0} = -0.5\gamma_5[\gamma_z, \gamma_0]} from the right */
static inline void
prop_G_g5siz0(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_g5siz0_nva; k++)
           out[i][j] += (prop_G_g5siz0_val[i][j][k])*in[prop_G_g5siz0_idx[i][j][k][0]][prop_G_g5siz0_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by \gamma_5\sigma_{y,x} = -0.5\gamma_5[\gamma_y, \gamma_x] from the left */
static inline void
prop_g5siyx_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_g5siyx_G_nva; k++)
           out[i][j] += (prop_g5siyx_G_val[i][j][k])*in[prop_g5siyx_G_idx[i][j][k][0]][prop_g5siyx_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{\gamma_5\sigma_{y,x} = -0.5\gamma_5[\gamma_y, \gamma_x]} from the right */
static inline void
prop_G_g5siyx(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_g5siyx_nva; k++)
           out[i][j] += (prop_G_g5siyx_val[i][j][k])*in[prop_G_g5siyx_idx[i][j][k][0]][prop_G_g5siyx_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by \gamma_5\sigma_{z,x} = -0.5\gamma_5[\gamma_z, \gamma_x] from the left */
static inline void
prop_g5sizx_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_g5sizx_G_nva; k++)
           out[i][j] += (prop_g5sizx_G_val[i][j][k])*in[prop_g5sizx_G_idx[i][j][k][0]][prop_g5sizx_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{\gamma_5\sigma_{z,x} = -0.5\gamma_5[\gamma_z, \gamma_x]} from the right */
static inline void
prop_G_g5sizx(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_g5sizx_nva; k++)
           out[i][j] += (prop_G_g5sizx_val[i][j][k])*in[prop_G_g5sizx_idx[i][j][k][0]][prop_G_g5sizx_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by \gamma_5\sigma_{z,y} = -0.5\gamma_5[\gamma_z, \gamma_y] from the left */
static inline void
prop_g5sizy_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_g5sizy_G_nva; k++)
           out[i][j] += (prop_g5sizy_G_val[i][j][k])*in[prop_g5sizy_G_idx[i][j][k][0]][prop_g5sizy_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{\gamma_5\sigma_{z,y} = -0.5\gamma_5[\gamma_z, \gamma_y]} from the right */
static inline void
prop_G_g5sizy(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_g5sizy_nva; k++)
           out[i][j] += (prop_G_g5sizy_val[i][j][k])*in[prop_G_g5sizy_idx[i][j][k][0]][prop_G_g5sizy_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by (1+\gamma_0) from the left */
static inline void
prop_1pg0_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_1pg0_G_nva; k++)
           out[i][j] += (prop_1pg0_G_val[i][j][k])*in[prop_1pg0_G_idx[i][j][k][0]][prop_1pg0_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{(1+\gamma_0)} from the right */
static inline void
prop_G_1pg0(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_1pg0_nva; k++)
           out[i][j] += (prop_G_1pg0_val[i][j][k])*in[prop_G_1pg0_idx[i][j][k][0]][prop_G_1pg0_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by (1+\gamma_x) from the left */
static inline void
prop_1pgx_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_1pgx_G_nva; k++)
           out[i][j] += (prop_1pgx_G_val[i][j][k])*in[prop_1pgx_G_idx[i][j][k][0]][prop_1pgx_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{(1+\gamma_x)} from the right */
static inline void
prop_G_1pgx(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_1pgx_nva; k++)
           out[i][j] += (prop_G_1pgx_val[i][j][k])*in[prop_G_1pgx_idx[i][j][k][0]][prop_G_1pgx_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by (1+\gamma_y) from the left */
static inline void
prop_1pgy_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_1pgy_G_nva; k++)
           out[i][j] += (prop_1pgy_G_val[i][j][k])*in[prop_1pgy_G_idx[i][j][k][0]][prop_1pgy_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{(1+\gamma_y)} from the right */
static inline void
prop_G_1pgy(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_1pgy_nva; k++)
           out[i][j] += (prop_G_1pgy_val[i][j][k])*in[prop_G_1pgy_idx[i][j][k][0]][prop_G_1pgy_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by (1+\gamma_z) from the left */
static inline void
prop_1pgz_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_1pgz_G_nva; k++)
           out[i][j] += (prop_1pgz_G_val[i][j][k])*in[prop_1pgz_G_idx[i][j][k][0]][prop_1pgz_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{(1+\gamma_z)} from the right */
static inline void
prop_G_1pgz(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_1pgz_nva; k++)
           out[i][j] += (prop_G_1pgz_val[i][j][k])*in[prop_G_1pgz_idx[i][j][k][0]][prop_G_1pgz_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by (1+\gamma_t) from the left */
static inline void
prop_1pgt_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_1pgt_G_nva; k++)
           out[i][j] += (prop_1pgt_G_val[i][j][k])*in[prop_1pgt_G_idx[i][j][k][0]][prop_1pgt_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{(1+\gamma_t)} from the right */
static inline void
prop_G_1pgt(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_1pgt_nva; k++)
           out[i][j] += (prop_G_1pgt_val[i][j][k])*in[prop_G_1pgt_idx[i][j][k][0]][prop_G_1pgt_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by (1-\gamma_0) from the left */
static inline void
prop_1mg0_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_1mg0_G_nva; k++)
           out[i][j] += (prop_1mg0_G_val[i][j][k])*in[prop_1mg0_G_idx[i][j][k][0]][prop_1mg0_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{(1-\gamma_0)} from the right */
static inline void
prop_G_1mg0(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_1mg0_nva; k++)
           out[i][j] += (prop_G_1mg0_val[i][j][k])*in[prop_G_1mg0_idx[i][j][k][0]][prop_G_1mg0_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by (1-\gamma_x) from the left */
static inline void
prop_1mgx_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_1mgx_G_nva; k++)
           out[i][j] += (prop_1mgx_G_val[i][j][k])*in[prop_1mgx_G_idx[i][j][k][0]][prop_1mgx_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{(1-\gamma_x)} from the right */
static inline void
prop_G_1mgx(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_1mgx_nva; k++)
           out[i][j] += (prop_G_1mgx_val[i][j][k])*in[prop_G_1mgx_idx[i][j][k][0]][prop_G_1mgx_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by (1-\gamma_y) from the left */
static inline void
prop_1mgy_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_1mgy_G_nva; k++)
           out[i][j] += (prop_1mgy_G_val[i][j][k])*in[prop_1mgy_G_idx[i][j][k][0]][prop_1mgy_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{(1-\gamma_y)} from the right */
static inline void
prop_G_1mgy(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_1mgy_nva; k++)
           out[i][j] += (prop_G_1mgy_val[i][j][k])*in[prop_G_1mgy_idx[i][j][k][0]][prop_G_1mgy_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by (1-\gamma_z) from the left */
static inline void
prop_1mgz_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_1mgz_G_nva; k++)
           out[i][j] += (prop_1mgz_G_val[i][j][k])*in[prop_1mgz_G_idx[i][j][k][0]][prop_1mgz_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{(1-\gamma_z)} from the right */
static inline void
prop_G_1mgz(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_1mgz_nva; k++)
           out[i][j] += (prop_G_1mgz_val[i][j][k])*in[prop_G_1mgz_idx[i][j][k][0]][prop_G_1mgz_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by (1-\gamma_t) from the left */
static inline void
prop_1mgt_G(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_1mgt_G_nva; k++)
           out[i][j] += (prop_1mgt_G_val[i][j][k])*in[prop_1mgt_G_idx[i][j][k][0]][prop_1mgt_G_idx[i][j][k][1]];
     }

  return;
}

/* multiply prop by ar{(1-\gamma_t)} from the right */
static inline void
prop_G_1mgt(_Complex double out[NS*NC][NS*NC], _Complex double in[NS*NC][NS*NC])
{
  for(int i=0; i<NC*NS; i++)
     for(int j=0; j<NC*NS; j++) {
        out[i][j] = 0.;
        for(int k=0; k<prop_G_1mgt_nva; k++)
           out[i][j] += (prop_G_1mgt_val[i][j][k])*in[prop_G_1mgt_idx[i][j][k][0]][prop_G_1mgt_idx[i][j][k][1]];
     }

  return;
}



#endif /* _QHG_PROP_GAMMAS_H */
