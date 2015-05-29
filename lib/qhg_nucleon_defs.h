#ifndef _QHG_NUCLEON_DEFS_H
#define _QHG_NUCLEON_DEFS_H 1

#define NFLAV 2
#define NCOMP (4)
#define NCHAN (4)
#define SITE_SIZE (NFLAV*NCOMP*NCOMP*NCHAN)
#define NIDX(v, f, chi, si0, si1) (si1 + NCOMP*(si0 + NCOMP*(chi + NCHAN*(f+NFLAV*v))))

static char flav_tags[NFLAV][256] = {"ppm\0","pmm\0"};
static char chan_tags[NCHAN][256] = {"1-1\0","1-2\0","2-1\0","2-2\0"};

#endif /* _QHG_NUCLEON_DEFS_H */
