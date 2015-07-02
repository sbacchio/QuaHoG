#ifndef _QHG_MESON_DEFS_H
#define _QHG_MESON_DEFS_H 1

#define NGAMMAS 10			/* 1,gx,...,gt,g5,gxg5,...,gtg5 */
#define NFLAV 2				/* up/down */
#define NCHANNELS ((NGAMMAS)*(NFLAV))
#define VGF(v, g, f) ((v)*NCHANNELS + (g)*NFLAV + (f))

static char gamma_tags[NGAMMAS][256] = {"=1=\0",
					"=g5=\0",
					"=gx=\0",
					"=gy=\0",
					"=gz=\0",
					"=gt=\0",
					"=g5gx=\0",
					"=g5gy=\0",
					"=g5gz=\0",
					"=g5gt=\0"};

static char flav_tags[NFLAV][256] = {"up\0", "dn\0"};

#endif /* _QHG_MESON_DEFS_H */
