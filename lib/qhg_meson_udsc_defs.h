#ifndef _QHG_MESON_UDSC_DEFS_H
#define _QHG_MESON_UDSC_DEFS_H 1

#define NGAMMAS 10								/* 1,gx,...,gt,g5,gxg5,...,gtg5 */
#define NFLAV 6									/* up,down,s+,s-,c+,c- */
#define NCHANNELS ((NGAMMAS)*(NFLAV*NFLAV))					/* all combinations of two of up,down,s+,s-,c+,c- */
#define VGF(v, g, f0, f1) ((f1) + NFLAV*((f0) + NFLAV*((g) + NGAMMAS*(v))))

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

static char flav_tags[NFLAV][256] = {"up\0", "dn\0", "s+\0", "s-\0", "c+\0", "c-\0"};

#endif /* _QHG_MESON_UDSC_DEFS_H */
