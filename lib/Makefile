include ../Makefile.in

ifndef V
       V = 0
endif

ifeq ($(V), 1) 
	E = @echo > /dev/null
	C = 
else
	E = @echo
	C = @
endif

all: libqhg.a

CFLAGS += -I./
SOURCES=\
	qhg_lattice		\
	qhg_comms		\
	qhg_alloc		\
	qhg_gauge_field		\
	qhg_spinor_field	\
	qhg_correlator		\
	qhg_xchange_gauge	\
	qhg_xchange_spinor	\
	qhg_import_gauge_field	\
	qhg_point_spinor_field	\
	qhg_plaquette		\
	qhg_mom_list		\
	qhg_mesons		\
	qhg_nucleons		\
	qhg_ape_smear		\
	qhg_su3_project		\
	qhg_gauss_smear		\
	qhg_io_utils		\
	qhg_spinor_field_io	\
	qhg_spinor_field_bc	\
	qhg_spinor_field_linalg	\
	qhg_nn_thrp		\
	qhg_prop_field_linalg	\
	qhg_ft			\
	qhg_nn_sequential_sink	\
	lines_utils

include $(wildcard *.d)

libqhg.a: ${addsuffix .o, $(SOURCES)}
	$(E) LD $@
	$(C)$(AR) cr $@  ${addsuffix .o, $(SOURCES)}
	$(C)$(RANLIB) $@

%.o: %.c
	$(E) CC $<
	$(C)$(CC) $(CFLAGS) -c $< -o $@
	$(C)$(CC) -MM $(CFLAGS) $< > $*.d

clean:
	$(E) CLEAN in $(PWD)
	$(C)$(RM) *.o *.d

cleanall: clean
	$(E) CLEANALL in $(PWD)
	$(C)$(RM) libqhg.a
