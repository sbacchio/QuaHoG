include ../../Makefile.in

ifndef V
       V = 0
endif

ifeq ($(V), 1) 
	E = @echo > /dev/null
	C = 
else
	E = @echo
	C = @
	MAKE += --no-print-directory
endif

TARGET = main

all: $(TARGET) lib

LIBDIR = ../../lib

CFLAGS += -I$(LIBDIR)
LDFLAGS += -L$(LIBDIR)
LIBS := -lqhg $(LIBS)

SOURCES = main 

$(TARGET): ${addsuffix .o, $(SOURCES)} lib
	$(E) LD $@
	$(C)$(CC) $(CFLAGS) -o $@  ${addsuffix .o, $(SOURCES)} $(LDFLAGS) $(LIBS)

lib:
	$(E) MAKE LIB
	$(C)$(MAKE) -C $(LIBDIR)

%.o: %.c
	$(E) CC $< in $(shell pwd)
	$(C)$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(E) CLEAN in $(PWD)
	$(C)$(RM) *.o

cleanall: clean
	$(E) CLEANALL in $(PWD)
	$(C)$(RM) $(TARGET) 
	$(C)$(MAKE) -C $(LIBDIR) cleanall
