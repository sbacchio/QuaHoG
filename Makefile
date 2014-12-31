include ./Makefile.in

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

DIRS = mainprogs lib

all: 
	-$(C)for d in $(DIRS); do ($(MAKE) -C $$d); done

clean:
	-$(C)for d in $(DIRS); do ($(MAKE) -C $$d clean); done

cleanall: clean
	-$(C)for d in $(DIRS); do ($(MAKE) -C $$d cleanall); done
