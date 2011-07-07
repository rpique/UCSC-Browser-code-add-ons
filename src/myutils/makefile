# Build all directories in utils

include ../inc/common.mk

all:    
	@for D in $(USE_DIRS) x; do \
	  if test "$$D" != "x" ; then \
	    ( cd $$D && echo $$D && $(MAKE) ) ;\
	    x=$$? ; if [ $$x -ne 0 ]; then exit $$x ; fi \
	    fi ;\
	done

scripts:
	@for F in ${SCRIPT_FILES}; do \
	    chmod +x $$F; \
	    rm -f ${SCRIPTS}/$${F}; \
	    echo cp -p $$F ${SCRIPTS}; \
	    cp -p $$F ${SCRIPTS}; \
	 done
	chmod +x genePredToBed/genePredToBed;
	rm -f ${SCRIPTS}/genePredToBed;
	echo cp -p genePredToBed/genePredToBed ${SCRIPTS};
	@cp -p genePredToBed/genePredToBed ${SCRIPTS};
	cd cvsTools && ${MAKE} scripts

alpha: scripts all

test:
	cd faToTwoBit && ${MAKE} test
	cd twoBitToFa && ${MAKE} test
	cd twoBitInfo && ${MAKE} test

clean:
	touch non-empty-rm.o
	- find . -name \*.o -print | xargs rm

SCRIPT_FILES = \
        agpToLift \
        iSync \
        lodToBedScore \
	mkMafFrames.pl

# List of directories to build
USE_DIRS = $(DIRS)

DIRS = `find ./ -mindepth 1 -maxdepth 1 -type d`

BAD_DIRS= \
	ccCpLock \
	est2genomeToPsl \
	geneStarts \
	gtfToGenePred \
	jimgrep \
	jkShell \
	dnsInfo \
	wigAsciiToBinary \
	cCp
