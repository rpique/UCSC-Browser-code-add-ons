include ../../inc/common.mk

#include ../myinc/common.mk

O = $(patsubst %.c,%.o,$(wildcard *.c)) 
##Use patsubst

$(MACHTYPE)/rprlib.a: $(O) $(MACHTYPE)
	ar rcus $(MACHTYPE)/rprlib.a $(O)

$(MACHTYPE):
	mkdir $(MACHTYPE)

#test:
#	cd tests && ${MAKE} test

clean:
	rm -f ${O} $(MACHTYPE)/rprlib.a
	cd tests && ${MAKE} clean

#tags:
#	etags ../inc/*.h ../lib/*.h ../lib/*.c  ../hg/inc/*.h ../hg/lib/*.h ../hg/lib/*.c ../hg/hgTracks/*.c ../hg/hgc/*.c ../hg/hgTrackUi/*.c

