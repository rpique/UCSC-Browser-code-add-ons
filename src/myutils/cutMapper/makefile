COPT = -O2

include ../../inc/common.mk

##BINDIR= .

L += -lm ${SOCKETLIB} -lz
MYLIBDIR = ../../lib/${MACHTYPE}
MYLIBS = ${MYLIBDIR}/jkweb.a

##L += -lm
##MYLIBDIR = ../../lib/${MACHTYPE}
##MYLIBS =  ${MYLIBDIR}/jkweb.a


A = cutMapperMapReads  \
    cutMapperMakeIndex \
    cutMapperMakeMapp  \
    cutMapperFastQtoJfd \
    cutMapperMap

O = cutMapperCore.o \
    cutMapperResampler.o 

O_files=$(patsubst %.c,%.o,$(wildcard *.c))

##%.o:

cutMapperMap: cutMapperMap.o $O ${MYLIBS}
	${CC} ${COPT} ${CFLAGS} -o ${DESTDIR}${BINDIR}/$@${EXE} $^ ${L}
##	${STRIP} ${DESTDIR}${BINDIR}/${A}${EXE}

cutMapperMakeIndex: cutMapperMakeIndex.o $O ${MYLIBS}
	${CC} ${COPT} ${CFLAGS} -o ${DESTDIR}${BINDIR}/$@${EXE} $^ ${L}
##	${STRIP} ${DESTDIR}${BINDIR}/${A}${EXE}

cutMapperMakeMapp: cutMapperMakeMapp.o $O ${MYLIBS}
	${CC} ${COPT} ${CFLAGS} -o ${DESTDIR}${BINDIR}/$@${EXE} $^ ${L} 
##	${STRIP} ${DESTDIR}${BINDIR}/${A}${EXE}


cutMapperMapReads: cutMapperMapReads.o $O ${MYLIBS}
	${CC} ${COPT} ${CFLAGS} -o ${DESTDIR}${BINDIR}/$@${EXE} $^ ${L} 
##	${STRIP} ${DESTDIR}${BINDIR}/${A}${EXE}

cutMapperFastQtoJfd: cutMapperFastQtoJfd.o $O ${MYLIBS}
	${CC} ${COPT} ${CFLAGS} -o ${DESTDIR}${BINDIR}/$@${EXE} $^ ${L} 


##compile:: ${O}
##	${CC} ${COPT} ${CFLAGS} -o ${A}${EXE} ${O} ${MYLIBS} $L

all: $(O_files) $(A)


clean::
	rm -f ${A}${EXE} ${O}
	rm -f *.o

