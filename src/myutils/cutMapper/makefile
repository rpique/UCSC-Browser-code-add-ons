COPT = -O2

include ../../inc/common.mk

##BINDIR= .

L += -lm ${SOCKETLIB} -lz
MYLIBDIR = ../../lib/${MACHTYPE}
MYLIBS = ${MYLIBDIR}/jkweb.a 

##L += -lm
##MYLIBDIR = ../../lib/${MACHTYPE}
##MYLIBS =  ${MYLIBDIR}/jkweb.a

##CC:= /wsu/arch/x86_64/mpich2-1.2/bin/mpicc

##CC:= /wsu/arch/x86_64/mvapich-1.2rc1/bin/mpicc

COPT= -fopenmp -O4

##CC:= mpicc

A = cutMapperMapReads  \
    cutMapperMakeIndex \
    cutMapperMakeMapp  \
    cutMapperFastQtoJfd \
    cutMapperMap \
    cutMapperMapSam \
    cutMapperOmpMapSam

O = cutMapperCore.o \
    cutMapperResampler.o 

O_files=$(patsubst %.c,%.o,$(wildcard *.c))

##%.o:

cutMapperOmpMapSam:  cutMapperOmpMapSam.o $O ${MYLIBS}
	${CC} ${COPT} ${CFLAGS} -o ${DESTDIR}${BINDIR}/$@${EXE} $^ ${L}


##cutMapperMpiMapReadsSam:  cutMapperMpiMapReadsSam.o $O ${MYLIBS}
##	${CC} ${COPT} ${CFLAGS} -o ${DESTDIR}${BINDIR}/$@${EXE} $^ ${L}


cutMapperMap: cutMapperMap.o $O ${MYLIBS}
	${CC} ${COPT} ${CFLAGS} -o ${DESTDIR}${BINDIR}/$@${EXE} $^ ${L}
##	${STRIP} ${DESTDIR}${BINDIR}/${A}${EXE}

cutMapperMapSam: cutMapperMapSam.o $O ${MYLIBS}
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

