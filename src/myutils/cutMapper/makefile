COPT = -O2

include ../../inc/common.mk

##BINDIR= .

L += -lm 
MYLIBDIR = ../../lib/${MACHTYPE}
MYLIBS =  ${MYLIBDIR}/jkweb.a


A = cutMapperMapReads  \
    cutMapperMakeIndex \
    cutMapperMakeMapp  \
    cutMapperFastQtoJfd \
    cutMapperMap

O = cutMapperCore.o \
    cutMapperResampler.o 

O_files=$(patsubst %.c,%.o,$(wildcard *.c))

%.o:

cutMapperMap: cutMapperMap.o $O ${MYLIBS}
	${CC} ${COPT} -o ${DESTDIR}${BINDIR}/$@${EXE} $^
##	${STRIP} ${DESTDIR}${BINDIR}/${A}${EXE}

cutMapperMakeIndex: cutMapperMakeIndex.o $O ${MYLIBS}
	${CC} ${COPT} -o ${DESTDIR}${BINDIR}/$@${EXE} $^
##	${STRIP} ${DESTDIR}${BINDIR}/${A}${EXE}

cutMapperMakeMapp: cutMapperMakeMapp.o $O ${MYLIBS}
	${CC} ${COPT} -o ${DESTDIR}${BINDIR}/$@${EXE} $^ 
##	${STRIP} ${DESTDIR}${BINDIR}/${A}${EXE}


cutMapperMapReads: cutMapperMapReads.o $O ${MYLIBS}
	${CC} ${COPT} -o ${DESTDIR}${BINDIR}/$@${EXE} $^ 
##	${STRIP} ${DESTDIR}${BINDIR}/${A}${EXE}

cutMapperFastQtoJfd: cutMapperFastQtoJfd.o $O ${MYLIBS}
	${CC} ${COPT} -o ${DESTDIR}${BINDIR}/$@${EXE} $^ 


##compile:: ${O}
##	${CC} ${COPT} -o ${A}${EXE} ${O} ${MYLIBS} $L

all: $(O_files) $(A)


clean::
	rm -f ${A}${EXE} ${O}
	rm -f *.o

