PROGRAM=MIM

#VERSION=dev
VERSION=0.33r1#   ### rewrite lagmain.f90

INSTALL=~/bin


##### Fujitsu f90 #####
#FC=f90
#OPTION=-Am -X9 -Kfast -Et#                  # optimized
#OPTION=-Am -X9 -O0#                         # no optimized, safe
#OPTION=-Am -X9 -Et -O0 -Haesux#             # for debug


##### Intel ifort (for 10.0) #####
FC=ifort
#OPTION=-assume byterecl -O3 -warn all#    # optimized
#OPTION=-assume byterecl -O0 -warn all -heap-arrays#    # no optimized
OPTION=-assume byterecl -O0 -warn all -g -traceback -heap-arrays#    # no optimized
#OPTION=-assume byterecl -C -warn all#     # for debug


OBJS =
OBJS_MODULE =

include src/Makefile
OBJS += ${OBJS_SRC}
OBJS_MODULE += ${OBJS_MODULE_SRC}



${PROGRAM} : ${OBJS_MODULE} ${OBJS}
	${FC} ${OPTION} ${OBJS_MODULE} ${OBJS} -o $@


.SUFFIXES : .o .f90 
.f90.o :
	${FC} ${OPTION} -c $< -o $@


all : ${PROGRAM}


clean:
	rm -rf ${OBJS_MODULE} *.mod ${OBJS} ${PROGRAM}

release:
	rm -rf ${OBJS_MODULE} ${OBJS_MODULE:.o=.f90~} *.mod \
               ${OBJS} ${OBJS:.o=.f90~} ${PROGRAM} *~

install:
	cp ${PROGRAM} ${PROGRAM}-${VERSION}
	mv ${PROGRAM}-${VERSION} ${INSTALL}
	rm -f ${INSTALL}/${PROGRAM}
	ln -s ${INSTALL}/${PROGRAM}-${VERSION} ${INSTALL}/${PROGRAM}

