CXX           = g++

CXXFLAGS      = -O2 -fopenmp $(INCLUDE)

CXX           = g++

INCLUDE       = -I. -I.. -I../BAMG_prol -I../common

DFLAGS        = -O0 -g -fopenmp -Wall $(INCLUDE)

DEST	      = .

EXTHDRS	      = ../BAMG_prol/Apply_HouHol_Rot.h ../BAMG_prol/BAMG.h\
		../BAMG_prol/BAMG_params.h ../BAMG_prol/ProlStripe_BAMG.h\
		../BAMG_prol/Reduce_IntSet.h ../BAMG_prol/SWAP.h\
		../BAMG_prol/add_new_neighs.h ../BAMG_prol/backsolve.h\
		../BAMG_prol/check_neigh_cond.h ../BAMG_prol/get_cond.h\
		../BAMG_prol/heapsort.h ../BAMG_prol/linsol_error.h\
		../BAMG_prol/maxVol.h ../BAMG_prol/maxVol_inner.h\
		../BAMG_prol/mk_HouHolVec.h ../BAMG_prol/move_row_front.h\
		../BAMG_prol/precision.h ../common/DebEnv.h\
		../common/inl_blas1.h ../cptBAMGProl.h

HDRS	      = prol_parm.h readMIS.h readMat.h

LDFLAGS	      = -O2 -fopenmp -L. -L/usr/local/lib

LDDFLAGS      = -O0 -g -fopenmp -L. -L/usr/local/lib

LIBS          = -llapacke -llapack -lblas

LINKER	      = $(CXX)

MAKEFILE      = Makefile

OBJS	      = Apply_HouHol_Rot.o BAMG.o DebEnv.o ProlStripe_BAMG.o\
		Reduce_IntSet.o add_new_neighs.o backsolve.o check_neigh_cond.o\
		get_cond.o heapsort.o linsol_error.o main.o maxVol.o\
		maxVol_inner.o mk_HouHolVec.o move_row_front.o readMIS.o\
		readMat.o

PRINT	      = pr

CAT	      = cat

PROGRAM	      = a.out

SRCS	      = Apply_HouHol_Rot.cpp BAMG.cpp DebEnv.cpp ProlStripe_BAMG.cpp\
		Reduce_IntSet.cpp add_new_neighs.cpp backsolve.cpp\
		check_neigh_cond.cpp get_cond.cpp heapsort.cpp linsol_error.cpp\
		main.cpp maxVol.cpp maxVol_inner.cpp mk_HouHolVec.cpp\
		move_row_front.cpp readMIS.cpp readMat.cpp

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) 
		$(LINKER) -o $(PROGRAM) $(LDFLAGS) $(OBJS) $(LIBS)

profile:;       @make "FFLAGS=-pg" "CFLAGS=-pg" "LDFLAGS=-pg" $(PROGRAM)

debug:;         @make "CXXFLAGS=$(DFLAGS)" "LDFLAGS=$(LDDFLAGS)" $(PROGRAM)

clean:;		rm -f $(OBJS) $(PROGRAM).cat $(PROGRAM).ftn

depend:;	mkmf -f $(MAKEFILE) PROGRAM=$(PROGRAM) DEST=$(DEST)

index:;		ctags -wx $(HDRS) $(SRCS)

install:	$(PROGRAM)
		install -s $(PROGRAM) $(DEST)

program:        $(PROGRAM)

profile:        $(PROFILE)

tags:           $(HDRS) $(SRCS); ctags $(HDRS) $(SRCS)

update:		$(DEST)/$(PROGRAM)
###
Apply_HouHol_Rot.o: ../BAMG_prol/Apply_HouHol_Rot.h ../BAMG_prol/precision.h \
	../common/inl_blas1.h
BAMG.o: ../BAMG_prol/precision.h ../common/DebEnv.h \
	../BAMG_prol/BAMG_params.h ../BAMG_prol/linsol_error.h \
	../BAMG_prol/ProlStripe_BAMG.h
DebEnv.o: ../common/DebEnv.h
ProlStripe_BAMG.o: ../common/DebEnv.h ../BAMG_prol/precision.h \
	../BAMG_prol/BAMG_params.h ../BAMG_prol/linsol_error.h \
	../BAMG_prol/check_neigh_cond.h ../BAMG_prol/add_new_neighs.h \
	../BAMG_prol/get_cond.h ../BAMG_prol/move_row_front.h \
	../common/inl_blas1.h ../BAMG_prol/mk_HouHolVec.h \
	../BAMG_prol/Apply_HouHol_Rot.h ../BAMG_prol/maxVol.h \
	../BAMG_prol/backsolve.h ../BAMG_prol/SWAP.h \
	../BAMG_prol/maxVol_inner.h ../BAMG_prol/heapsort.h \
	../BAMG_prol/Reduce_IntSet.h
Reduce_IntSet.o: ../BAMG_prol/precision.h ../common/DebEnv.h \
	../BAMG_prol/maxVol.h ../BAMG_prol/linsol_error.h \
	../BAMG_prol/move_row_front.h ../common/inl_blas1.h \
	../BAMG_prol/mk_HouHolVec.h ../BAMG_prol/Apply_HouHol_Rot.h \
	../BAMG_prol/backsolve.h ../BAMG_prol/SWAP.h \
	../BAMG_prol/maxVol_inner.h ../BAMG_prol/heapsort.h
add_new_neighs.o: ../BAMG_prol/add_new_neighs.h ../BAMG_prol/precision.h \
	../BAMG_prol/linsol_error.h
backsolve.o: ../BAMG_prol/backsolve.h ../BAMG_prol/precision.h
check_neigh_cond.o: ../BAMG_prol/check_neigh_cond.h ../BAMG_prol/precision.h \
	../BAMG_prol/linsol_error.h ../BAMG_prol/add_new_neighs.h \
	../BAMG_prol/get_cond.h ../BAMG_prol/move_row_front.h \
	../common/inl_blas1.h ../BAMG_prol/mk_HouHolVec.h \
	../BAMG_prol/Apply_HouHol_Rot.h
get_cond.o: ../BAMG_prol/get_cond.h ../BAMG_prol/precision.h \
	../BAMG_prol/move_row_front.h ../common/inl_blas1.h \
	../BAMG_prol/mk_HouHolVec.h ../BAMG_prol/Apply_HouHol_Rot.h
heapsort.o: ../BAMG_prol/heapsort.h ../BAMG_prol/precision.h \
	../BAMG_prol/SWAP.h
linsol_error.o: ../BAMG_prol/precision.h ../BAMG_prol/linsol_error.h
main.o: readMat.h readMIS.h prol_parm.h ../cptBAMGProl.h \
	../BAMG_prol/BAMG.h ../BAMG_prol/precision.h \
	../BAMG_prol/BAMG_params.h ../common/DebEnv.h
maxVol.o: ../BAMG_prol/maxVol.h ../BAMG_prol/precision.h \
	../BAMG_prol/linsol_error.h ../BAMG_prol/move_row_front.h \
	../common/inl_blas1.h ../BAMG_prol/mk_HouHolVec.h \
	../BAMG_prol/Apply_HouHol_Rot.h ../BAMG_prol/backsolve.h \
	../BAMG_prol/SWAP.h ../BAMG_prol/maxVol_inner.h
maxVol_inner.o: ../BAMG_prol/maxVol_inner.h ../BAMG_prol/precision.h \
	../BAMG_prol/linsol_error.h ../BAMG_prol/SWAP.h
mk_HouHolVec.o: ../BAMG_prol/mk_HouHolVec.h ../BAMG_prol/precision.h \
	../common/inl_blas1.h
move_row_front.o: ../BAMG_prol/move_row_front.h ../BAMG_prol/precision.h \
	../common/inl_blas1.h
