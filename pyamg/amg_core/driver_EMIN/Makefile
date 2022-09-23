CXX           = g++

CXXFLAGS      = -O2 -std=c++11 -fopenmp $(INCLUDE)
PCXXFLAGS     = -O1 -pg -std=c++11 $(INCLUDE)

#CPLUS_INCLUDE_PATH = /usr/include/x86_64-linux-gnu/c++/4.8

CXX           = g++

INCLUDE       = -I. -I.. -I../cptEMIN -I../common

DFLAGS        = -O0 -std=c++11 -g -fopenmp -Wall $(INCLUDE)

DEST	      = .

EXTHDRS	      = ../common/DebEnv.h ../common/inl_blas1.h ../cptEMIN.h\
		../cptEMIN/DEFL_PCG_matfree.h ../cptEMIN/EMIN_ImpProl.h\
		../cptEMIN/EMIN_matfree.h ../cptEMIN/EMIN_parm.h\
		../cptEMIN/KP_spmat.h ../cptEMIN/LinvP_spmat.h\
		../cptEMIN/Orth_Q.h ../cptEMIN/Prol_add_Cnodes.h\
		../cptEMIN/Transp_Patt.h ../cptEMIN/UinvDP_spmat.h\
		../cptEMIN/apply_perm.h ../cptEMIN/copy_Prol.h\
		../cptEMIN/count_rowterms.h ../cptEMIN/cpt_Trace_Acc.h\
		../cptEMIN/ddot_par.h ../cptEMIN/dnrm2_par.h\
		../cptEMIN/gather_B_QR.h ../cptEMIN/gather_B_dump.h\
		../cptEMIN/gather_f.h ../cptEMIN/load_Jacobi.h\
		../cptEMIN/mkiat_Tglo.h ../cptEMIN/mkiat_Tloc.h\
		../cptEMIN/mvjcol.h ../cptEMIN/print_Q.h ../cptEMIN/wrCSRmat.h

HDRS	      = readMat.h

LDFLAGS	      = -O2 -fopenmp -L. -L/usr/local/lib -L../lapack -L./ICHOL

LPFLAGS	      = -O1 -pg -L. -L/usr/local/lib -L../lapack -L./ICHOL

LDDFLAGS      = -O0 -g -fopenmp -L. -L/usr/local/lib -L../lapack -L./ICHOL

LIBRARIES     = -llapacke -llapack -lblas -lgfortran
#LIBRARIES     = -llapacke -llapack -lblas -lgfortran

LIBS          = $(LIBRARIES)

LINKER	      = $(CXX)

MAKEFILE      = Makefile

OBJS	      = DEFL_PCG_matfree.o EMIN_ImpProl.o EMIN_matfree.o KP_spmat.o\
		LinvP_spmat.o Orth_Q.o Prol_add_Cnodes.o Transp_Patt.o\
		UinvDP_spmat.o apply_perm.o copy_Prol.o count_rowterms.o\
		cpt_Trace_Acc.o ddot_par.o dnrm2_par.o gather_B_QR.o\
		gather_B_dump.o gather_f.o load_Jacobi.o main.o mkiat_Tglo.o\
		mkiat_Tloc.o mvjcol.o print_Q.o readMat.o wrCSRmat.o

PRINT	      = pr

CAT	      = cat

PROGRAM	      = a.out

SRCS	      = DEFL_PCG_matfree.cpp EMIN_ImpProl.cpp EMIN_matfree.cpp\
		KP_spmat.cpp LinvP_spmat.cpp Orth_Q.cpp Prol_add_Cnodes.cpp\
		Transp_Patt.cpp UinvDP_spmat.cpp apply_perm.cpp copy_Prol.cpp\
		count_rowterms.cpp cpt_Trace_Acc.cpp ddot_par.cpp dnrm2_par.cpp\
		gather_B_QR.cpp gather_B_dump.cpp gather_f.cpp load_Jacobi.cpp\
		main.cpp mkiat_Tglo.cpp mkiat_Tloc.cpp mvjcol.cpp print_Q.cpp\
		readMat.cpp wrCSRmat.cpp
all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) 
		$(LINKER) -o $(PROGRAM) $(LDFLAGS) $(OBJS) $(LIBS)

prof:;          @make "CXXFLAGS=$(PCXXFLAGS)" "LDFLAGS=$(LPFLAGS)" $(PROGRAM)

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
DEFL_PCG_matfree.o: ../cptEMIN/EMIN_parm.h ../cptEMIN/apply_perm.h \
	../cptEMIN/Orth_Q.h ../cptEMIN/ddot_par.h ../cptEMIN/dnrm2_par.h \
	../cptEMIN/KP_spmat.h ../cptEMIN/LinvP_spmat.h \
	../cptEMIN/UinvDP_spmat.h
EMIN_ImpProl.o: ../cptEMIN/EMIN_parm.h ../cptEMIN/EMIN_matfree.h
EMIN_matfree.o: ../cptEMIN/EMIN_parm.h ../cptEMIN/wrCSRmat.h \
	../cptEMIN/Transp_Patt.h ../cptEMIN/load_Jacobi.h \
	../cptEMIN/copy_Prol.h ../cptEMIN/gather_B_QR.h \
	../cptEMIN/gather_B_dump.h ../cptEMIN/gather_f.h \
	../cptEMIN/cpt_Trace_Acc.h ../cptEMIN/print_Q.h \
	../cptEMIN/Prol_add_Cnodes.h ../cptEMIN/DEFL_PCG_matfree.h
LinvP_spmat.o: ../cptEMIN/LinvP_spmat.h
Transp_Patt.o: ../cptEMIN/count_rowterms.h ../cptEMIN/mkiat_Tloc.h \
	../cptEMIN/mkiat_Tglo.h ../cptEMIN/mvjcol.h
UinvDP_spmat.o: ../cptEMIN/UinvDP_spmat.h
count_rowterms.o: ../cptEMIN/count_rowterms.h
gather_B_QR.o: ../common/inl_blas1.h ../common/DebEnv.h \
	../cptEMIN/EMIN_parm.h
gather_B_dump.o: ../common/DebEnv.h ../cptEMIN/EMIN_parm.h
main.o: ../common/DebEnv.h readMat.h ../cptEMIN/wrCSRmat.h ../cptEMIN.h \
	../cptEMIN/EMIN_ImpProl.h
mkiat_Tglo.o: ../cptEMIN/mkiat_Tglo.h
mkiat_Tloc.o: ../cptEMIN/mkiat_Tloc.h
