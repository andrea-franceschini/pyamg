CXX           = g++

CXXFLAGS      = -O2 -std=c++11 -fopenmp $(INCLUDE)
PCXXFLAGS     = -O1 -pg -std=c++11 $(INCLUDE)

#CPLUS_INCLUDE_PATH = /usr/include/x86_64-linux-gnu/c++/4.8

CXX           = g++

INCLUDE       = -I. -I.. -I../cptEMIN -I../common

DFLAGS        = -O0 -std=c++11 -g -fopenmp -Wall $(INCLUDE)

DEST	      = .

EXTHDRS	      = ../common/DebEnv.h ../cptEMIN.h ../cptEMIN/EMIN_ImpProl.h\
		../cptEMIN/wrCSRmat.h

HDRS	      = readMat.h

LDFLAGS	      = -O2 -fopenmp -L. -L/usr/local/lib -L../lapack -L./ICHOL

LPFLAGS	      = -O1 -pg -L. -L/usr/local/lib -L../lapack -L./ICHOL

LDDFLAGS      = -O0 -g -fopenmp -L. -L/usr/local/lib -L../lapack -L./ICHOL

LIBRARIES     = -llapacke -llapack -lblas -lgfortran

LIBS          = $(LIBRARIES)

LINKER	      = $(CXX)

MAKEFILE      = Makefile

OBJS	      = main.o readMat.o

PRINT	      = pr

CAT	      = cat

PROGRAM	      = a.out

SRCS	      = main.cpp readMat.cpp
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
main.o: ../common/DebEnv.h readMat.h ../cptEMIN/wrCSRmat.h ../cptEMIN.h \
	../cptEMIN/EMIN_ImpProl.h
