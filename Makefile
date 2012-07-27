SHELL = /bin/bash -O extglob

#--------------------------------------------------
#variables
#--------------------------------------------------
#compilation flag
OMP = yes
FC = ifort
FCFLAGS = -module $(BIN) -O3
OMPFLAG = -openmp -parallel -fpp

ifeq ($(FC),gfortran)
    FCFLAGS = -J$(BIN) -O3
    OMPFLAG = -fopenmp
endif

ifeq ($(OMP),yes)
    FCFLAGS += $(OMPFLAG)
endif

#directory for executables
BIN = bin

#--------------------------------------------------
#compiling
#--------------------------------------------------
all: checkdir 1D 2D

#mkdir
checkdir:
	mkdir -p $(BIN)

#build executables
1D : checkdir
	$(FC) $(FCFLAGS) -o $(BIN)/UGKS1D src/UGKS1D.f90

2D : checkdir
	$(FC) $(FCFLAGS) -o $(BIN)/UGKS2D src/UGKS2D.f90

#build manual
manual: 
	cd doc; latex  manual.tex
	cd doc; bibtex manual
	cd doc; latex  manual
	cd doc; latex  manual
	cd doc; dvips  manual
	cd doc; ps2pdf manual.ps

#clean
clean:
	rm -f bin/*
	cd doc; rm -f !(*.tex|*.bib|*.pdf)
