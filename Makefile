#
# This is the makefile for diablo.
# To compile the code, just type make.  Such an approach makes
# recompilation of the code easy, recompiling only as necessary
# to account for recent changes to the code.
#

COMPILER = mpif90 -vec-report0 -mcmodel=medium -shared-intel -traceback -O3 -ipo -no-prec-div -heap-arrays 99999 -march=core-avx-i -axcore-avx-i
#COMPILER = mpif90 -vec-report0 -mcmodel=medium -shared-intel -C -warn all -fp-stack-check -g -debug all -traceback
COMPOPTS =   -I/opt/intel/fc/10.1.015/include  -i-dynamic
LINKOPTS =   -L/home/bob/FFTW2_ifort/lib -lrfftw -lfftw
DLDFLAGS =   -lm /home/bob/tecplot10/lib/tecio64.a -lstdc++

# Executable file Name (Should Match with NAME_EXE in RUN file)
NAME_EXE = unstratified

# LES Option
LES = TRUE

ifeq ($(LES),TRUE)
LES_CHAN = les_chan.o les_chan_th.o
else
LES_CHAN = no_les.o
endif

$(NAME_EXE): diablo.F90 modules.o ALLOCATION.o $(LES_CHAN) \
	duct.o fft.o mpi_duct.o  dstretch.o dmgd9v.o mg_solver2d.o mudcom.o mud2cr.o flow_statistics.o boundary.o \
	grid_def
	$(COMPILER) $(COMPOPTS) diablo.F90 -o $(NAME_EXE) \
	$(LES_CHAN) \
	duct.o fft.o mpi_duct.o dstretch.o dmgd9v.o mg_solver2d.o mudcom.o mud2cr.o modules.o ALLOCATION.o flow_statistics.o boundary.o $(LINKOPTS) ${DLDFLAGS}

ifeq ($(LES),TRUE) 
les_chan.o: les_chan.F90 fft.o  grid_def
	$(COMPILER) $(COMPOPTS) -c les_chan.F90

les_chan_th.o: les_chan_th.F90 fft.o  grid_def
	$(COMPILER) $(COMPOPTS) -c les_chan_th.F90
else
no_les.o: no_les.f
	$(COMPILER) $(COMPOPTS) -c no_les.f
endif

mpi_duct.o: mpi_duct.F90  grid_def 
	$(COMPILER) $(COMPOPTS) -c mpi_duct.F90

duct.o: duct.F90 fft.o mpi_duct.o  
	$(COMPILER) $(COMPOPTS) -c duct.F90

fft.o:  fft.F90   grid_def
	$(COMPILER) $(COMPOPTS) -c fft.F90

dmgd9v.o: dmgd9v.f    grid_def
	$(COMPILER) $(COMPOPTS) -c dmgd9v.f

dstretch.o: dstretch.F90    grid_def
	$(COMPILER) $(COMPOPTS) -c dstretch.F90

mudcom.o:  mudcom.f   
	$(COMPILER) $(COMPOPTS) -c mudcom.f
      
mud2cr.o:  mud2cr.f    
	$(COMPILER) $(COMPOPTS) -c mud2cr.f     
      
mg_solver2d.o:  mg_solver2d.F90   modules.o    grid_def 
	$(COMPILER) $(COMPOPTS) -c mg_solver2d.F90
	
ALLOCATION.o:   ALLOCATION.F90  modules.o grid_def
	$(COMPILER) $(COMPOPTS) -c ALLOCATION.F90

modules.o:      modules.F90    grid_def
	$(COMPILER) $(COMPOPTS) -c modules.F90

boundary.o:      boundary.F90  grid_def
	$(COMPILER) $(COMPOPTS) -c boundary.F90

flow_statistics.o:      flow_statistics.F90  grid_def
	$(COMPILER) $(COMPOPTS) -c flow_statistics.F90

clean:
	rm -f *.o* fort.* *~ $(NAME_EXE) output.txt out_screen.txt core *.mod *.phy 

# Compiler specific notes:
#
# Compilation with Absoft Linux Fortran 77 appears to be impossible, as it
# cannot handle the INTEGER*8 option required by FFTW.  If someone finds
# a way around this, please let me know.
# 
# Compilation with Absoft Linux Fortran 90 is possible, but the option
# -YEXT_NAMES=LCS must be used as one of the link options so the compiler
# can find the lowercase external library function names.
#
# Compilation with Lahey Fortran 95 (lf95) is possible, but there is an
# underscore incompatability with the FFTW libraries, which are compiled
# with g77.  To get around this, you need to go into fft.f and add 
# trailing underscores to the name of every fftw function where they
# appear throughout the code.

