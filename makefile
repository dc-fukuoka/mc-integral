all:
	ifort -c $(MKLROOT)/include/mkl_vsl.f90
	ifort -stand f08 -warn all -g -O3 -march=core-avx2 -fopenmp mc_gauss.F90 -mkl

clean:
	rm -f a.out *.o *.mod *~
