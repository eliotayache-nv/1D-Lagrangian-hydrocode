FC = gfortran
FOPT = -g -Wall -O0 -fbounds-check
# FOPT = -Ofast -fno-second-underscore -falign-loops=16 -march=corei7 -ffast-math -fassociative-math -fno-ipa-cp-clone

all: lh1

# commons.o: commons.f
# 	$(FC) $(FOPT) -c commons.f

lh1: lh1.f
	$(FC) $(FOPT) lh1.f -o lh1

clean:
	\rm -f -r *.o *mod lh1 *.dSYM *.out

