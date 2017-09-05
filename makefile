all : nlts prc lts

lts.o : lts.h lts.c
	gcc  -c lts.c 
driver.o : lts.h driver.c lts.c
	gcc  -c driver.c 
nlts.o : nlts.h nlts.c
	gcc  -c nlts.c 
ndriver.o : nlts.h ndriver.c nlts.c
	gcc  -c ndriver.c 
out.o : lts.h out.c
	gcc  -c out.c 
nout.o : nlts.h nout.c
	gcc  -c nout.c 
radau5.o: lts.h radau5.f
	gfortran    -c radau5.f
decsol.o: lts.h decsol.f
	gfortran  -c decsol.f
lts:  radau5.o decsol.o out.o driver.o lts.o lts.h cblock.h
	gfortran    -o lts radau5.o decsol.o  out.o  driver.o -lm 
prcdriver.o : lts.h prcdriver.c prclts.c
	gcc  -c prcdriver.c 
prclts.o : lts.h prclts.c
	gcc  -c prclts.c 
prcout.o : lts.h prcout.c
	gcc  -c prcout.c 
prc:  radau5.o decsol.o prcout.c prcdriver.c prclts.c lts.h cblock.h prcout.o prcdriver.o 
	gfortran   -o prc radau5.o decsol.o  prcout.o  prcdriver.o -lm 
nlts:  radau5.o decsol.o nout.o ndriver.o nlts.o nlts.h cblock.h
	gfortran    -o nlts radau5.o decsol.o  nout.o  ndriver.o -lm 

cleanall :
	rm -fR *.o nlts prc lts
