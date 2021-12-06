CC = gcc

MUMPSroot = /Users/bkorb/SparseSolver/Multifrontal/mumps-main
CFlags_MUMPS = -O3
INCLUDE_MUMPS = -I$(MUMPSroot)/cache/MUMPS-5.4.1-src/include -I$(MUMPSroot)/build/_deps/mumps-build/include -I$(MUMPSroot)/cache/MUMPS-5.4.1-src/src -I$(MUMPSroot)/cache/MUMPS-5.4.1-src/PORD/lib/../include -I/opt/homebrew/Cellar/open-mpi/4.1.1_2/include -I/opt/homebrew/Cellar/open-mpi/4.1.1_2/lib
LINKING_MUMPS = -L/opt/homebrew/Cellar/libevent/2.1.12/lib -L/opt/homebrew/Cellar/libevent/2.1.12/lib -L/opt/homebrew/Cellar/gcc/11.2.0_3/lib/gcc/11/gcc/aarch64-apple-darwin20/11  -L/opt/homebrew/Cellar/gcc/11.2.0_3/lib/gcc/11  $(MUMPSroot)/build/libsmumps.a $(MUMPSroot)/build/libdmumps.a $(MUMPSroot)/build/libmumps_common.a $(MUMPSroot)/build/lib/libscalapack.a $(MUMPSroot)/build/lib/libblacs.a /Library/Developer/CommandLineTools/SDKs/MacOSX11.3.sdk/usr/lib/liblapack.tbd /Library/Developer/CommandLineTools/SDKs/MacOSX11.3.sdk/usr/lib/libblas.tbd /opt/homebrew/Cellar/open-mpi/4.1.1_2/lib/libmpi_usempif08.dylib /opt/homebrew/Cellar/open-mpi/4.1.1_2/lib/libmpi_usempi_ignore_tkr.dylib /opt/homebrew/Cellar/open-mpi/4.1.1_2/lib/libmpi_mpifh.dylib /opt/homebrew/Cellar/open-mpi/4.1.1_2/lib/libmpi.dylib $(MUMPSroot)/build/libpord.a -lgfortran -lgcc_s.1.1 -lgcc -lm -lgcc_s.1.1 -lgcc -lgcc

SuperLUroot     = /Users/bkorb/SparseSolver/Multifrontal/superlu-5.3.0/build
CFlags_SuperLU = -O3 -DNDEBUG -DPRNTlevel=0 -DDEBUGlevel=0
INCLUDE_SuperLU = -I/Users/bkorb/SparseSolver/Multifrontal/superlu-5.3.0/SRC
LINKING_SuperLU = -lm $(SuperLUroot)/SRC/libsuperlu.a $(SuperLUroot)/CBLAS/libblas.a

UMFroot = /Users/bkorb/SparseSolver/Multifrontal/SuiteSparse-5.10.1
CFLAGS_UMF = -O3 -fexceptions -fPIC -fno-common
INCLUDE_UMF = -I$(UMFroot)/include
LINKING_UMF = -lm -rpath $(UMFroot)/lib -L$(UMFroot)/lib -lumfpack -lamd -lsuitesparseconfig -lccolamd -lcamd $(UMFroot)/lib/libmetis.dylib -llapack -lblas

all: MumpsTest SuperLUTest UMFTest

MumpsTest: MumpsTest.o
	$(CC) MumpsTest.o $(LINKING_MUMPS) -o MumpsTest

MumpsTest.o:
	$(CC) $(CFlags_MUMPS) $(INCLUDE_MUMPS) -c MumpsTest.c -o MumpsTest.o


SuperLUTest: SuperLUTest.o
	$(CC) SuperLUTest.o $(LINKING_SuperLU) -o SuperLUTest

SuperLUTest.o:
	$(CC) $(CFLAGS_SuperLU) $(INCLUDE_SuperLU) -c SuperLUTest.c -o SuperLUTest.o


UMFTest: UMFTest.o
	$(CC) UMFTest.o $(LINKING_UMF) -o UMFTest

UMFTest.o: UMFTest.c
	$(CC) $(CFLAGS_UMF) $(INCLUDE_UMF) -c UMFTest.c -o UMFTest.o

clean:
	rm UMFTest SuperLUTest MumpsTest *.o