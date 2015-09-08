EXEC = tssp
CC = icpc
SRC_FILES = createPsi.cpp createPotential.cpp diagonalize.cpp main.cpp
MKLROOT = /opt/intel/parallel/mkl
MKLINCLUDE = $(MKLROOT)/include
LOCALINCLUDE = /usr/local/include
LOCALLIB = /use/local/lib

CFLAGS = -c -ggdb -Wall -DMKL_ILP64 -openmp -fast -O3 -xhost -ip -qopt-report=5 -fbuiltin -ipo -no-ftz -shared-intel -qopt-report-phase=par,vec,openmp -std=c++11 -mkl=parallel -I$(MKLINCLUDE) -I$(LOCALINCLUDE)
LFLAGS = -L$(MKLROOT)/lib/intel64 -L$(LOCALLIB) -parallel -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lcfitsio -lm

O_FILES = $(SRC_FILES:.cpp=.o)

print-%  : ; @echo $* = $($*)

$(EXEC): $(O_FILES)
	$(CC) -o $@ $(O_FILES) $(LFLAGS)
	
%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@

