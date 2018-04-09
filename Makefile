CC = gcc
CFLAGS = -m64
#linker paths to BLAS and LAPACK for use with Mac OS X 10.10 or later
LIBS = -lm /usr/lib/libblas.dylib /usr/lib/liblapack.dylib

DEPS = prototypes.h
OBJ = arc.o brick.o frame.o fsi.o main.o memory.o misc.o model.o shell.o solve.o truss.o

all: ben.exe

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)


ben.exe: $(OBJ)  
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm $(OBJ) ben.exe results*.txt
