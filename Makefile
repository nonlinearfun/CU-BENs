CC = gcc
CFLAGS = -m64
LIBS = -lm /usr/lib/libblas.dylib /usr/lib/liblapack.dylib

DEPS = prototypes.h
OBJ = arc.o brick.o frame.o fsi.o main.o memory.o misc.o model.o shell.o solve.o truss.o

all: ben.exe

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)


ben.exe: $(OBJ)  
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm $(OBJ) ben.exe