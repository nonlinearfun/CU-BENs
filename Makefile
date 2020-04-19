CC = gcc

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S), Linux)
	CFLAGS = -m64 -I/usr/include/suitesparse
	LIBS = -lm /usr/lib/x86_64-linux-gnu/libblas.so /usr/lib/x86_64-linux-gnu/liblapack.so /usr/lib/x86_64-linux-gnu/libumfpack.so
endif 

ifeq ($(UNAME_S), Darwin)
	CFLAGS = -m64
	LIBS = -lm /usr/lib/libblas.dylib /usr/lib/liblapack.dylib -lumfpack
endif

DEPS = prototypes.h
OBJ = arc.o brick.o frame.o fsi.o main.o memory.o misc.o model.o shell.o solve.o truss.o

all: ben.exe

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)


ben.exe: $(OBJ)  
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm $(OBJ) ben.exe results*.txt
