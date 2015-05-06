CC=mpif90

CFLAGS=-lm -lgfortran -m64 -Ofast -Wall

pipe.out: $(wildcard *.for)
	$(LINK.c) -I. \
	$^ $(LOADLIBS) $(LDLIBS) -o $@
	@rm -f *.o

all:  pipe.out

clean:
	rm -rf pipe.out *.o
