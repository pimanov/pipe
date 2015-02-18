CC=mpif90

CFLAGS=-O3 

pipe.out: $(wildcard *.for)
	$(LINK.c) -I. \
	$^ $(LOADLIBS) $(LDLIBS) -o $@
	@rm -f *.o

all:  pipe.out

clean:
	rm -rf pipe.out
