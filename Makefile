CC=mpif90

CFLAGS=-lm -lgfortran -m64 -Ofast

pipe: $(wildcard *.for)
	$(LINK.c) -I. \
	$^ $(LOADLIBS) $(LDLIBS) -o $@
	@rm -f *.o

all:  pipe

clean:
	rm -rf pipe
