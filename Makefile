CC=mpif90

CFLAGS=-O3 

pipe: $(wildcard *.for)
	$(LINK.c) -I. \
	$^ $(LOADLIBS) $(LDLIBS) -o $@
	@rm -f *.o

all:  pipe

clean:
	rm -rf pipe
