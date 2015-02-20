CC=mpif90

CFLAGS=-lm -lgfortran -m64 -Ofast
SUBROUTINE_FILES=com.for prg3cycl.for prt.for rp.for serv0.for fft8.for lin1.for pres.for prog3.for rrt8.for step.for

all:  pipe.out init.out

init.out: ${SUBROUTINE_FILES} init.for
	${CC} ${SUBROUTINE_FILES} init.for -o $@ ${CFLAGS}

pipe.out: ${SUBROUTINE_FILES} pipe.for
	${CC} ${SUBROUTINE_FILES} pipe.for -o $@ ${CFLAGS}

clean:
	rm -rf *.out *.o
