FC = gfortran
compile_flags = -Wall -Ofast
link_flags = -Wall

programs = init pipe
subprograms = com fft8 lin pres prog3 prt rp rrt8 serv0 step
objects = $(addsuffix .o, $(subprograms))


all: $(addsuffix .out, $(programs))


%.out: %.o $(objects) 
	$(FC) $< $(objects) $(link_flags) -o $@


%.o: %.for
	$(FC) -c $< $(compile_flags) -o $@


clean: 
	rm -rf *.o *.out
