	objects = sph.o input.o output.o time_integration.o\
		  virt_part.o single_step.o direct_find.o\
		  kernel.o internal_force.o eos.o density.o\
		  hsml.o external_force.o art_visc.o viscosity.o\
		  art_heat.o av_vel.o time_print.o time_elapsed.o\
		  link_list.o grid_geom.o init_grid.o

# compiler
FC = f77

# compiler flags
FCFLAGS = -g -c

SRCS = $(patsubst %.f, %.o, $(wildcard *.f))

PROGRAM = edit

all: $(PROGRAM)

$(PROGRAM): $(SRCS)
	$(FC) $(FCFLAGS) -o $@ $^

%.o: %.f
	$(FC) $(FCFLAGS) -o $@ $<


edit:$(objects)
	 f77 -g -o edit $(objects)

sph.o:param.inc
input.o:param.inc
output.o:param.inc
time_integration.o:param.inc
virt_part.o:param.inc
single_step.o:param.inc
direct_find.o:param.inc
link_list.o:param.inc
grid_geom.o:param.inc
init_grid.o:param.inc
kernal.o:param.inc
internal_force.o:param.inc
eos.o:param.inc
density.o:param.inc
hsml.o:param.inc
external_force.o:param.inc
art_visc.o:param.inc
external_force.o:param.inc
viscosity.o:param.inc
art_heat.o:param.inc
av_vel.o:param.inc
time_print.o:time_print.f90
	f77 -c time_print.f90
time_elapsed.o:time_elapsed.f90
	f77 -c time_elapsed.f90

.PHONY:clean
clean:
	rm edit $(objects)

