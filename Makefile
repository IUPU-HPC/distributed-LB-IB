########################################
#   -- Scalable LB-IB --
#     Indiana University Purdue University Indianapolis, USA

#     @file

#     @date

#     @author Yuankun Fu
#
#########################################

## Compiler, tools and options
CC      = cc
FC      = ftn
LINK    = cc
OPT     = -g -craympich-mt

CCFLAGS = $(OPT)
LDFLAGS = $(OPT)

## Libraries
LIBS = -lpthread -lm
INC  = -I.

## FILES
OBJECTS = main.o print_info.o gen_fluid_grid.o gen_fiber_shape.o init_gv.o fiber2thread.o cube2thread.o cube2thread_and_machine.o init_eqlbrmdistrfuncDF0.o init_df1.o init_df_inout.o do_thread.o compute_bendingforce.o compute_elasticforce.o compute_stretchingforce.o fiber_SpreadForce.o fluid_get_SpreadForce.o compute_eqlbrmdistrfuncDF1.o stream_distrfunc.o bounceback_rigidwalls.o compute_rho_and_u.o fluid_SpreadVelocity.o fiber_get_SpreadVelocity.o copy_inout_to_df2.o replace_old_DF.o periodicBC.o
TARGET  = mpilbmib

##Implicit rules
.SUFFIXES: .c
.c.o:
	$(CC) -c $(CCFLAGS) $(INC) $<

## Build rules
all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(LINK) -o $@ $(OBJECTS) $(LDFLAGS) $(LIBS)

clean:
	rm -f $(OBJECTS) $(TARGET)
	rm -f *~ core
