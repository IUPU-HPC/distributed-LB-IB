## -- Distributed-LB-IB --
 # Copyright 2018 Indiana University Purdue University Indianapolis
 #
 # Licensed under the Apache License, Version 2.0 (the "License");
 # you may not use this file except in compliance with the License.
 # You may obtain a copy of the License at
 #
 #     http://www.apache.org/licenses/LICENSE-2.0
 #
 # Unless required by applicable law or agreed to in writing, software
 # distributed under the License is distributed on an "AS IS" BASIS,
 # WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 # See the License for the specific language governing permissions and
 # limitations under the License.
 #
 # @author: Yuankun Fu (Purdue University, fu121@purdue.edu)
 #
 # @file:
 #
 # @date:
##

## Compiler, tools and options
DEBUG ?= 0

MPICXX ?= mpicxx
CXXFLAGS ?=
CXXFLAGS += -std=c++11 -fPIC

CC      ?= mpicc
CFLAGS ?=
CFLAGS += -std=c11 -fPIC

LDFLAGS ?=

ifeq ($(strip $(DEBUG)),0)
	CXXFLAGS += -O3 -DNDEBUG
else
	CXXFLAGS += -O0 -ggdb
endif

LINK    = mpicxx

## Libraries
LIBS = -lpthread -lm
INC  = -I.

## FILES
OBJECTS = main.o print_info.o gen_fluid_grid.o gen_fiber_sheet.o init_gv.o \
	timer.o fiber2thread.o cube2thread.o cube2thread_and_task.o \
	init_eqlbrmdistrfuncDF0.o init_df1.o init_df_inout.o do_thread.o \
	compute_bendingforce.o compute_elasticforce.o compute_stretchingforce.o \
	fiber_SpreadForce.o fluid_get_SpreadForce.o \
	compute_eqlbrmdistrfuncDF1.o stream_distrfunc.o \
	bounceback_rigidwalls.o compute_rho_and_u.o \
	fluid_SpreadVelocity.o fiber_get_SpreadVelocity.o \
	copy_inout_to_df2.o replace_old_DF.o periodicBC.o
TARGET  = distributed-lb-ib

##Implicit rules
.SUFFIXES: .c
.c.o:
	$(MPICXX) -c $(CXXFLAGS) $(INC) $<

## Build rules
all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(LINK) -o $@ $(OBJECTS) $(LDFLAGS) $(LIBS)

clean:
	rm -f $(OBJECTS) $(TARGET)
	rm -f *~ core
