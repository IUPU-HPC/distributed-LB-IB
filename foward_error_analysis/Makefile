## Compiler, tools and options
DEBUG ?= 0

CXX ?= gcc
CXXFLAGS ?=
CXXFLAGS += -std=c++11 -fPIC

CC      ?= mpicc
CFLAGS ?=
CFLAGS += -std=c11 -fPIC

LDFLAGS ?=

ifeq ($(strip $(DEBUG)),0)
	CXXFLAGS += -O3
else
	CXXFLAGS += -O0 -ggdb
endif

LINK    = mpicxx

## Libraries
LIBS = -lpthread -lm
INC  = -I.

## FILES
OBJECTS = fwa.o
TARGET  = fwa

##Implicit rules
.SUFFIXES: .c
.c.o:
	$(CXX) -c $(CXXFLAGS) $(INC) $<

## Build rules
all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(LINK) -o $@ $(OBJECTS) $(LDFLAGS) $(LIBS)

clean:
	rm -f $(OBJECTS) $(TARGET)
	rm -f *~ core