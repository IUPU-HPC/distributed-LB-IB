########################################
#
# IUPUI ScalableIB (integrated withLBM)
#
#########################################

## Compiler, tools and options
CC      = gcc
LINK    = gcc
OPT     = -O3

CCFLAGS = $(OPT) -Wall 
LDFLAGS = $(OPT)

## Libraries
LIBS = -lpthread -lm  
INC  = -I.

## FILES
OBJECTS1 = lbm_ib.o
TARGET1  = scalable_IB_v4.1


##Implicit rules
.SUFFIXES: .c
.c.o:
	$(CC) -c $(CCFLAGS) $(INC) $<

## Build rules
all: $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4)

$(TARGET1): $(OBJECTS1)
	$(LINK) -o $@ $(OBJECTS1) $(LDFLAGS) $(LIBS)

clean:
	rm -f $(OBJECTS1) $(TARGET1)  
	rm -f *~ core
