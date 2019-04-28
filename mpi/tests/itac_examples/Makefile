## Compiler, tools and options
CC      =  cc
FC      =  ftn 
LINK    =  ftn 
#INSTRUM = -finstrument-functions

##TAU, not used!
#CC       = tau_cc.sh
#FC       = tau_f90.sh 
#LINK     = tau_f90.sh

OPT     =  -O3
#ADDCFLAGS = -Wall

CCFLAGS = $(OPT) 

## Grig LAPACK
#BLASLIB    =  /home/song/GotoBLAS/libgoto.a
#LAPACKLIB  = /home/song/lapack-3.1.1/lapack_64.a -lpthread

LIBS    = -lpthread -lm $(LAPACKLIB) $(BLASLIB)

## Files
OBJECTS = qr.o semaphore.o tblas.o bins.o lv.o mailbox.o do_generate_thread.o tasknode.o do_compute_thread.o do_communicate_thread.o distribute_data.o blocknode.o resource.o debug.o core_dgeqrt.o core_dormqr.o core_dtsqrt.o core_dtsssmqr.o core_dttqrt.o core_dttssmqr.o core_dttssrfb.o core_dtsssrfb.o

TARGET  = qr


## Implicit rules
.SUFFIXES: .c
.c.o:
	$(CC) -c $(CCFLAGS) $(ADDCFLAGS) $(INC) $<
.f.o:
	$(FC) -c $(CCFLAGS) $(INC) $<
	
## Build rules
all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(LINK) -o $@ $(OBJECTS) $(LDFLAGS) $(LIBS)

clean:
	rm -f $(OBJECTS) $(TARGET)
	rm -f *~ core
