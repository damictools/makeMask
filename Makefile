CFITSIO = $(FITSIOROOT)
CPP = g++
CC = gcc
CFLAGS = -Wall -I$(CFITSIO) 
LIBS = -L$(CFITSIO) -lcfitsio -lm
GLIBS = 
GLIBS += 
OBJECTS = makeMask.o 
HEADERS = globalConstants.h

ALL : makeMask.exe
	@echo "Listo!"

makeMask.exe : $(OBJECTS)
	$(CPP) $(OBJECTS) -o makeMask.exe $(LIBS) $(GLIBS) $(CFLAGS)

makeMask.o : makeMask.cc $(HEADERS)
	$(CPP) -c makeMask.cc -o makeMask.o $(CFLAGS)

clean:
	rm -f *~ *.o *.exe
