CC          = icpc
CLINKER     = icpc

#CFLAGS      =   -Wall -O3 -march=pentium
CFLAGS      =   -Wall -O3 -xHost
#CFLAGS      = -i-fast -lm  
#LIBS        = -lm
DEPEND= makedepend

SRC        = voronoi.c system.c readinput.c write.c bubblesort.c work.c
OBJS       = voronoi.o system.o readinput.o write.o bubblesort.o work.o
EXECS      = voronoi

default: voronoi

all: $(EXECS)

voronoi:$(OBJS)
	$(CLINKER) $(OPTFLAGS) -o voronoi $(OBJS) $(LIBS)

clean:
	/bin/rm -f *.o *~ $(EXECS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c

voronoi.o: system.h
system.o: system.h
readinput.o: system.h
write.o: system.h
bubblesort.o: system.h
work.o: system.h
