
CXX=		/usr/bin/g++
#CXX=		/opt/local/bin/gcc
#CXXFLAGS=	-Wall -O3 -D_WITH_DEBUG 
CXXFLAGS=	-Wall -g3 -ggdb -D_WITH_DEBUG -fno-inline
#CXXFLAGS=	-Wall -pg -g -D_WITH_DEBUG 

PROG=		smorgas

LIBS=		-lz

OBJS=		smorgas.o PileupParser.o

HEAD_COMM=  smorgas.h smorgas_util.h SimpleOpt.h PileupParser.h

HEAD=		$(HEAD_COMM)


#---------------------------  Main program


all: $(PROG)

smorgas: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LIBS)


#---------------------------  Individual object files


# SimpeOpt.h is from http://code.jellycan.com/simpleopt and processes command-line args

# rebuild everything if the common headers change
$(OBJS): $(HEAD_COMM)

# rebuild the main file if any header changes
smorgas.o: $(HEAD)

PileupParser.o: PileupParser.h


#---------------------------  Other targets


clean:
	rm -f gmon.out *.o $(PROG)

clean-all: clean


#---------------------------  Obsolete and/or waiting for cleanup/reuse


