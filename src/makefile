CC = g++

CCFLAGS = -O3 -I.
OBJS = bucky.o

all	: bucky mbsum

clean	:
	rm *.o bucky mbsum

bucky   : ${OBJS}
	${CC} -o $@ ${CCFLAGS} bucky.o TaxonSet.o mbsumtree.o

bucky.o : bucky.h bucky.C TaxonSet.o mbsumtree.o
	${CC} -c ${CCFLAGS} bucky.C

TaxonSet.o : TaxonSet.C TaxonSet.h
	${CC} -c ${CCFLAGS} TaxonSet.C

mbsumtree.o : mbsumtree.C mbsumtree.h
	${CC} -c ${CCFLAGS} mbsumtree.C

mbsum	: mbsum.C mbsumtree.o
	${CC} -o $@ ${CCFLAGS} mbsum.C mbsumtree.o

