CC = g++

CCFLAGS = -O3 -I.
OBJS = bucky.o

all	: bucky mbsum

clean	:
	rm *.o bucky mbsum

bucky   : ${OBJS}
	${CC} -o $@ ${CCFLAGS} bucky.o TaxonSet.o

bucky.o : bucky.h bucky.C TaxonSet.o
	${CC} -c ${CCFLAGS} bucky.C

TaxonSet.o : TaxonSet.C TaxonSet.h
	${CC} -c ${CCFLAGS} TaxonSet.C

mbsum	: mbsum.C
	${CC} -o $@ ${CCFLAGS} mbsum.C

