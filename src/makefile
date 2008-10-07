CC = g++

CCFLAGS = -O3
OBJS = bucky.o

all	: bucky mbsum

clean	:
	rm *.o bucky mbsum

bucky   : ${OBJS}
	${CC} -o $@ ${CCFLAGS} ${OBJS}

bucky.o : bucky.h bucky.C
	${CC} -c ${CCFLAGS} bucky.C

mbsum	: mbsum.C
	${CC} -o $@ ${CCFLAGS} mbsum.C

