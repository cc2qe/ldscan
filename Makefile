CFLAGS = -O3

ldscan: ldscan.c
	gcc -o ldscan ldscan.c ${CFLAGS}
