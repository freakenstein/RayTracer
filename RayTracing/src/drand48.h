#ifndef DRAND48H
#define DRAND48H

#include <stdlib.h>

#define m 0x100000000LL
#define c 0xB16
#define arbit 0x5DEECE66DLL

static unsigned long long seed = 1;

double drand48(void)
{
	seed = (arbit * seed + c) & 0xFFFFFFFFFFFFLL;
	unsigned int x = seed >> 16;
	return 	((double)x / (double)m);
}

void srand48(unsigned int i)
{
	seed = (((long long int)i) << 16) | rand();
}


#endif
