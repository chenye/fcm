#include "rand_gen.h"
#include "twister/dSFMT.h"
#include <stdlib.h>

#if USE_SYS_RAND == 1


void init_random(long int seed)
{
	int i;
	srand(seed);
	for (i=0; i<100; i++)
		rand();
}

float rand0to1()
{
	return (float)rand()/RAND_MAX;
}

#else

void init_random(long int seed)
{
	int i;
	srand(seed);
	for (i=0; i<100; i++)
		rand();
	dsfmt_gv_init_gen_rand(seed);
	for (i=0; i<100; i++)
		rand0to1();
}

float rand0to1()
{
	return dsfmt_gv_genrand_close_open();
}

#endif
