#ifndef RAND_GEN_H

#ifndef USE_SYS_RAND
#define USE_SYS_RAND	0
#endif

void init_random(long int seed);
float rand0to1();

#endif
