#ifndef CUT_MAPPER_RESAMPLER_H
#define CUT_MAPPER_RESAMPLER_H


#include "cutMapperCore.h"

/*
typedef struct {
	unsigned char chr;
	int pos;
} chrpos_t;
*/



void initrand();

float randfloat();

double **read_sens(char **chrNames,int numChr, char *folder );

chrpos_t resample_reads(double **sens, chrpos_t *locs, int len);

#endif //CUT_MAPPER_RESAMPLER_H
