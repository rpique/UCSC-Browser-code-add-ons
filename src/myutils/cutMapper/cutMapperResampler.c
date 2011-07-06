#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>


#include "common.h"
#include "cutMapperCore.h"


void initrand() 
{ 
  srand((unsigned)(time(0))); 
} 

float randfloat() 
{
  double max = RAND_MAX; 
  float numb = rand()/max;
  return(numb); 
} 

double **read_sens(char **chrNames,int numChr, char *folder ) {
  FILE *my_stream;
  char buff[strlen(folder)+200];
  int i;
  int j;
  
  double **vect=malloc(numChr*sizeof(double *));
  double aux[1000000];
 
  verbose(1,"#Opening sensitivity files \n");
  
  for(j=1; j < numChr; j++) {
    sprintf(buff,"%s/%s.sensitivites.txt",folder,chrNames[j]);
    verbose(2," Loading %s\n",buff); 
    my_stream = mustOpen(buff, "r");
    for (i=0; i < 494500; i++) {
      if(feof(my_stream)) break;
      fscanf(my_stream, "%lf\n", &aux[i]);
    }
    vect[j] = malloc((i) * sizeof(double));
    memcpy(vect[j], aux, (i) * sizeof(double));
    verbose(2,"#%d;1)%lf,%d)%lf\n",j,vect[j][0],i-1,vect[j][i-1]);
    fclose(my_stream);
  }
  initrand();
  return(vect);
}

chrpos_t resample_reads(double **sens, chrpos_t *locs, int len) {
  int i;
  double tot_dens = 0;
  int index;
  //  double test;
  double prob_low=0;
  double prob_high=0;
  float random1 = randfloat();
  chrpos_t retloc;
  
  retloc.chr=0;
  retloc.pos=0xFFFFFFFF;


  for(i=0; i < len; i++) {
    index = abs(locs[i].pos/500);
    tot_dens = tot_dens + sens[locs[i].chr][index];
  }
  for(i=0; i < len; i++) {
    index = abs(locs[i].pos/500);
    prob_high = prob_low + (sens[locs[i].chr][index])/tot_dens;
    if(random1 <= prob_high && random1 >= prob_low) {
      retloc = locs[i];
      break;
    }
    prob_low = prob_high;		
  }
  
  if(retloc.pos==0xFFFFFFFF){
    verbose(1,"#RESAMPLING PROBLEM %i,%f\n",i,tot_dens);
    for(i=0; i < len; i++){
      index = abs(locs[i].pos/500);
      verbose(1,"    %d,%d,%f\n",locs[i].chr,index,sens[locs[i].chr][index]);
    }
    i = rand()%len;
    retloc = locs[i];
  }  

  return(retloc);
}


/*
int main()
{
  initrand();
  int i;
  int len=3;
  chrpos_t locs[3];
  locs[1].chr = 1;
  locs[2].chr = 3;
  locs[0].chr = 23;
  locs[1].pos = 10000;
  locs[2].pos = -32222;
  locs[0].pos = 2312344;
  double **vect = read_sens();
  for(i=1; i < 1000; i++) {
    chrpos_t chosen = resample_reads(vect, locs, len);
    printf("%i\t%i\n", chosen.chr, chosen.pos);
  }
  return(0);
}

*/
