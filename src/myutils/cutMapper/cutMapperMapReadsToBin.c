/* cutMapperMapReads - Maps reads to an indexed genome.. */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"

#include "cutMapperCore.h"


static char const rcsid[] = "$Id: newProg.c,v 1.30 2010/03/24 21:18:33 hiram Exp $";

void usage()
/* Explain usage and exit. */
{
errAbort(
  "cutMapperMapReads - Maps reads to an indexed genome.\n"
  "usage:\n"
  "   cutMapperMapReads XXX\n"
  "options:\n"
  "   -xxx=XXX\n"
  );
}

static struct optionSpec options[] = {
   {NULL, 0},
};


void cutMapperMapReads(char *inFastq, char *indexFolder)
/* cutMapperMapReads - Maps reads to an indexed genome.. */
{
  //  FILE *f;
  char cbuff[strlen(indexFolder)+256];
  char kmerString[21];
  
  int i,j,k;
  // int n;
  int TotalMappedReads=0;
  int UniquelyMappedReads=0;
  int UnMappableReads=0;
  int LessThan10Reps=0;
  khint_t khit,indexSize;  

  clock_t t;
  clock_t t0=clock();


  khash_t(hCount_t) *hReads[256];
  khash_t(hashPos_t) *hIndex;  
  
  repvec_t *repvecp;  
  replist_t replist;
  
  //Initializing hash maps.
  for(j=0;j<256;j++) 
    hReads[j]=kh_init(hCount_t);   
  
  //Hashing input reads.
  t = clock();
  hashFastqFile(inFastq, hReads);
  verbose(1,"# Hashed fastq reads in %f seconds (%f total)\n",
	  (double)(clock() - t)/CLOCKS_PER_SEC,
	  (double)(clock() - t0)/CLOCKS_PER_SEC);
  

  for(j=0;j<256;j++){ //256
    //Load repeats!!
    sprintf(cbuff,"%s/repeats%03d.bin",indexFolder,j);
    t = clock();
    repeatFileOpen(cbuff,&replist);
    verbose(1,"# Opened repeats in %f seconds (%f total)\n",
	    (double)(clock() - t)/CLOCKS_PER_SEC,
	    (double)(clock() - t0)/CLOCKS_PER_SEC);
    
    //Load Hash! 
    sprintf(cbuff,"%s/Hash%03d.bin",indexFolder,j);
    t = clock();
    hIndex = hashFileOpen(cbuff);
    verbose(1,"# Opened hash in %f seconds (%f total)\n",
	    (double)(clock() - t)/CLOCKS_PER_SEC,
	    (double)(clock() - t0)/CLOCKS_PER_SEC);
    
    //Now we can map the reads of prefix j, traversing the hReads[j] hash
    // Repeats can be resolved now, if we have a sensitivity index available, or later
    // making a second pass, or by storing the repeats in a separate hash. 
    // If we store the repeats, we could probably do it multiple times, if it fits in memory.
    
    t = clock();
    //Traverse the count hash
    indexSize=kh_size(hIndex);
    for (k = kh_begin(hReads[j]); k < kh_end(hReads[j]); ++k)
      if (kh_exist(hReads[j], k)){
	unsigned int Suffix=kh_key(hReads[j], k);
	unsigned int num=kh_val(hReads[j], k);
	//Look up the key, on the Suffix hash table hIndex
	khit = kh_get(hashPos_t, hIndex , Suffix);
	if(kh_exist(hIndex, khit)){ // I found it!
	  TotalMappedReads+=num;
	  if(kh_value(hIndex, khit).chr==255){//Non uniquely mappable
	    //print locations... if less than XXX locations...`
	    unPackKmer((((unsigned long int)j)<<32) + Suffix, 20, kmerString);
	    fprintf(stdout,"%s %d ",kmerString,num);
	    //Go to repeat table
	    repvecp=&kv_A(replist,kh_val(hIndex, khit).pos);
	    fprintf(stdout,"repeats %d", (int) kv_size(*repvecp));
	    //Traverse repeats:
	    if(kv_size(*repvecp)<=10)
	      LessThan10Reps+=num;
	    for(i=0;i<kv_size(*repvecp);i++){
	      fprintf(stdout," (%d,%d)",kv_A(*repvecp,i).chr,kv_A(*repvecp,i).pos);
	      //fprintf("  %d\n",i);
	    }
	    fprintf(stdout,"\n");
	  }
	  else{ //Uniquely mappable
	    UniquelyMappedReads+=num;
	    unPackKmer((((unsigned long int)j)<<32) + Suffix, 20, kmerString);
	    fprintf(stdout,"%s %d at %d %d\n",kmerString,num,
		    kh_value(hIndex, khit).chr,
		    kh_value(hIndex, khit).pos);
	    //print location ... 
	  }
	}
	else{
	  UnMappableReads+=num;
	  unPackKmer((((unsigned long int)j)<<32) + Suffix, 20, kmerString);
	  fprintf(stdout,"%s %d Not mapped\n",kmerString,num);
	}
      }
    verbose(1, "#\tIndex %d: Cum. total mapped %d,"
	    " Cum Total unique %d, Cum Total <10rep %d,"
	    " Cum Total Unmappable %d (%f %%)\n"
	    ,j,TotalMappedReads,UniquelyMappedReads
	    ,LessThan10Reps,UnMappableReads
	    ,UnMappableReads/(float)(UnMappableReads+TotalMappedReads)*100);
    verbose(1,"#\tIndex %d: Mapped %d kmers in %f seconds (%f total)\n",
	    j,indexSize,
	    (double)(clock() - t)/CLOCKS_PER_SEC,
	    (double)(clock() - t0)/CLOCKS_PER_SEC);

    //DESTROYING OBJECTS 
    for (k = 0; k < kv_size(replist); ++k)
      kv_destroy(kv_A(replist,k));
    kv_destroy(replist);
    kh_destroy(hashPos_t,hIndex);    
  }/* End loop through the index file */   

  /* --------------------------------------------------------------------- */

  verbose(1,"# Total time: %f seconds = %f minutes\n",
	  (double)(clock() - t0)/CLOCKS_PER_SEC,
	  (double)(clock() - t0)/CLOCKS_PER_SEC/60);
  verbose(1,
	  "# Total mapped %d\n"
	  "# Total unique %d (%f %%)\n"
	  "# Total <=10 times %d (%f %%)\n"
	  "# Total Unmappable %d (%f %%)\n"
	  ,TotalMappedReads,
	  UniquelyMappedReads,UniquelyMappedReads/(float)TotalMappedReads*100,
	  LessThan10Reps,LessThan10Reps/(float)TotalMappedReads*100,
	  UnMappableReads,UnMappableReads/(float)(UnMappableReads+TotalMappedReads)*100);
}

/* ******************************************************************************************** */
/* ******************************************************************************************** */

int main(int argc, char *argv[])
/* Process command line. */
{
  optionInit(&argc, argv, options);
  if (argc != 3)
    usage();
  cutMapperMapReads(argv[1],argv[2]);
  return 0;
}
