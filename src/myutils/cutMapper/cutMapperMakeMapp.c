/* cutMapperMakeMapp - Uses the hash 20bp indexed genome to make a mappability map. */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"

#include "cutMapperCore.h"

#include "kvec.h"

#include "xbfiles.h"


static char const rcsid[] = "$Id: newProg.c,v 1.30 2010/03/24 21:18:33 hiram Exp $";

void usage()
/* Explain usage and exit. */
{
errAbort(
  "cutMapperMakeMapp - Uses the hash 20bp indexed genome to make a mappability map\n"
  "usage:\n"
  "   cutMapperMakeMapp indexFolder mappability.b8\n"
  "options:\n"
  "   -xxx=XXX\n"
  );
}

/* */

static struct optionSpec options[] = {
   {NULL, 0},
};


void cutMapperMakeMapp(char *indexFolder, char *outBinFile)
/* cutMapperMakeMapp - Uses the hash 20bp indexed genome to make a mappability map. */
{
  
  char cbuff[strlen(indexFolder)+256];

  int i,j;
  khiter_t k;
  int n,num,pos,sum;

  unsigned long long hist[256];

  khash_t(hashChr_t) *hChr; 
  char **chromNames;
  unsigned *chromSizes;
  unsigned int rloc;

  int AllZero;
  int AllEqual;

  xbList_t *xbl;

  clock_t t;
  clock_t t0=clock();
     
  khash_t(hashPos_t) *h;  
  repvec_t *repvecp;  
  replist_t replist;

  
  //Open chromosome structure.
  verbose(1,"# Loading Chromosome sizes\n");
  sprintf(cbuff,"%s/chromSizes.txt",indexFolder);
  hChr = loadChromSizes(cbuff, &chromNames, &chromSizes);
  for(j=0;j<kh_end(hChr);j++)
    if(kh_exist(hChr,j)){
      rloc=kh_val(hChr,j);
      verbose(2,"## %d %s %s %d\n",rloc,kh_key(hChr,j),chromNames[rloc],chromSizes[rloc]);
    }
  
  //Make file structure. 
  verbose(1,"# Initializing mappability map for %d chromosomes\n",kh_size(hChr));
  xbl=xbInit(0,kh_size(hChr),chromNames,chromSizes);
  
  verbose(2,"## xbl Stranded %d\n",xbl->isStranded); 

  //mmap each chromosome?
  
  //cycle throught the 256 hashes, 
  // increase +1 all the positons for each kmer.
 
  verbose(1,"# Calculating unique mappability map\n");  
  //for(j=57;j<58;j++){ 
  for(j=0;j<256;j++){ 

    //Load Hash! 
    sprintf(cbuff,"%s/Hash%03d.bin",indexFolder,j);
    t = clock();
    h = hashFileOpen(cbuff);
    verbose(1,"# Opened hash in %f seconds (%f total)\n",
	    (double)(clock() - t)/CLOCKS_PER_SEC,
	    (double)(clock() - t0)/CLOCKS_PER_SEC);

    //Traverse hash and repeats
    t = clock();
    for (k = kh_begin(h); k < kh_end(h); ++k)
      if (kh_exist(h, k)){
	if (kh_value(h, k).chr==255){	
	  continue; //We will handle repeats later
	  //verbose(2,"##  rpeat\n");
	}else{//Kmer-maps uniquely 
	  pos=abs(kh_val(h, k).pos)-1;
	  //verbose(2,"##  unique %d %d\n",(int)kh_val(h, k).chr,pos);
	  n=xbl->vec[kh_val(h, k).chr].a[pos]; 
	  if(n>1)
	    xbl->vec[kh_val(h, k).chr].a[pos]=135;
	  else
	    xbl->vec[kh_val(h, k).chr].a[pos]=1;

	  //fprintf(f,"%x\t%s\t%s\t%d\n", kh_key(h, k), buff ,chrNames[kh_val(h, k).chr], kh_val(h, k).pos );	
	}
      }
    verbose(1,"# Calculated mappability for index %d in %f seconds (%f total)\n",j,
	    (double)(clock() - t)/CLOCKS_PER_SEC,
	    (double)(clock() - t0)/CLOCKS_PER_SEC);

    kh_destroy(hashPos_t,h);    
  }
  
  verbose(1,"Calculating extended mappability map\n");  
  //for(j=57;j<58;j++){ 
  for(j=0;j<256;j++){ 
    //Load repeats!!
    sprintf(cbuff,"%s/repeats%03d.bin",indexFolder,j);
    t = clock();
    repeatFileOpen(cbuff,&replist);
    verbose(1,"#1) Opened repeats index %d in %f seconds (%f total)\n",j,
	    (double)(clock() - t)/CLOCKS_PER_SEC,
	    (double)(clock() - t0)/CLOCKS_PER_SEC);
    
    //Loop through repeat table!
    for(k = 0; k < kv_size(replist);k++){
      repvecp=&kv_A(replist,k);
      n=kv_size(*repvecp);
      //Traverse repeats: if any has not 0 
      num=0;
      AllZero=1;
      AllEqual=1;
      for(i=0;i<n;i++){
	//Check if 255 or more locations would map at this position.
	pos=abs(kv_A(*repvecp,i).pos)-1;
	AllZero=(AllZero && (xbl->vec[kv_A(*repvecp,i).chr].a[pos]==0));
	AllEqual=(AllEqual && (xbl->vec[kv_A(*repvecp,i).chr].a[pos]==n));
      }
      if(AllZero)
	for(i=0;i<n;i++){
	  pos=abs(kv_A(*repvecp,i).pos)-1;
	  xbl->vec[kv_A(*repvecp,i).chr].a[pos]=n;
	}
      else if(!AllEqual){
	for(i=0;i<n;i++){
	  pos=abs(kv_A(*repvecp,i).pos)-1;
	  xbl->vec[kv_A(*repvecp,i).chr].a[pos]=130; //Flag SNP overlapping
	}	  
      }
    } 
    //Destroy repeats
    for (k = 0; k < kv_size(replist); ++k)
      kv_destroy(kv_A(replist,k));
    kv_destroy(replist);
  }

  verbose(1,"Calculating extended mappability map 2\n");  
  //for(j=57;j<58;j++){ 
  for(j=0;j<256;j++){ 
    //Load repeats!!
    sprintf(cbuff,"%s/repeats%03d.bin",indexFolder,j);
    t = clock();
    repeatFileOpen(cbuff,&replist);
    verbose(1,"#2) Opened repeats index %d in %f seconds (%f total)\n",j,
	    (double)(clock() - t)/CLOCKS_PER_SEC,
	    (double)(clock() - t0)/CLOCKS_PER_SEC);
    
    //Loop through repeat table!
    for(k = 0; k < kv_size(replist);k++){
      repvecp=&kv_A(replist,k);
      n=kv_size(*repvecp);
      //Traverse repeats: if any has not 0 
      num=0;
      AllZero=1;
      AllEqual=1;
      for(i=0;i<n;i++){
	//Check if 255 or more locations would map at this position.
	pos=abs(kv_A(*repvecp,i).pos);
	AllZero=(AllZero && (xbl->vec[kv_A(*repvecp,i).chr].a[pos]==0));
	AllEqual=(AllEqual && (xbl->vec[kv_A(*repvecp,i).chr].a[pos]==n));
      }
      if(AllZero)
	for(i=0;i<n;i++){
	  //Check if 255 or more locations would map at this position.
	  pos=abs(kv_A(*repvecp,i).pos)-1;
	  xbl->vec[kv_A(*repvecp,i).chr].a[pos]=n;
	}
      else if(!AllEqual){
	for(i=0;i<n;i++){
	  //Check if 255 or more locations would map at this position.
	  pos=abs(kv_A(*repvecp,i).pos)-1;
	  xbl->vec[kv_A(*repvecp,i).chr].a[pos]=131; //Flag SNP overlapping
	}	  
      }
    } 
    //Destroy repeats
    for (k = 0; k < kv_size(replist); ++k)
      kv_destroy(kv_A(replist,k));
    kv_destroy(replist);
  }

  
  //
  t = clock();
  verbose(1,"Calculating mappability histogram \n");
  for(j=0;j<256;j++)
    hist[j]=0;
  for(j=1;j<kh_size(hChr);j++){
    verbose(2,"## for %d) %s \n",j,chromNames[j]);
    for(k=0;k<chromSizes[j];k++)
      hist[xbl->vec[j].a[k]]++;
  }
  for(j=0;j<140;j++)
    fprintf(stdout,"%d\t",j);
  fprintf(stdout,"\n");
  for(j=0;j<140;j++)
    fprintf(stdout,"%llu\t",hist[j]);
  fprintf(stdout,"\n");
  verbose(1,"# Time for this step = %f seconds (%f total)\n",
	    (double)(clock() - t)/CLOCKS_PER_SEC,
	    (double)(clock() - t0)/CLOCKS_PER_SEC);


  verbose(2,"## xbl Stranded %d\n",xbl->isStranded);

  //Write file..
  t = clock();
  xbDump(xbl,outBinFile);
  verbose(1,"# Time for this step = %f seconds (%f total)\n",
	    (double)(clock() - t)/CLOCKS_PER_SEC,
	    (double)(clock() - t0)/CLOCKS_PER_SEC);

  verbose(1,"# Freeing memory\n");
  xbDestroy(xbl);

  t = clock();
  verbose(1,"# Reloading by mmaping the mappability file to double-check everything is OK.\n");
  xbl=xbLoadMmap(outBinFile);
  verbose(1,"# Time for this step = %f seconds (%f total)\n",
	    (double)(clock() - t)/CLOCKS_PER_SEC,
	    (double)(clock() - t0)/CLOCKS_PER_SEC);

  t = clock();
  verbose(1,"Calculating mappability histogram \n");
  for(j=1;j<kh_size(hChr);j++)
    for(k=0;k<chromSizes[j];k++){
      //if((j==24)&&((k % 1000000)==0))verbose(2,"%d\n",k);
      hist[xbl->vec[j].a[k]]--;
    }

  sum=0;
  for(j=0;j<256;j++)
    sum+=abs(hist[j]);
  //    verbose(1,"%llu\t",hist[j]);
  verbose(1,"# Errors %d\n",sum);
  verbose(1,"# Time for this step = %f seconds (%f total)\n",
	    (double)(clock() - t)/CLOCKS_PER_SEC,
	    (double)(clock() - t0)/CLOCKS_PER_SEC);
}

int main(int argc, char *argv[])
/* Process command line. */ 
{
  optionInit(&argc, argv, options);
  if (argc != 3)
    usage(); 
  cutMapperMakeMapp(argv[1],argv[2]);
  return 0;
}
