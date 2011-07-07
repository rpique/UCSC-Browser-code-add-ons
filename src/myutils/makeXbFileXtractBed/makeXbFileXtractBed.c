/* makeXbFileXtractBed - Makes a temporal binary file and extracts the regions sorrounding a list of beds. */
//#include <stdlib.h>

#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"

static char const rcsid[] = "$Id: newProg.c,v 1.30 2010/03/24 21:18:33 hiram Exp $";

#include "khash.h"
#include "kvec.h"

#include "xbfiles.h"

KHASH_MAP_INIT_STR(hashChr_t, unsigned char)


void usage()
/* Explain usage and exit. */
{
errAbort(
  "makeXbFileXtractBed - Makes a temporal binary file and extracts the regions sorrounding a list of beds\n"
  "usage:\n"
  "   makeXbFileXtractBed chromSizes.txt input.txt regions.bed output.txt\n"
  "options:\n"
  //  "   -stranded = Make a stranded xb file (default yes)\n"
  "   -window = Window around the motif to expand \n"
  "   -readSize = Size of the reads (20bp default) so the minus strand is shifted\n"
  //  "   -binOut Makes the output in binary format rather than text"
  //  "   -numBytes = Number of bytes used to count reads at each position (default 1), counts up to 255\n"
  "chromSizes.txt File with the size of the chromosomes\n"
  "   chr1\t122143214\n"
  "   chr2\t133978974\n"
  "   ...\n"
  "input.txt File with the read locations\n"
  "   chr1\t10000\t+\n"
  "   chr1\t13123\t-\n"
  "   ...\n"
  "regions.bed is a tab delimited file with at least four columns: \n"
  "   Chr1\t1000100\t1000150\t+\t[Other Columns]\n"
  "   .... \n"
  "   ChrY\t1000100\t1000150\t-\t[Other Columns]\n"
  "output.txt is in the following form: \n"
  "   Instance1: Forward reads (2*Window+MotifLen columns), Reverse strand reads (2*Window+MotifLen columns).\n"
  " .... \n\n"      
  );
}

int window = 0;      /* size of the window in which we extend the bed region. */
int readSize = 20;   

static struct optionSpec options[] = {
   {"window", OPTION_INT},
   {"readSize", OPTION_INT},
   {NULL, 0},
};

/* Repeated from cutMapperCore.c */
khash_t(hashChr_t) *loadChromSizes(char *chromFile, char ***chromNames, unsigned **chromSizes)
{
  FILE *f=mustOpen(chromFile,"r");
  int chrnum=0;
  char buff[256];
  char empty[]="---";
  //char *buff2;
  unsigned size; 
  khash_t(hashChr_t) *hChr= kh_init(hashChr_t); 
  
  khiter_t k;
  int hret;

  kvec_t(unsigned) kvChromSizes;
  kvec_t(char *)   kvChromNames;
  
  kv_init(kvChromSizes); 
  kv_init(kvChromNames);

  //Pushing the empty crhomosome.
  //buff=cloneString(empty);
  k = kh_put(hashChr_t, hChr , cloneString(empty), &hret);
  kv_push(char *,kvChromNames, (char *)kh_key(hChr, k));
  kv_push(unsigned,kvChromSizes,0);
    
  while(!feof(f)){
    chrnum++; 
    assert(chrnum<256);
    fscanf(f,"%s\t%u\t", buff, &size);
    verbose(2,"#%d) %s\t%u\n",chrnum, buff, size);
    //buff2=cloneString(buff);
    //kv_push(char *,kvChromNames, buff2 );
    k = kh_put(hashChr_t, hChr , cloneString(buff), &hret);
    assert(hret==1);
    kh_val(hChr, k) = chrnum;
    kv_push(char *,kvChromNames,(char *)kh_key(hChr, k));
    kv_push(unsigned,kvChromSizes,size);
  }
  kv_trim(char *,kvChromNames);
  kv_trim(unsigned,kvChromSizes);
  
  *chromNames=kvChromNames.a;
  *chromSizes=kvChromSizes.a;

  return hChr;
}


void extractBedFromXbFile(xbList_t *xbl, char *bedFileName, char *outFileName)
/* extractBedFromXbFile - Extract the region sorrounding a list of beds from an Xb binary file. */
{
  //FILE *bedF=mustOpen(bedFileName,"r");
  struct lineFile *lf = lineFileOpen(bedFileName, TRUE);

  FILE *outF=mustOpen(outFileName,"w");

  //  xbList_t *xbl=xbLoadMmap(xbFileName);

  int skipLine=0;
  long int cF,cR;
  //char buff[1024];
  //char chr_str[512];

  int iChr,j;
  //  int MaxCutSites;
  int left,right;
  int count=0;
  //  long int offset;
  //  double offsetCheck;
  char cStrand;

  char *row[10];
  int wordCount;
  char *chr_str;
  
  //Hash chromosome names? instead of strcmp??  
  cF=0;cR=0;
  while ((wordCount = lineFileChop(lf, row)) != 0){
  //while(!feof(bedF)){
    //fscanf(bedF,"%s\t%d\t%d\t%c\%[^\n]\n",chr_str,&left,&right,&cStrand,buff);
    assert(wordCount>=6);
    chr_str = row[0];
    left = lineFileNeedNum(lf, row, 1);
    right = lineFileNeedNum(lf, row, 2);
    cStrand = row[5][0];
  

    // I should use the hash
    for(iChr=0;iChr<xbl->count;iChr++)
      if(strcmp(xbl->names[iChr],chr_str)==0)
	break;
    if(iChr<xbl->count){      

      left=left-window;
      right=right+window;
      
      if((left >= 0) && (right < (xbl->sizes[iChr]))){
	if(cStrand=='+'){
	  for(j=left;j<=right;j++) 
	    fprintf(outF,"\t%d",xbl->vec[(iChr<<1)].a[j]);
	  for(j=left;j<=right;j++)
	    fprintf(outF,"\t%d",xbl->vec[(iChr<<1)+1].a[j]);
	  fprintf(outF,"\n");
	}
	else if(cStrand=='-'){
	  for(j=right;j>=left;j--)
	    fprintf(outF,"\t%d",xbl->vec[(iChr<<1)+1].a[j]);
	  for(j=right;j>=left;j--)
	    fprintf(outF,"\t%d",xbl->vec[(iChr<<1)].a[j]);	    
	  fprintf(outF,"\n");
	}
	else{
	  verbose(1,"# Unrecognized strand parameter!!!\n");	
	  skipLine=1;
	} 
	
	if(skipLine==0){
	  for(j=left;j<=right;j++)
            cF+=xbl->vec[iChr].a[j];
          for(j=left;j<=right;j++)
            cR+=xbl->vec[iChr].a[j];
	}
      }
      else{
	verbose(1,"# Skipping segment off-limits\n");
	skipLine=1;
      }      

      if((count%1000)==0)
	verboseDot();
      //verbose(1,"# Processed lines: %d\r",count);
    }
    else{
      verbose(1,"# Unreq %s:%d-%d in line %d !!!\n",chr_str,left,right,count);
      skipLine=1;
    }
    
    if(skipLine==1){
      fprintf(outF,"NA\n");
      skipLine=0;
    }    
    count++;
  }//end while(!feof(bedF)){

  
  //  carefulClose(&bedF);
  lineFileClose(&lf);

  carefulClose(&outF);
  verbose(1,"# Mapped %ld on F, %ld on R strands \n",cF,cR);
  verbose(1,"# Succesfully processed %d segments\n",count);

  //free xbl?

}


khash_t(hashChr_t) *makeXbFile(char *chromSizesFile, char *inputReadsFile, char *outXbFile, xbList_t **pxbl)
/* makeXbFile - Makes a binary file with chromosome sizes and data type XXX. */
{
  int j;
  unsigned int rloc;
  struct lineFile *lf = lineFileOpen(inputReadsFile, TRUE);
  //  FILE *inF=mustOpen(inputReadsFile,"r");
  xbList_t *xbl;

  khash_t(hashChr_t) *hChr; 
  char **chromNames;
  unsigned *chromSizes;

  //  char *chr_str;

  char *row[3];
  int wordCount;
  char *chr_str;


  char cStrand;
  int left;
  int iChr;
  long int count;

  khiter_t khit;

  
  /* --------------------------------------------------------------------- */
  
  verbose(1,"Loading Chromosome sizes \n");
  hChr = loadChromSizes(chromSizesFile, &chromNames, &chromSizes);
  for(j=0;j<kh_end(hChr);j++)
    if(kh_exist(hChr,j)){  
      rloc=kh_val(hChr,j);
      verbose(2,"%d %s %s %d\n",rloc,kh_key(hChr,j),chromNames[rloc],chromSizes[rloc]);
    }
  
  /* --------------------------------------------------------------------- */

  //Initializing memmory mapped file. 
  xbl=xbInitMmap(1,kh_size(hChr),chromNames,chromSizes,outXbFile);
  
  /* --------------------------------------------------------------------- */

  verbose(2,"Reading in %s ...\n",inputReadsFile);
  count=0;
  while ((wordCount = lineFileChop(lf, row)) != 0){
  //while(!feof(inF)){
    //fscanf(inF,"%s\t%d\t%c\%*[^\n]\n",chr_str,&left,&cStrand);
    assert(wordCount>=3);
    chr_str = row[0];
    left = lineFileNeedNum(lf, row, 1);
    cStrand = row[2][0];
    
    khit = kh_get(hashChr_t, hChr , chr_str);
    count++;
    if(kh_exist(hChr,khit)){
      iChr=kh_val(hChr,khit);
      if(cStrand=='+' && (left>0) && (left<=chromSizes[iChr])){
	if(xbl->vec[(iChr<<1)].a[left-1]<255)
	  xbl->vec[(iChr<<1)].a[left-1]++;
      }
      else if(cStrand=='-' && (left>0) && (left+readSize-1<=chromSizes[iChr])){
	if(xbl->vec[(iChr<<1)+1].a[left-1+readSize-1]<255)
	  xbl->vec[(iChr<<1)+1].a[left-1+readSize-1]++;
      }
      else
	verbose(2,"Unrecognized entry: %s,%d,%c\n",chr_str,left,cStrand);
    }
    else
      verbose(2,"Unrecognized entry: %s,%d,%c\n",chr_str,left,cStrand);
  }
  verbose(2,"... completed with %ld reads\n",count);
  
  //carefulClose(&inF);
  lineFileClose(&lf);


  // free mmap. make sure it syncronizes with disk. 
  //  verbose(2,"Syncronizing with disk ...");
  //  msync(xbl->pmmap,xbl->mmapLength,MS_SYNC);
  // verbose(2,"... completed\n");
  //  munmap(xbl->pmmap,xbl->mmapLength);
  *pxbl=xbl;
  return hChr;  
}


void makeXbFileXtractBed(char *chromSizesFile, char *inputReadsFile, char *bedFileName, char *outFileName)
/* makeXbFileXtractBed - Makes a temporal binary file and extracts the regions sorrounding a list of beds. */
{
  khash_t(hashChr_t) *hChr; 
  xbList_t *xbl;
  char *outXbFile="/tmp/tmp.x8b.XXXXXX"; // make tmp file...
  //  int fd;
  
  //if((fd = mkstemp(outXbFile)) < 0)
  //  errnoAbort("Problems opening a temp file\n");
  
  hChr = makeXbFile(chromSizesFile, inputReadsFile, outXbFile, &xbl);


  extractBedFromXbFile(xbl,bedFileName,outFileName);
  

}


int main(int argc, char *argv[])
/* Process command line. */
{
  optionInit(&argc, argv, options);
  if (argc != 5)
    usage();
  
  window = optionInt("window", window);
  readSize = optionInt("readSize", readSize);
  
  makeXbFileXtractBed(argv[1],argv[2],argv[3],argv[4]);
  return 0;
}

