/* makeXbFile - Makes a binary file with chromosome sizes and data type XXX. */
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
  "makeXbFile - Makes a binary file with chromosome sizes and data type XXX\n"
  "usage:\n"
  "   makeXbFile chromSizes.txt input.txt output.xb\n"
  "options:\n"
  //  "   -stranded = Make a stranded xb file (default yes)\n"
  "   -readSize = Size of the reads (20bp default) so the minus strand is shifted\n"
  //  "   -numBytes = Number of bytes used to count reads at each position (default 1), counts up to 255\n"
  "chromSizes.txt File with the size of the chromosomes\n"
  "   chr1\t122143214\n"
  "   chr2\t133978974\n"
  "   ...\n"
  "input.txt File with the read locations\n"
  "   chr1\t10000\t+\n"
  "   chr1\t13123\t-\n"
  "   ...\n"
  "output.xb "
  );
}

int readSize = 20;   

static struct optionSpec options[] = {
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




void makeXbFile(char *chromSizesFile, char *inputReadsFile, char *outXbFile)
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
  long int TotalClippedPositions=0;

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
	if(xbl->vec[(iChr<<1)].a[left-1]==254)
	  TotalClippedPositions++;
	if(xbl->vec[(iChr<<1)].a[left-1]<255)
	  xbl->vec[(iChr<<1)].a[left-1]++;
      }
      else if(cStrand=='-' && (left>0) && (left+readSize-1<=chromSizes[iChr])){
	if(xbl->vec[(iChr<<1)+1].a[left-1+readSize-1]==254)
	  TotalClippedPositions++;
	if(xbl->vec[(iChr<<1)+1].a[left-1+readSize-1]<255)
	  xbl->vec[(iChr<<1)+1].a[left-1+readSize-1]++;
      }
      else
	verbose(2,"Unrecognized entry: %s,%d,%c\n",chr_str,left,cStrand);
    }
    else
      verbose(2,"Unrecognized entry: %s,%d,%c\n",chr_str,left,cStrand);

    if((count%20000000)==0)
      verbose(2,"\n%ldM",count/1000000);
    if((count%1000000)==0)
      verboseDot();
  }
  verbose(2,"\n... completed with %ld reads\n",count);
  verbose(2,"Total clipped positions = %ld\n",TotalClippedPositions);
  
  //carefulClose(&inF);
  lineFileClose(&lf);

  // free mmap. make sure it syncronizes with disk. 
  verbose(2,"Syncronizing with disk ...");
  msync(xbl->pmmap,xbl->mmapLength,MS_SYNC);
  verbose(2,"... completed\n");
  munmap(xbl->pmmap,xbl->mmapLength);
  
}

int main(int argc, char *argv[])
/* Process command line. */
{
  optionInit(&argc, argv, options);
  if (argc != 4)
    usage();
  
  readSize = optionInt("readSize", readSize);
  
  makeXbFile(argv[1],argv[2],argv[3]);
  return 0;
}
