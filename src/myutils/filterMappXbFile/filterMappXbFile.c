/* filterMappXbFile - Uses a mappability Xb file and outputs the mappability of the aligment. */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"

#include "khash.h"
#include "xbfiles.h"

//KHASH_MAP_INIT_STR(hashChr_t, unsigned char)

static char const rcsid[] = "$Id: newProg.c,v 1.30 2010/03/24 21:18:33 hiram Exp $";

void usage()
/* Explain usage and exit. */
{
errAbort(
  "filterMappXbFile - Uses a mappability Xb file and outputs the mappability of the aligment\n"
  "usage:\n"
  "   filterMappXbFile mapp.x8b aligment.bed output.txt\n"
  " * supports stdin for regions.bed and stdout for output.bed\n"
  "   \n"
  "aligments.bed is a tab delimited file with at least three columns: \n"
  "   Chr1\t1000100\t+\t[Other Columns]\n"
  "   .... \n"
  "   ChrY\t1000100\t+\t[Other Columns]\n"
  "   1-based coordinates are assumed \n"  
  "options:\n"
  "   -verbose=XXX\n"
  );
}

static struct optionSpec options[] = {
   {NULL, 0},
};


void filterMappXbFile(char *xbFileName, char *bedFileName, char *outFileName)
/* filterMappXbFile - Uses a mappability Xb file and outputs the mappability of the aligment. */
{
  struct lineFile *lf = lineFileOpen(bedFileName, TRUE);

  FILE *outF=mustOpen(outFileName,"w");

  xbList_t *xbl=xbLoadMmap(xbFileName);

  int skipLine=0;
  //long int cF,cR;

  int iChr,j;
  int left;//,right;
  int count=0;
  char cStrand;

  char *row[10];
  int wordCount;
  char *chr_str;

  khiter_t k;
  khiter_t khit;
  int hret;

  khash_t(hashChr_t) *hChr= kh_init(hashChr_t); 
  
  //Hash chromosome names
  for(iChr=0;iChr<xbl->count;iChr++){
    k = kh_put(hashChr_t, hChr , cloneString(xbl->names[iChr]), &hret);
    assert(hret==1);
    kh_val(hChr, k) = iChr;
  }
  
  // cF=0;cR=0;
  while ((wordCount = lineFileChop(lf, row)) != 0){
    assert(wordCount>=3);
    chr_str = row[0];
    left = lineFileNeedNum(lf, row, 1);
    //    right = lineFileNeedNum(lf, row, 2);
    cStrand = row[2][0];
    
    for(j=0;j<wordCount;j++)
      fprintf(outF,"%s\t",row[j]);
    
    /*    for(iChr=0;iChr<xbl->count;iChr++)
      if(strcmp(xbl->names[iChr],chr_str)==0)
	break;
	if(iChr<xbl->count){
    */
    /* SHOULD I DO SOMETHING ABOUT STRAND?*/
    /* No; because I use the first position or the forward even if the match is in the reverse. */
    khit = kh_get(hashChr_t, hChr , chr_str);
    if(kh_exist(hChr,khit)){
      iChr=kh_val(hChr,khit);
      if((left>0) && (left<=xbl->sizes[iChr])){      
	fprintf(outF,"%d\n",xbl->vec[iChr].a[left-1]);
      }else{
	verbose(2,"# Skipping segment off-limits\n");
	//verbose(1,"# Unreq %s:%d-%d in line %d !!!\n",chr_str,left,right,count);
	skipLine=1;
      }
    }else{
      verbose(2,"# Skipping segment off-limits\n");
      //verbose(1,"# Unreq %s:%d-%d in line %d !!!\n",chr_str,left,right,count);
      skipLine=1;
    }             
    if(verboseLevel()==2){
      if((count%10000)==0){
	verboseDot();
	if((count%800000)==0)
	  verbose(2,"\n");
      }
    }
    if(skipLine==1){
      fprintf(outF,"NA\n");
      skipLine=0;
    }    
    count++;
  }
  
  lineFileClose(&lf);

  carefulClose(&outF);
  //verbose(1,"# Mapped %ld on F, %ld on R strands \n",cF,cR);
  verbose(1,"# Succesfully processed %d segments\n",count);

  //free xbl?
}


int main(int argc, char *argv[])
/* Process command line. */
{
  optionInit(&argc, argv, options);
  if (argc != 4)
    usage();

  //  window = optionInt("window", window);
  //  readsize = optionInt("readsize", readsize);
  
  filterMappXbFile(argv[1],argv[2],argv[3]);
  return 0;
}


