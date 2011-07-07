/* bedXbFileXtract - Extracts the data on a Xbfile using the bed defined regions. */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"

#include "khash.h"
#include "xbfiles.h"

static char const rcsid[] = "$Id: newProg.c,v 1.30 2010/03/24 21:18:33 hiram Exp $";

void usage()
/* Explain usage and exit. */
{
  errAbort(
	   "bedXbFileXtract - Extracts the data on a Xbfile using the bed defined regions\n"
	   "usage:\n"
	   "   bedXbFileXtract file.x8b regions.bed output.txt\n"
	   "   supports stdin for regions.bed and stdout for output.bed\n"
	   "options:\n"
	   "   -window = Window around the motif to expand in bp (0 default)\n"
	   //	   "   -shiftR = X  shift Reverse strand X bp (0 default), only used if xb isStranded\n"
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
int shiftR = 0;   

static struct optionSpec options[] = {
   {"window", OPTION_INT},
   {"shiftR", OPTION_INT},
   {NULL, 0},
};

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
}

void bedXbFileXtract(char *xbFileName, char *bedFileName, char *outFileName)
/* bedXbFileXtract - Extracts the data on a Xbfile using the bed defined regions. */
{
  xbList_t *xbl;

  xbl=xbLoadMmap(xbFileName);  
  extractBedFromXbFile(xbl,bedFileName,outFileName);

}

int main(int argc, char *argv[])
/* Process command line. */
{
  optionInit(&argc, argv, options);
  if (argc != 4)
    usage();

  window = optionInt("window", window);
  shiftR = optionInt("shiftR", shiftR);
  
  bedXbFileXtract(argv[1],argv[2],argv[3]);
  return 0;
}

