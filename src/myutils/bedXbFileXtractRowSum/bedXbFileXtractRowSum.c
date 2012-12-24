/* bedXbFileXtractRowSum - Extracts/row-summarizes the data on a Xbfile using the bed defined regions. */
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
	   "bedXbFileXtractRowSum - Extracts/row-summarizes the data on a Xbfile using the bed defined regions\n"
	   "usage:\n"
	   "   bedXbFileXtractRowSum file.x8b regions.bed output.txt\n"
	   "   supports stdin for regions.bed and stdout for output.bed\n"
	   "options:\n"
	   "   -window = Window around the motif to expand in bp (0 default)\n"
	   "   -shiftF = Shift window on the Forward (upper) strand in BP to the left (0 default)\n"
	   "   -shiftR = Shift window on the Reverse (bottom) strand in BP to the right (0 default)\n"
	   "    shiftR should only be used if xb isStranded\n"
	   "   -mappFile = Mappability x8b file \n"
	   "   -readSize = 20 (default) Shift for the reverse strand on mapp. \n" 
	   "regions.bed is a tab delimited file with at least four columns: \n"
	   "   Chr1\t1000100\t1000150\t+\t[Other Columns]\n"
	   "   .... \n"
	   "   ChrY\t1000100\t1000150\t-\t[Other Columns]\n"
	   "output.txt is in the following form: \n"
	   "   Instance1: Sum Of Forward reads (2*Window+MotifLen columns) and Reverse strand reads (2*Window+MotifLen columns).\n"
	   " .... \n\n"    
	   );
}

int window = 0;      /* size of the window in which we extend the bed region. */
int shiftR = 0;   
int shiftF = 0;   

char * mappFile = (char *) NULL;
int readSize= 0;
int useMappFile = FALSE;



static struct optionSpec options[] = {
   {"window", OPTION_INT},
   {"shiftF", OPTION_INT},
   {"shiftR", OPTION_INT},
   {"mappFile", OPTION_STRING},
   {"readSize", OPTION_INT},
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
  int rowsum;
  //  long int offset;
  //  double offsetCheck;
  //char cStrand;

  char *row[10];
  int wordCount;
  char *chr_str;
  
  //Hash chromosome names? instead of strcmp??  
  cF=0;cR=0;
  while ((wordCount = lineFileChop(lf, row)) != 0){
  //while(!feof(bedF)){
    //fscanf(bedF,"%s\t%d\t%d\t%c\%[^\n]\n",chr_str,&left,&right,&cStrand,buff);
    assert(wordCount>=3);
    chr_str = row[0];
    left = lineFileNeedNum(lf, row, 1);
    right = lineFileNeedNum(lf, row, 2);
    // cStrand = row[5][0];

    // I should use the hash instead
    for(iChr=0;iChr<xbl->count;iChr++)
      if(strcmp(xbl->names[iChr],chr_str)==0)
	break;
    if(iChr<xbl->count){      
      rowsum=0;
      left=left-window;
      right=right+window;
      
      if((left-shiftF >= 0) && (right+shiftR < (xbl->sizes[iChr]))){
	if(xbl->isStranded){
	  for(j=left-shiftF;j<=right-shiftF;j++)
	    rowsum += xbl->vec[(iChr<<1)].a[j];
	  for(j=left+shiftR;j<=right+shiftR;j++)
	    rowsum += xbl->vec[(iChr<<1)+1].a[j];
	  fprintf(outF,"%d\n",rowsum);
	}else{
	  for(j=left-shiftF;j<=right-shiftF;j++)
	    rowsum += xbl->vec[iChr].a[j];	  
	}
      }else{
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
  verbose(1,"# Succesfully processed %d segments\n",count);
}



void extractBedFromXbFileWithMapp(xbList_t *xbl,xbList_t *xblMapp, char *bedFileName, char *outFileName)
/* extractBedFromXbFile - Extract the region sorrounding a list of beds from an Xb binary file. */
{
  //FILE *bedF=mustOpen(bedFileName,"r");
  struct lineFile *lf = lineFileOpen(bedFileName, TRUE);

  FILE *outF=mustOpen(outFileName,"w");

  int skipLine=0;
  long int cF,cR;

  int iChr,j;
  int left,right;
  int count=0;
  char cStrand;
  int rowSum;
  int rowSumMapp;


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
      rowSum=0;
      rowSumMapp=0;
      left=left-window;
      right=right+window;
      
      if((left >= 0) && (right < (xbl->sizes[iChr]))){
	if(cStrand=='+'){	  
	  for(j=left;j<=right;j++)
	    if(xblMapp->vec[iChr].a[j]==1){
	      rowSum+=xbl->vec[(iChr<<1)].a[j];
	      rowSumMapp++;
	    }
	  for(j=left;j<=right;j++)
	    if(xblMapp->vec[iChr].a[j-readSize+1]==1){
	      rowSum+=xbl->vec[(iChr<<1)+1].a[j];
	      rowSumMapp++;
	    }
	  fprintf(outF,"%d\t%d\n",rowSum,rowSumMapp);	  
	}else if(cStrand=='-'){
	  for(j=right;j>=left;j--)
	    if(xblMapp->vec[iChr].a[j-readSize+1]==1){
	      rowSum+=xbl->vec[(iChr<<1)+1].a[j];
	      rowSumMapp++;
	    }
	  for(j=right;j>=left;j--)
	    if(xblMapp->vec[iChr].a[j]==1){
	      rowSum+=xbl->vec[(iChr<<1)].a[j];	    
	      rowSumMapp++;
	    }	  
	  fprintf(outF,"%d\t%d\n",rowSum,rowSumMapp);	  
	}
	else{
	  verbose(1,"# Unrecognized strand parameter!!!\n");	
	  skipLine=1;
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
      fprintf(outF,"NA\t0\n");
      skipLine=0;
    }    
    count++;
  }//end while(!feof(bedF)){

  //  carefulClose(&bedF);
  lineFileClose(&lf);

  carefulClose(&outF);
  //  verbose(1,"# Mapped %ld on F, %ld on R strands \n",cF,cR);
  verbose(1,"# Succesfully processed %d segments\n",count);
}




void bedXbFileXtractRowSum(char *xbFileName, char *bedFileName, char *outFileName)
/* bedXbFileXtractRowSum - Extracts/row-summarizes the data on a Xbfile using the bed defined regions. */
{
  xbList_t *xbl;
  xbList_t *xblMapp;

  xbl=xbLoadMmap(xbFileName);  

  if(sameString(mappFile,"0")){    
    extractBedFromXbFile(xbl,bedFileName,outFileName);
  }else{
    xblMapp=xbLoadMmap(mappFile);  
    extractBedFromXbFileWithMapp(xbl,xblMapp,bedFileName,outFileName);
  }
}

int main(int argc, char *argv[])
/* Process command line. */
{
  optionInit(&argc, argv, options);
  if (argc != 4)
    usage();

  window = optionInt("window", window);
  shiftR = optionInt("shiftR", shiftR);
  shiftF = optionInt("shiftF", shiftF);

  mappFile = optionVal("mappFile", "0");
  readSize = optionInt("readSize", readSize);
  
  bedXbFileXtractRowSum(argv[1],argv[2],argv[3]);
  return 0;
}



