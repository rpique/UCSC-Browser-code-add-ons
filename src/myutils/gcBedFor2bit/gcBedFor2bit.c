/* gcBedFor2bit - Calculate g/c percent for regions covered by a bed. */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"

#include "dnaseq.h"
#include "twoBit.h"


static char const rcsid[] = "$Id: newProg.c,v 1.30 2010/03/24 21:18:33 hiram Exp $";

void usage()
/* Explain usage and exit. */
{
errAbort(
  "gcBedFor2bit - Calculate g/c percent for regions covered by a bed\n"
  "usage:\n"
  "   gcBedFor2bit seq.2bit regions.bed output.bed\n"
  "options:\n"
  "   -verbose=1 Default verbose level=1\n"
  "details:\n"
  "   Only the first 3 columns are necessary. Returns on the output the\n"
  "   first 3 columns defining the bed regions plus the GC percent of the region.\n"
  );
}

static struct optionSpec options[] = {
   {NULL, 0},
};



void processSeqsFromBed(struct twoBitFile *tbf, char *bedFileName, FILE *outFile)
/* Get sequences defined by beds.   */
{
  //struct bed *bed, *bedList = bedLoadAll(bedFileName);
  //FILE *inFile = mustOpen(bedFileName, "r");
  struct lineFile *lf = lineFileOpen(bedFileName, TRUE);

  char *row[10];
  int wordCount;
  char *chr_str;


  //  char chr_str[100];
  int left;
  bits32 seqSize;
  int right,j;
  //  char cStrand;
  //  char buff[500];
  struct dnaSeq *seq;

  int countGC=0;
  int countACGT=0;

  //while(!feof(inFile)){
  //fscanf(inFile,"%s\t%d\t%d\t%*s\t%*s\t%c\%[^\n]\n",chr_str,&left,&right,&cStrand,buff);
  while ((wordCount = lineFileChop(lf, row)) != 0){
    assert(wordCount>=3);
    chr_str = row[0];
    left = lineFileNeedNum(lf, row, 1);
    right = lineFileNeedNum(lf, row, 2);
    //cStrand = row[5][0];
    
    //    left=left-Window;
    // right=right+Window+1; // Remember to use the right coordinates. UCSC no +1 

    fprintf(outFile,"%s\t%d\t%d\t",chr_str,left,right);
    
    seqSize =  twoBitSeqSize(tbf,chr_str);
    if (right > seqSize){// Not valid region...
      fprintf(outFile,"NA");
      fprintf(outFile,"\n");      
      /* errAbort("twoBitReadSeqFrag in %s end (%d) >= seqSize (%d)", name, fragEnd, seqSize); */
    }else{      
      seq = twoBitReadSeqFrag(tbf, chr_str, left, right);      
      //if (noMask)
      toUpperN(seq->dna, seq->size);
      
      for(j=0;j<seq->size;++j){
	if((seq->dna[j]=='C')||(seq->dna[j]=='G')){
	  countGC++;
	  countACGT++;
	}
	else if((seq->dna[j]=='A')||(seq->dna[j]=='T')){
	  countACGT++;
	}
      }
      if(countACGT>0)
	fprintf(outFile,"%f\n",(double)countGC/(double)countACGT*100.00);
      else 
	fprintf(outFile,"NA\n");
	
      //Print GC content on the outfile on outfile
      
      dnaSeqFree(&seq);
    }
  }
}


void gcBedFor2bit(char *inName, char *bedName, char *outName)
/* gcBedFor2bit - Calculate g/c percent for regions covered by a bed. */
{
  struct twoBitFile *tbf;
  FILE *outFile = mustOpen(outName, "w");
  struct twoBitSpec *tbs;

  tbs = twoBitSpecNew(inName);
  tbf = twoBitOpen(tbs->fileName);
  
  processSeqsFromBed(tbf, bedName, outFile);
  
  twoBitSpecFree(&tbs);
  carefulClose(&outFile);
  twoBitClose(&tbf);
}

int main(int argc, char *argv[])
/* Process command line. */
{
  optionInit(&argc, argv, options);
  if (argc != 4)
    usage();
  //  Window = optionInt("window", Window);
  //  noMask = optionExists("noMask");
  dnaUtilOpen();
  gcBedFor2bit(argv[1], argv[2],argv[3]);
  return 0;
}
