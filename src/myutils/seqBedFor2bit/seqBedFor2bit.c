/* seqBedFor2bit - Extract sequences from 2bit file on Bed regions. */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "dnaseq.h"
#include "twoBit.h"


static char const rcsid[] = "$Id: newProg.c,v 1.28 2009/07/02 17:14:39 angie Exp $";

void usage()
/* Explain usage and exit. */
{
errAbort(
  "seqBedFor2bit - Extract sequences from 2bit file on Bed regions\n"
  "usage:\n"
  "   seqBedFor2bit input.2bit regions.bed output.fa \n"
  "options:\n"
  "   -window=0 Window around the bed\n"
  "   -noMask - convert sequence to all upper case\n"
  );
}

int Window = 0;	/* Window size from command line. */
bool noMask = FALSE;   /* convert seq to upper case */
// char *clBed = NULL;	/* Bed file that specifies bounds of sequences. */


static struct optionSpec options[] = {
   {"noMask", OPTION_BOOLEAN},
   {"window", OPTION_INT},
   {NULL, 0},
};


void processSeqsFromBed(struct twoBitFile *tbf, char *bedFileName, FILE *outFile)
/* Get sequences defined by beds.   */
{
  //struct bed *bed, *bedList = bedLoadAll(bedFileName);
  FILE *inFile = mustOpen(bedFileName, "r");
  char chr_str[100];
  int left;
  bits32 seqSize;
  int right,j;
  char cStrand;
  char buff[500];
  struct dnaSeq *seq;

  while(!feof(inFile)){
    fscanf(inFile,"%s\t%d\t%d\t%*s\t%*s\t%c\%[^\n]\n",chr_str,&left,&right,&cStrand,buff);
    
    left=left-Window;
    right=right+Window+1; // Remember to use the right coordinates. UCSC no +1 

    seqSize =  twoBitSeqSize(tbf,chr_str);
    if (right > seqSize){
      for(j=left;j<right;j++)
	fprintf(outFile,"N");
      fprintf(outFile,"\n");
      
      /* errAbort("twoBitReadSeqFrag in %s end (%d) >= seqSize (%d)", name, fragEnd, seqSize); */
    }else{

      seq = twoBitReadSeqFrag(tbf, chr_str, left, right);
      if (cStrand == '-')
	reverseComplement(seq->dna, seq->size);
      
      if (noMask)
	toUpperN(seq->dna, seq->size);
      
      //Print sequence on outfile
      writeSeqWithBreaks(outFile, seq->dna, seq->size,right-left);
      
      dnaSeqFree(&seq);
    }
  }
}
 
void seqBedFor2bit(char *inName, char *bedName, char *outName)
/* seqBedFor2bit - Extract sequences from 2bit file on Bed regions. */
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
  Window = optionInt("window", Window);
  noMask = optionExists("noMask");
  dnaUtilOpen();
  seqBedFor2bit(argv[1], argv[2],argv[3]);
  return 0;
}

