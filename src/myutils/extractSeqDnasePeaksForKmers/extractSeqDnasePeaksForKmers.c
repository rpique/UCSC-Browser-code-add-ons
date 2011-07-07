/* extractSeqDnasePeaksForKmers - Extracts sequence from 2bit file from a bin file with DNase counts on each base.. */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"

#include "dnaseq.h"
#include "twoBit.h"

#include <assert.h>
#include <unistd.h>
#include <sys/mman.h> 

// The following is set for hg18 
char *chromosomes[]={"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
  "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
  "chr21","chr22","chrX","chrY"};
int chromSizes[]={247249719,242951149,199501827,191273063,180857866,170899992,
                  158821424,146274826,140273252,135374737,134452384,132349534,
                  114142980,106368585,100338915,88827254 ,78774742 ,76117153 ,
                  63811651 ,62435964 ,46944323 ,49691432 ,154913754,57772954 };
#define NUMCHR 24

typedef unsigned short int myCutSiteIntType; //2 byte
//typedef unsigned char myCutSiteIntType;

static char const rcsid[] = "$Id: newProg.c,v 1.28 2009/07/02 17:14:39 angie Exp $";

void usage()
/* Explain usage and exit. */
{
errAbort(
  "extractSeqDnasePeaksForKmers - Extracts sequence from 2bit file from a bin file with DNase counts on each base.\n"
  "usage:\n"
  "   extractSeqDnasePeaksForKmers input.2bit DnaseCounts.bin output.txt\n"
  "options:\n"
  "   -minReads=200 Minimum number of DNase reads to get a k-mer"
  "   -window=10 Use half of the maximum k-mer for which the produced file will be used \n"
  "   -remMask - Remove Ns and masked genome (not yet implemented)\n"
  );
}

int Window = 10;	/* Window size from command line. */
int minReads= 200;      /* Minimum number of reads to call. */
bool remMask = FALSE;   /* remove N and masked sequence  */


static struct optionSpec options[] = {
   {"remMask", OPTION_BOOLEAN},
   {"window", OPTION_INT},
   {"minReads", OPTION_INT},
   {NULL, 0},
};


void processSeqsFromBinFile(struct twoBitFile *tbf, char *binFileName, FILE *outFile)
/* Get sequences defined by beds.   */
{
  //struct bed *bed, *bedList = bedLoadAll(bedFileName);
  //  char chr_str[100];
  int iChr,j;
  int numCounts;
  int left=0;
  int right=0;
  int LastUp;
  // char cStrand;
  // char buff[500];
  struct dnaSeq *seq;

  FILE *inF; // = mustOpen(binFileName, "r");
  myCutSiteIntType *Counts[NUMCHR];
  myCutSiteIntType *mmapall;

  long int offset;
  double offsetCheck=0;
  int fd,aux;

  // Opening the input memory mapped file //
  aux=0;
  offset=0;
  inF=fopen(binFileName,"r+b");
  assert(inF!=NULL);

  for(iChr=0;iChr<NUMCHR;iChr++){
    offset=offset+chromSizes[iChr]*sizeof(myCutSiteIntType);
    offsetCheck=offsetCheck+chromSizes[iChr]*sizeof(myCutSiteIntType);
  }
  
  fd=fileno(inF);
  fprintf(stderr,"# Mapping %s with %ld=%g bytes in memory\n",binFileName,offset,offsetCheck);
  
  mmapall=(myCutSiteIntType *)mmap(0,offset,PROT_READ,MAP_PRIVATE,fd,0);   
  assert(mmapall!=NULL);

  offset=0;
  for(iChr=0;iChr<NUMCHR;iChr++){
    Counts[iChr]=&(mmapall[offset]);
    offset=offset+chromSizes[iChr];
  }
  fprintf(stderr,"# offset %ld \n",offset);

  // Getting reads
  for(iChr=0;iChr<NUMCHR;iChr++){
    fprintf(stderr,"#Processing %s\n",chromosomes[iChr]);
    LastUp=-1000;
    for(j=Window;j<(chromSizes[iChr]-Window);j++){
      numCounts=Counts[iChr][j];
      if(numCounts>minReads){
	if((j-LastUp)>Window){
	  if(LastUp>0){//Print
	    fprintf(stderr,"%s:%d-%d LastUp%d j%d\n",chromosomes[iChr], left, right, LastUp,j);
	      seq = twoBitReadSeqFrag(tbf, chromosomes[iChr], left, right);
	      writeSeqWithBreaks(outFile, seq->dna, seq->size,right-left);
	      dnaSeqFree(&seq);
	  }
	  left=j-Window;
	}
	right=j+Window;
	LastUp=j;
      }
    }
  }  


}


void extractSeqDnasePeaksForKmers(char *inName, char *binName, char *outName)
/* extractSeqDnasePeaksForKmers - Extracts sequence from 2bit file from a bin file with DNase counts on each base.. */
{
  struct twoBitFile *tbf;
  FILE *outFile = mustOpen(outName, "w");
  struct twoBitSpec *tbs;

  tbs = twoBitSpecNew(inName);
  tbf = twoBitOpen(tbs->fileName);

  processSeqsFromBinFile(tbf, binName, outFile);
 
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
  minReads = optionInt("minReads", minReads);
  remMask = optionExists("remMask");

  dnaUtilOpen();
  extractSeqDnasePeaksForKmers(argv[1], argv[2],argv[3]);
  return 0;
}
