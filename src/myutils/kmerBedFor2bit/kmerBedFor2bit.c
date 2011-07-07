/* kmerBedFor2bit - Extract kmers at every base of a Bed region. */
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
  "kmerBedFor2bit - Extract kmers at every base of a Bed region\n"
  "usage:\n"
  "   kmerBedFor2bit input.2bit regions.bed output.txt\n"
  "options:\n"
  "options:\n"
  "   -window=0 Window around the bed\n"
  "   -kmerSize= (8 default)\n"
  "   -xxx=XXX\n"
  );
}

int kmerSize=8;
int window=0;

static struct optionSpec options[] = {
   {"kmerSize", OPTION_INT},
   {"window", OPTION_INT},
   {NULL, 0},
};


void extractKmers(FILE *outFile,  struct dnaSeq *seq, char cStrand){
  unsigned int maxKmerNum=(1<<((unsigned long int)kmerSize<<1));
  unsigned long int ForKmer=0,RevKmer=0;
  unsigned int kmer=0;
  unsigned int mask=0xFFFFFFFF;  //0x000FFFFF
  char Base;
  unsigned long int b;
  unsigned int kmerMask=maxKmerNum-1; // (0x00FFFFF)
  unsigned int cutMask=(1<<((unsigned long int)kmerSize))-1; // (0x00FFFFF)
  int j;

  for(j=0;j<seq->size;j++){
	Base=seq->dna[j];
	Base=Base&0xDF; // 11011111 sets letter to uppercase
	mask<<=1;       // This bitMask monitors how many letters are valid 
	if(Base=='A') {
	  b=0;
	}else if(Base=='C'){
	  b=1;
	}else if(Base=='G'){
	  b=2;
	}else if(Base=='T'){
	  b=3;	}else{
	  b=0;
	  mask|=1;
	}
	ForKmer = (ForKmer<<2) + b;
	RevKmer = (RevKmer>>2) + ((((~b) & 0x03)) << ((kmerSize << 1) - 2));

	if(j>=kmerSize){
	  if(!( mask & cutMask )){
	    if (cStrand == '+')
	      kmer=ForKmer&kmerMask;
	    else if (cStrand == '-')
	      kmer=RevKmer&kmerMask;
     	    //fprintf(outFile,"\t%c-%x\n",Base,kmer);
	    fprintf(outFile,"\t%d",kmer);	    
	  }
	  else{
	    //fprintf(outFile,"\t%c-NA %x\n",Base,mask);
	    fprintf(outFile,"\tNA");
	  }
	}
      }
}





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

  int hKmerSize=kmerSize>>1;


  //  verbose(1,"%u %d %x %x\n",maxKmerNum,kmerSize,cutMask,kmerMask);


  while(!feof(inFile)){ // Read region... 
    fscanf(inFile,"%s\t%d\t%d\t%*s\t%*s\t%c\%[^\n]\n",chr_str,&left,&right,&cStrand,buff);
    
    left=left-window;
    right=right+window+1; // Remember to use the right coordinates. UCSC no +1 

    seqSize =  twoBitSeqSize(tbf,chr_str);
    if ((right+hKmerSize > seqSize) || (left-hKmerSize < 0) ){ // Segment off limits... 
      for(j=left;j<right;j++)
	fprintf(outFile,"\tNA");
      fprintf(outFile,"\n");
      
      /* errAbort("twoBitReadSeqFrag in %s end (%d) >= seqSize (%d)", name, fragEnd, seqSize); */
    }else{

      if (cStrand == '+'){ // -1, -1 for forward, 0 0 rev
	seq = twoBitReadSeqFrag(tbf, chr_str, left-hKmerSize-1, right+hKmerSize-1);
      }else{ // +1 +1 for forward, 0 0 reverse...
	seq = twoBitReadSeqFrag(tbf, chr_str, left-hKmerSize+1, right+hKmerSize+1);
	reverseComplement(seq->dna, seq->size);
      }
      extractKmers(outFile,seq,'+');
      dnaSeqFree(&seq);

      if (cStrand == '+'){ // -1, -1 for forward, 0 0 rev
	seq = twoBitReadSeqFrag(tbf, chr_str, left-hKmerSize, right+hKmerSize);
      }else{ // +1 +1 for forward, 0 0 reverse...
	seq = twoBitReadSeqFrag(tbf, chr_str, left-hKmerSize, right+hKmerSize);
	reverseComplement(seq->dna, seq->size);
      }
      extractKmers(outFile,seq,'-');
      dnaSeqFree(&seq);
      
      fprintf(outFile,"\n");


      //      if (noMask)
      //	toUpperN(seq->dna, seq->size);      
      //Print sequence on outfile  
      //      writeSeqWithBreaks(outFile, seq->dna, seq->size,right-left);

      //Loop sequence to get kmers... 

      
      dnaSeqFree(&seq);
    }
  }
}



void kmerBedFor2bit(char *inName, char *bedName, char *outName)
/* kmerBedFor2bit - Extract kmers at every base of a Bed region. */
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
  window = optionInt("window", window);
  kmerSize = optionInt("kmerSize", kmerSize);

  if(kmerSize&0x1)errAbort("kmerSize must be an even number\n");

  dnaUtilOpen();
  kmerBedFor2bit(argv[1], argv[2],argv[3]);
  return 0;
}
