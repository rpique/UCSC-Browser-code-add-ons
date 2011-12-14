/* cutSiteKmerExtract - Extract the kmer surrounding the cut-site. */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"

#include "twoBit.h"
#include "dnautil.h"
#include "dnaseq.h"
#include "dnaLoad.h"

static char const rcsid[] = "$Id: newProg.c,v 1.30 2010/03/24 21:18:33 hiram Exp $";

void usage()
/* Explain usage and exit. */
{
errAbort(
  "cutSiteKmerExtract - Extract the kmer surrounding the cut-site\n"
  "usage:\n"
  "   cutSiteKmerExtract hg18.2bit mappedReads.gz out.txt\n"
  "options:\n"
  "   -readSize = Size of the reads (20bp default) so the minus strand is shifted\n"
  "   -kmerSize= (8 default)\n"
  "   -xxx=XXX\n"
  );
}

int readSize = 20;   
int kmerSize = 8;
int window = 10;


static struct optionSpec options[] = {
   {"kmerSize", OPTION_INT},
   {"readSize", OPTION_INT},
   {NULL, 0},
};

int encodeACGT(char c){
  c=c&0xDF; // 11011111 sets letter to uppercase
  switch(c){
  case 'A':
    return 0;
  case 'C':
    return 1;
  case 'G':
    return 2;
  case 'T':
    return 3;
  default:
    return -1;
  }
}

unsigned long int packKmer(char *s,int k){
  int Base,j;
  unsigned long int kmer=0;
  //  assert((k<<1)<sizeof(unigned long int));
  for(j=0;j<k;j++){
    Base=encodeACGT(s[j]);
    if(Base==-1)
      return 0xFFFFFFFFFFFFFFFF; //
    kmer=(kmer<<2)+Base;
  }
  return kmer;
}

void processLinesFromBed(struct twoBitFile *tbf, char *bedFileName, FILE *oF)
/* Get Cut-site aligments defined by beds.   */
{
  //struct bed *bed, *bedList = bedLoadAll(bedFileName);
  //FILE *inFile = mustOpen(bedFileName, "r");
  struct lineFile *lf = lineFileOpen(bedFileName, TRUE);
  char *row[5];
  int wordCount;
  char *chr_str;

  char cStrand;
  int leftStart;
  // int numGC=0;
  long int count;
 
  bits32 seqSize;
  int right=0,left=0,j;
  // char *adapter="TCGTATGCCGTCTTCTGCTT";


  unsigned long int kmer;
  
  // int hks=(kmerSize>>1);

  struct dnaSeq *seq;


  verbose(2,"Reading in %s ...\n",bedFileName);
  count=0;
  while ((wordCount = lineFileChop(lf, row)) != 0){
    assert(wordCount>=4);
    chr_str = row[0];
    leftStart = lineFileNeedNum(lf, row, 1)-1;
    cStrand = row[2][0];
    
    count++;

    //if(xbl->vec[(iChr<<1)+1].a[left-1+readSize-1]==254)

    if(cStrand=='+'){
      left=leftStart-(kmerSize>>1)-window;
      right=leftStart+readSize-1+kmerSize+window+3;//+(kmerSize>>1)+3;
    }
    else if(cStrand=='-'){
      left=leftStart-kmerSize-3+1-window; //-(kmerSize>>1)-6+1;
      right=leftStart+readSize-1+(kmerSize>>1)+1+window;
    }
    else{
      errAbort("Strand bad format???\n");
    }

    seqSize =  twoBitSeqSize(tbf,chr_str);
    if ((right < seqSize)&&(left>=0)){
      
      seq = twoBitReadSeqFrag(tbf, chr_str, left, right);
      if (cStrand == '-')
	reverseComplement(seq->dna, seq->size);
   
      toUpperN(seq->dna, seq->size);

      kmer=packKmer(seq->dna+window,kmerSize);

      // kmer=packKmer(seq->dna+hks,readSize);
      //fprintf(oF,"\t%d",cStrand,(int)kmer);

      fprintf(oF,"%c\t%d",cStrand,(int)kmer);
      //fprintf(oF,"\t%d\t%s\t%s\t",wordCount,row[3],seq->dna);

      fprintf(oF,"\t%d",(int)packKmer(seq->dna+readSize+window,kmerSize));
      /*
      for(j=0;j<readSize+kmerSize+window*2;j++){
	fprintf(oF,"\t%c",seq->dna[j]);	
      }
      */
      for(j=0;j<(readSize+kmerSize+window*2+2-1);j++){
	fprintf(oF,"\t%d",(int)packKmer(seq->dna+j,5));	
      }
      fprintf(oF,"\n");

      dnaSeqFree(&seq);
    }//if seq
  }//while line
}



void cutSiteKmerExtract(char *inName, char *bedFileName, char *outName)
/* cutSiteKmerExtract - Extract the kmer surrounding the cut-site. */
{
  struct twoBitFile *tbf;
  FILE *outFile = mustOpen(outName, "w");
  struct twoBitSpec *tbs;

  tbs = twoBitSpecNew(inName);
  tbf = twoBitOpen(tbs->fileName);

  processLinesFromBed(tbf, bedFileName, outFile);
 
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
  readSize = optionInt("readSize", readSize);
  kmerSize = optionInt("kmerSize", kmerSize);

  //Window = optionInt("window", Window);
  //noMask = optionExists("noMask");

  dnaUtilOpen();
  cutSiteKmerExtract(argv[1], argv[2],argv[3]);

  return 0;
}


