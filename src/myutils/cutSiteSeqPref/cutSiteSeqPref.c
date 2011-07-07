/* cutSiteSeqPref - Finds DNase-I cut-site preference.. */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"

#include "twoBit.h"
#include "dnautil.h"
#include "dnaseq.h"
#include "dnaLoad.h"

#include "khash.h"
#include "xbfiles.h"

static char const rcsid[] = "$Id: newProg.c,v 1.30 2010/03/24 21:18:33 hiram Exp $";

void usage()
/* Explain usage and exit. */
{
errAbort(
  "cutSiteSeqPref - Finds DNase-I cut-site preference.\n"
  "usage:\n"
  "   cutSiteSeqPref cuts.x8b mapp.x8b hg18.2bit \n"
  "options:\n"
  "   -kmerSize= (8 default)\n"
  "   -window= window at each side of the cut to look (100 default)\n"
  "   -ThreshLow= min number of reads of the window (0)\n"
  "   -ThreshUpper= max number of reads of the window (1E100)\n"
  );
}

int kmerSize=8;
int window=100;
double ThreshLow=0.0;
double ThreshUpper=1E100;

static struct optionSpec options[] = {
   {"kmerSize", OPTION_INT},
   {"window", OPTION_INT},
   {"ThreshLow", OPTION_DOUBLE},
   {"ThreshUpper", OPTION_DOUBLE},
   {NULL, 0},
};


unsigned int rev32bitKmer(unsigned int kmer, int kmerSize){
    unsigned int n = (~kmer);
    //n = ((n >> 1) & 0x55555555) | ((n << 1) & 0xaaaaaaaa);
    n = ((n >> 2) & 0x33333333) | ((n << 2) & 0xcccccccc);
    n = ((n >> 4) & 0x0f0f0f0f) | ((n << 4) & 0xf0f0f0f0);
    n = ((n >> 8) & 0x00ff00ff) | ((n << 8) & 0xff00ff00);
    n = ((n >> 16) & 0x0000ffff) | ((n << 16) & 0xffff0000);
    n = n >> (32-(kmerSize<<1));
    return n;
}

void unPackKmer(unsigned long int kmer,int k,char *s){
  int j;
  char ACGT[]={"ACGT"};
  for(j=k-1;j>=0;j--){
    s[j]=ACGT[(kmer & 0x03)];
    kmer=(kmer>>2);
  }
  s[k]='\0';
}


void cutSiteSeqPref(char *xbCutsFileName,char *xbMappFileName, char *twoBitFileName)
/* cutSiteSeqPref - Finds DNase-I cut-site preference.. */
{
  xbList_t *xblCuts=xbLoadMmap(xbCutsFileName);
  xbList_t *xblMapp=xbLoadMmap(xbMappFileName);
  FILE *outFile = stdout; //mustOpen(outName, "w");

  char buff[kmerSize+1];

  unsigned int maxKmerNum=(1<<((unsigned long int)kmerSize<<1));
  int iChr,j,k,ii;

  unsigned int kmerOcc[maxKmerNum];
  unsigned int kmerCutsRev[maxKmerNum];
  unsigned int kmerCutsFor[maxKmerNum];

  


  struct twoBitFile *tbf;
  struct twoBitSpec *tbs;

  struct dnaSeq *seq;

  unsigned long int ForKmer=0,RevKmer=0;
  
  unsigned int mask=0xFFFFFFFF;  //0x000FFFFF
  char Base;
  unsigned long int b;
  unsigned int cutMask=maxKmerNum-1; // (0x00FFFFF)

  int cCutsFor,cMappFor;


  tbs = twoBitSpecNew(twoBitFileName);
  tbf = twoBitOpen(tbs->fileName);

  /*  */
  for(k=0;k<maxKmerNum;++k)
    kmerCutsFor[k]=0;
  for(k=0;k<maxKmerNum;++k)
    kmerCutsRev[k]=0;
  for(k=0;k<maxKmerNum;++k)
    kmerOcc[k]=0;

  

  for(iChr=1;iChr<xblCuts->count;iChr++){
    //   for(iChr=1;iChr<6;iChr++){
    verbose(1,"# Processing %s\n",xblCuts->names[iChr]);
    
    //LOAD SEQUENCE FOR 2BIT FILE
    //seqF = twoBitReadSeqFrag(tbf, xblCuts->names[iChr], left, right);
    //kmerF=packKmer(seqF->dna,kmerSize);
    seq = twoBitReadSeqFragLower(tbf, xblCuts->names[iChr],0,0);
    
    assert( seq->size == xblCuts->sizes[iChr] );
    assert( seq->size == xblMapp->sizes[iChr] );

    cMappFor=0;
    cCutsFor=0;


    for(j=0;j<xblCuts->sizes[iChr];j++){
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
	b=3;
      }else{
	b=0;
	mask|=1;
      }
      ForKmer = (ForKmer<<2) + b;
      RevKmer = (RevKmer>>2) + ((((~b) & 0x03)) << ((kmerSize << 1) - 2));
    
      //  (i + 1 - kmerSize + 1) //Put the cut-point in the middle
      k=j + 1 -(kmerSize>>1);

      // Moving window... 
      if(k==(window)){
	for(ii=(k-window);ii<=(k+window);ii++){
	  if(xblMapp->vec[(iChr)].a[ii]==1){
	    cMappFor++;
	    cCutsFor+= xblCuts->vec[(iChr<<1)].a[ii];
	    //cRev+= xblCuts->vec[(iChr<<1)+1].a[ii];
	    }
	}
      }
      else if((k>window) & (k < xblCuts->sizes[iChr]-window)){
	if(xblMapp->vec[(iChr)].a[k+window]==1){
	    cMappFor++;
	    cCutsFor+= xblCuts->vec[(iChr<<1)].a[k+window];
	}
	if(xblMapp->vec[(iChr)].a[k-window-1]==1){
	  cMappFor--;
	    cCutsFor-= xblCuts->vec[(iChr<<1)].a[k-window-1];
	}
      }


      if(!( mask & cutMask )){
	/* I don't know if the best would be xblMapp at -19,-20,-21 for the - strand? */
	if((xblMapp->vec[(iChr)].a[k]==1) && 
	   (xblMapp->vec[(iChr)].a[k-20]==1) && 
	   ((double)cCutsFor/(double)cMappFor >= ThreshLow) &&
	   ((double)cCutsFor/(double)cMappFor <= ThreshUpper)){
	  int kmer=ForKmer&cutMask;
	  kmerOcc[kmer]++;
	  kmerCutsFor[kmer]+=xblCuts->vec[(iChr<<1)].a[k];
	  kmerCutsRev[kmer]+=xblCuts->vec[(iChr<<1)+1].a[k-1];
	}
      }
    }
    dnaSeqFree(&seq);	
  }
  
  for(k=0;k<maxKmerNum;++k){
    unsigned int n = rev32bitKmer(k,kmerSize); //   (~k);
    n = n & cutMask;
    fprintf(outFile,"%d\t%d",k,n);
    unPackKmer(k,kmerSize,buff);
    fprintf(outFile,"\t%s",buff);
    fprintf(outFile,"\t%d",kmerOcc[k]);//,kmerCutsFor[j][k]=0;
    fprintf(outFile,"\t%d",kmerCutsFor[k]);
    fprintf(outFile,"\t%d",kmerCutsRev[k]);    
    //    for(j=2;j<kmerSize;j=j+2)
    //  fprintf(outFile,"\t%d",(k>>(j))&(cutMask>>(j<<1)));
    
    //Extract central 4-mer .. 
    n = (~n);
    j = (kmerSize-4); 
    //Center = (k>>(j))&(cutMask>>(j<<1)))
    fprintf(outFile,"\t%d", (k>>(j)) & (cutMask >> (j<<1)) );
    for(j=0;j < (kmerSize - 4);j=j+2){
      int R = (k>>j) & 0x3F;
      int L = (n>>j) & 0x3F;
      fprintf(outFile,"\t%d\t%d",L,R); 
    }	             
    fprintf(outFile,"\n");
  }
  



  //  verbose(1,"# Forw reads %d, Revr read %d\n",cF,cR);
  //  extractBedFromXbFile(xbl,bedFileName,outFileName);

}

int main(int argc, char *argv[])
/* Process command line. */
{
  optionInit(&argc, argv, options);
  if (argc != 4)
    usage();

  kmerSize = optionInt("kmerSize", kmerSize);
  window = optionInt("window", window);

  ThreshLow = optionDouble("ThreshLow", ThreshLow);
  ThreshUpper = optionDouble("ThreshUpper", ThreshUpper);

  dnaUtilOpen();


  cutSiteSeqPref(argv[1],argv[2],argv[3]);
  return 0;
}
