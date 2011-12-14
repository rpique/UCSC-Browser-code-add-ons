/* makeXbCrossMapp - Uses a chainFile to edit target genome mappability from a query genome.. */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"

#include "chain.h"

#include "khash.h"

#include "xbfiles.h"

static char const rcsid[] = "$Id: newProg.c,v 1.30 2010/03/24 21:18:33 hiram Exp $";

void usage()
/* Explain usage and exit. */
{
errAbort(
  "makeXbCrossMapp - Uses a chainFile to edit target genome mappability from a query genome.\n"
  "usage:\n"
  "   makeXbCrossMapp chainFile targetXbFileName queryXbFileName\n"
  "options:\n"
  "   -readSize=20 (default) Important for calculating mappability accross inversions\n"
  "Notes:\n"
  "   1) Writting permission are necessary for both XbFiles\n"
  "   2) The TargetXbFile will be modified the other will not be changed\n"
  );
}

int readsize = 20;   

static struct optionSpec options[] = {
   {"readsize", OPTION_INT},
   {NULL, 0},
};

void makeXbCrossMapp(char *chainFileName,char *targetXbFileName,char *queryXbFileName)
/* makeXbCrossMapp - Uses a chainFile to edit target genome mappability from a query genome.. */
{
  struct lineFile *lf = lineFileOpen(chainFileName, TRUE);

  xbList_t *xbTarget=xbLoadMmap(targetXbFileName);
  xbList_t *xbQuery=xbLoadMmap(queryXbFileName);

  char *row[13];
  //int wordCount;
  struct chain *chain;

  int q,t;
  long int mappTotal=0;
  long int totalBases=0;

  khash_t(hashChr_t) *tHashChr= xbChrNamesKhash(xbTarget); 
  khash_t(hashChr_t) *qHashChr= xbChrNamesKhash(xbQuery); 

  //khash_t(hashChr_t) *qHashChr= kh_init(hashChr_t); 

  
  // AllocVar(chain);  

  /* Parse chain file... see /lib/chain.c chainReadChainLine */
  // while ((wordCount = lineFileChop(lf, row)) != 0){
  while ( (chain = chainReadChainLine(lf)) != NULL){
    int mappReset=0;
    int iChrQ=getChrNumber(qHashChr,chain->qName);
    int iChrT=getChrNumber(tHashChr,chain->tName);
    int offset=0;
    // verbose this
    chainWriteHead(chain,stderr);
    verbose(1,"##Target chr=%s, Query chr=%s\n",xbTarget->names[iChrT],xbQuery->names[iChrQ]);
    // Reading the blocks chainReadBlocks
    /* Now read in block list. */
    t = chain->tStart;
    q = chain->qStart;
    if(chain->qStrand=='-')
      offset= chain->qSize-readsize+1;
    for (;;){
      int j=0;
      int wordCount = lineFileChop(lf, row);
      int size = lineFileNeedNum(lf, row, 0);
      if((chain->qStrand == '-') && (q - chain->qStart < 100)){
	verbose(1,"%s:%i-%i\t%s:%i-%i\n",chain->tName,t,t+size,chain->qName,chain->qSize-(q+size),chain->qSize-q);
      }
      for(j=0;j<size;j++){
	totalBases++;
	if(xbTarget->vec[iChrT].a[t+j]==1){
	  if(chain->qStrand=='+'){
	    if(xbQuery->vec[iChrQ].a[q+j+offset]!=1){
	      xbTarget->vec[iChrT].a[t+j]=0;
	      mappReset++;
	    }else{
	      mappTotal++;
	    }
	  }else{
	    if(xbQuery->vec[iChrQ].a[offset-(q+j)]!=1){
	      xbTarget->vec[iChrT].a[t+j]=0;
	      mappReset++;
	    }else{
	      mappTotal++;
	    }
	  }
	}
      }
      q += size;
      t += size;
      if (wordCount == 1)
        break;
      else if (wordCount < 3)
        errAbort("Expecting 1 or 3 words line %d of %s\n", 
		 lf->lineIx, lf->fileName);
      t += lineFileNeedNum(lf, row, 1); //Adds gap t
      q += lineFileNeedNum(lf, row, 2); //Adds gap q
    }//for Block
    verbose(1,"##Changed %d bases to unmappable\n",mappReset);
    verbose(1,"##Fraction mappable %ld/%ld %g\n",mappTotal,totalBases,(double)mappTotal/(double)totalBases);

    if (q != chain->qEnd)
      errAbort("q end mismatch %d vs %d line %d of %s\n",q, chain->qEnd, lf->lineIx, lf->fileName);
    if (t != chain->tEnd)
      errAbort("t end mismatch %d vs %d line %d of %s\n",t, chain->tEnd, lf->lineIx, lf->fileName);
    //slReverse(&chain->blockList);
    chainFree(&chain);
  }// While Chain  
  xbDestroy(xbTarget);
  xbDestroy(xbQuery);
}

int main(int argc, char *argv[])
/* Process command line. */
{
optionInit(&argc, argv, options);
if (argc != 4)
    usage();
  readsize = optionInt("readsize", readsize);
 makeXbCrossMapp(argv[1],argv[2],argv[3]);
return 0;
}
