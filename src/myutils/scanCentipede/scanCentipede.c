/* scanCentipede - Scans a pwm and a lambda profile over a 2bit file and a x8b file. */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"

#include "dnautil.h"
#include "dnaseq.h"
#include "dnaLoad.h"

#include "kvec.h"
#include "khash.h"
#include "xbfiles.h"

#include "pwm.h"

#include <omp.h>

KHASH_MAP_INIT_STR(hashChr_t, unsigned char)


static char const rcsid[] = "$Id: newProg.c,v 1.30 2010/03/24 21:18:33 hiram Exp $";

void usage()
/* Explain usage and exit. */
{
errAbort(
  "scanCentipede - Scans a pwm and a lambda profile over a 2bit file and a x8b file\n"
  "usage:\n"
  "   scanCentipede motif.pwm seq.2bit lambda.txt cuts.x8b \n"
  "options:\n"
  "   -j - use Jaspar *.pfm format (the default is Dan's TRANSFAC matrix)\n"
  "   -snp - discount most negative base in the pwm aligment if unkown snp hits the score \n"
  "   -t=0 - use threshold cut-off on the number of reads for the window \n"
  "   -base=2.71828183  - use log in that base (default natural log) \n"
  "   -p=0.0 - add p pseudo-counts to frequencies (0.0  automatic) \n"
  "            p<=0  If prob is 0 then p=0 if not p=min(count)/2 \n"
  "   -omp=1 - multi-core implementation (default 1 == one thread)"
  "   -mask - remove repeat masked bases \n"
  "   -window = Window around the motif to scan lambda in bp (0 default)\n"
  //  "   -combined -> if set, combine + and - strand into a single output\n"
  );
}

static struct optionSpec options[] = {
   {"j", OPTION_BOOLEAN},
   {"mask", OPTION_BOOLEAN},
   {"snp", OPTION_BOOLEAN},
   {"t", OPTION_DOUBLE},
   {"p", OPTION_DOUBLE},
   {"base", OPTION_DOUBLE},
   {"omp" , OPTION_INT},
   {"window", OPTION_INT},
   {"combined", OPTION_BOOLEAN},
   {NULL, 0},
};

double addPseudoCounts = 0.00;
double thresh = 0 ;
double base = 2.71828183 ;
boolean useJaspar = TRUE;
boolean useSnpRobust =FALSE;
boolean maskRep = FALSE;
int ompNumThreads= 1;

boolean combined=FALSE;
int window = 0;      /* size of the window in which we extend the bed region. */



int scanOne(struct dnaSeq *seq,struct pssm *pwm,struct pssm *rpwm,xbList_t *xbl,int iChr, double *lambda, int winsize)
{
  int i,Start,End;
    
  Start=0+window-1;
  End=seq->size-pwm->w+1-window;
  if((End-Start+1) < pwm->w)
    return 0; //Region too small to fit the motif.. private(j) reduction(+:count)
#pragma omp parallel for  
  for(i=Start;i<=End;i++){
    int j,k; 
    int rowsum=0;
    double multF,multR;
    double scoreF,scoreR;
    int left=i-window;
    int right=i+pwm->w-1+window;      
    int nonATGCbase=0;
    for(j=i;j<=(i+pwm->w-1);j++)
      if(ATGCbase(seq->dna[j],maskRep)){
	nonATGCbase++;
	//shift i to skip??
      }
    if(nonATGCbase==0){
      for(j=left;j<=right;j++)
	rowsum += xbl->vec[(iChr<<1)].a[j];
      for(j=left;j<=right;j++)
	rowsum += xbl->vec[(iChr<<1)+1].a[j];
      if(rowsum>=thresh){
	multF=0;
	k=0;
	for(j=left;j<=right;j++) 
	  multF+=((double)xbl->vec[(iChr<<1)].a[j])*lambda[k++];
	for(j=left;j<=right;j++) 
	  multF+=((double)xbl->vec[(iChr<<1)+1].a[j])*lambda[k++];     
	multF=(multF+log((double)winsize)*rowsum)/log(base);
	scoreF=compare_subseq_to_pssm(seq->dna+i,pwm,useSnpRobust)/log(base); //Forward strand     
	multR=0;
	k=0;
	for(j=right;j>=left;j--)
	  multR+=((double)xbl->vec[(iChr<<1)+1].a[j])*lambda[k++];
	for(j=right;j>=left;j--)
	  multR+=((double)xbl->vec[(iChr<<1)].a[j])*lambda[k++];
	multR=(multR+log((double)winsize)*rowsum)/log(base);
	scoreR=compare_subseq_to_pssm(seq->dna+i,rpwm,useSnpRobust)/log(base);  //Reverse strand
	
	if(multF>=multR)
	  fprintf(stdout,"%s\t%d\t%d\t%d\t%f\t+\t%f\n", seq->name,i,i+pwm->w-1,rowsum,multF,scoreF);
	else
	  fprintf(stdout,"%s\t%d\t%d\t%d\t%f\t-\t%f\n", seq->name,i,i+pwm->w-1,rowsum,multR,scoreR);

	//fprintf(stdout,"%s\t%d\t%d\t%d\t+\t%f\t%f\t-\t%f\t%f\n", seq->name,i,i+pwm->w-1,rowsum,scoreF,multF,scoreR,multR);      
      }
    }   	
  }  

  return 1;
}

int readLambda(char *lambdaFileName,int winsize,double **plambda){
  struct lineFile *lf = lineFileOpen(lambdaFileName, TRUE);
  kvec_t(double) lambda;
  char *row[winsize];
  int j,wordCount;
  
  kv_init_size(double,lambda,winsize);
  
  wordCount = lineFileChop(lf, row);
  for(j=0;j<winsize;j++)
    kv_A(lambda,j)= lineFileNeedDouble(lf,row,j);
  assert(wordCount==winsize); 
  lineFileClose(&lf);
  *plambda=lambda.a;
  return wordCount;
}


/* ******************************************************************************** */

void scanCentipede(char *fileMotif, char *fileSeq, char *lambdaFileName, char *xbFileName)
/* scanCentipede - Scans a pwm and a lambda profile over a 2bit file and a x8b file. */
{
  struct dnaLoad *dl = dnaLoadOpen(fileSeq);
  xbList_t *xbl=xbLoadMmap(xbFileName);
  struct dnaSeq *seq; 

  struct pssm pwm;
  struct pssm rpwm;

  khiter_t k;
  khiter_t khit;
  int hret;
  int iChr,j;

  int winsize;
  double *lambda;
  double aux=0;

  khash_t(hashChr_t) *hChr= kh_init(hashChr_t); 
  
  //Hash chromosome names
  for(iChr=0;iChr<xbl->count;iChr++){
    k = kh_put(hashChr_t, hChr , cloneString(xbl->names[iChr]), &hret);
    assert(hret==1);
    kh_val(hChr, k) = iChr;
  }

  initialise_pssm(&pwm,fileMotif,addPseudoCounts,useJaspar);
  printMatrix(&pwm);
  convertPSSMToLogs(&pwm);  
  printMatrix(&pwm);
  allocateMemoryToMatrix(&rpwm);
  rpwm.w=pwm.w;
  reversePWM(&rpwm,&pwm);

  winsize=(pwm.w*2)+window*4;
  verbose(1,"Reading Lambda(%d):",winsize);
  readLambda(lambdaFileName,winsize,&lambda);
  aux=0;
  for(j=1;j<winsize;j++){
    verbose(2,"%g,",lambda[j]);
    assert((lambda[j]>=0) & (lambda[j]<=1));
    aux=aux+lambda[j];
    lambda[j]=log(lambda[j]);
  }
  verbose(1,"\n");
  assert((aux-1)<1E-7);  

  while ((seq = dnaLoadNext(dl)) != NULL){
    khit = kh_get(hashChr_t, hChr , seq->name);
    if(kh_exist(hChr,khit)){
      iChr=kh_val(hChr,khit);
      //if((left>0) && (left<=xbl->sizes[iChr])){     
      verbose(2, "#\tprocessing: %s\n", seq->name);
      scanOne(seq,&pwm,&rpwm,xbl,iChr,lambda,winsize);   
    }
  }
}






int main(int argc, char *argv[])
/* Process command line. */
{
optionInit(&argc, argv, options);
if (argc != 5)
    usage();

 dnaUtilOpen();

 useJaspar = optionExists("j");
 useSnpRobust = optionExists("snp");
 maskRep = optionExists("mask");

 thresh= optionDouble("t",thresh);
 base= optionDouble("base",base);
 addPseudoCounts= optionDouble("p",addPseudoCounts);
 ompNumThreads= optionInt("omp",ompNumThreads);

 window = optionInt("window", window);

 combined = optionExists("combined");
 
 omp_set_num_threads(ompNumThreads); //Adjust this to env variable.... or option

 // wigOutput = optionExists("wigOutput");

 scanCentipede(argv[1],argv[2],argv[3],argv[4]);
 return 0;
}

