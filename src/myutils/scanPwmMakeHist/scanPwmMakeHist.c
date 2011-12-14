/* scanPwmMakeHist - Scan PWM and make a histogram of the scores. */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"

#include "dnautil.h"
#include "dnaseq.h"
#include "dnaLoad.h"
#include "twoBit.h"

#include "pwm.h"
#include <omp.h>


static char const rcsid[] = "$Id: newProg.c,v 1.30 2010/03/24 21:18:33 hiram Exp $";

void usage()
/* Explain usage and exit. */
{
errAbort(
  "scanPwm - scans a sequence for motif matches to a PWM\n"
  "usage:\n"
  "   scanPwm motif.pwm sequence output.txt\n"
  "where:\n"
  "   sequence is a .fa , .nib or .2bit file or a file which is a list of sequence files.\n"
  "options:\n"
  "   -j - use Jaspar *.pfm format (the default is Dan's TRANSFAC matrix)\n"
  "   -snp - discount most negative base in the pwm aligment if unkown snp hits the score \n"
  "   -t=6.91 - use threshold cut-off \n"
  "   -base=2.71828183  - use log in that base (default natural log) \n"
  "   -p=0.0 - add p pseudo-counts to frequencies (0.0  automatic) \n"
  "            p<=0  If prob is 0 then p=0 if not p=min(count)/2 \n"
  "   -omp=1 - multi-core implementation (default 1 == one thread)"
  "   -mask - remove repeat masked bases \n"
  "   -outpwm = pwmfile.dat outputs the pwm used to scan in TRANSFAC format \n"
  "   -bedFile = Scan only regions defined in a bed file genome file has to be 2bit \n"
  "   -bedWindow = 0 (default) Expand bed region if bedFile is used \n" 
  "Histogram options:\n"
  "   -binSize=0.05=1/20 - Size of bins, default 1\n"
  "   -maxBinCount=1000 - Maximum # of bins, default 25\n"
  "   -minVal=-10 - Minimum value to put in histogram, default 0\n"
  "   -xxx=XXX\n"
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
   {"bedWindow" , OPTION_INT},
   {"outpwm", OPTION_STRING},
   {"bedFile", OPTION_STRING},
   {"binSize", OPTION_DOUBLE},
   {"maxBinCount", OPTION_INT},
   {"minVal", OPTION_DOUBLE},
   {NULL, 0},
};

double binSizeR = 0.05;
int maxBinCount = 1000;
double minValR = -10.0;
char * outpwm = (char *) NULL;

char * bedFile = (char *) NULL;
int bedWindow= 0;

double addPseudoCounts = 0.00;
double thresh = 6.91 ;
double base = 2.71828183 ;
boolean useJaspar = TRUE;
boolean useSnpRobust =FALSE;
boolean maskRep = FALSE;
int ompNumThreads= 1;


int scanPwmOneSeqMakeHist(struct dnaSeq *seq,struct pssm *pwm,double *hist)
{
  int i,j,Start,End,count;
  char nonATGCbase;
  double score;
  struct pssm revp,*rpwm;

  rpwm=&revp;
  allocateMemoryToMatrix(rpwm);
  rpwm->w=pwm->w;
  reversePWM(rpwm,pwm);
  
  count=-1;
  
  Start=0;
  End=seq->size-pwm->w+1;
  if((End-Start+1) < pwm->w)
    return 0; //Region too small to fit the motif..
#pragma omp parallel for private(nonATGCbase,j,score) 
  for(i=Start;i<=End;i++){
    int x=0;
    nonATGCbase=0;
    for(j=i;j<=(i+pwm->w-1);j++)
      if(ATGCbase(seq->dna[j],maskRep)){
	nonATGCbase++;
	//shift i to skip??
      }
    if(nonATGCbase==0){
      score=compare_subseq_to_pssm(seq->dna+i,pwm,useSnpRobust)/log(base); //Forward strand
      x = (int) floor((score-minValR) / binSizeR);
      if (x<0) x=0;
      if (x>maxBinCount) x=(maxBinCount-1);
#pragma omp atomic
      hist[x] += 1.0;
      score=compare_subseq_to_pssm(seq->dna+i,rpwm,useSnpRobust)/log(base);  //Reverse strand
      x = (int) floor((score-minValR) / binSizeR);
      if (x<0) x=0;
      if (x>maxBinCount) x=(maxBinCount-1);
#pragma omp atomic
      hist[x] += 1.0;
    } 
  }/* for omp */
  return 0;
}

/* ******************************************************************************** */

int scanPwmOneShortSeqMakeHist(struct dnaSeq *seq,struct pssm *pwm, struct pssm *rpwm, double *hist)
{
  int i,j,Start,End;
  char nonATGCbase;
  double score;
  
  Start=0;
  End=seq->size-pwm->w+1;
  if((End-Start+1) < pwm->w)
    return 0; //Region too small to fit the motif..
  for(i=Start;i<=End;i++){
    int x=0;
    nonATGCbase=0;
    for(j=i;j<=(i+pwm->w-1);j++)
      if(ATGCbase(seq->dna[j],maskRep)){
	nonATGCbase++;
	//shift i to skip??
      }
    if(nonATGCbase==0){
      score=compare_subseq_to_pssm(seq->dna+i,pwm,useSnpRobust)/log(base); //Forward strand
      x = (int) floor((score-minValR) / binSizeR);
      if (x<0) x=0;
      if (x>maxBinCount) x=(maxBinCount-1);
      hist[x] += 1.0;
      score=compare_subseq_to_pssm(seq->dna+i,rpwm,useSnpRobust)/log(base);  //Reverse strand
      x = (int) floor((score-minValR) / binSizeR);
      if (x<0) x=0;
      if (x>maxBinCount) x=(maxBinCount-1);
      hist[x] += 1.0;
    } 
  }
  return 0;
}

/* ******************************************************************************** */

void scanPwmBedMakeHist(char *fileMotif, char *fileSeq, char *fileBed)
/* scanMaxPwmBedMakeHist - Get maximum PWM score on a bed region. */
{
  struct lineFile *lf = lineFileOpen(fileBed, TRUE);
  //  FILE *outFile = mustOpen(fileOut, "w");
  struct twoBitFile *tbf;
  struct twoBitSpec *tbs;
  struct dnaSeq *seq; 

  char *row[10];
  int wordCount;
  char *chr_str;
  int left,right;
  bits32 seqSize;
  // int j;

  int i;
  double binStartR;
  double *hist = NULL;
  
  struct pssm pwm;
  struct pssm rpwm;
  FILE *pwmfd;

  AllocArray(hist, maxBinCount);

  initialise_pssm(&pwm,fileMotif,addPseudoCounts,useJaspar);
  if(!sameString(outpwm,"0")){
    pwmfd=mustOpen(outpwm,"w");
    fprintMatrix(pwmfd,&pwm);
    carefulClose(&pwmfd);
  }
  printMatrix(&pwm);
  convertPSSMToLogs(&pwm);  
  printMatrix(&pwm);

  allocateMemoryToMatrix(&rpwm); // Where I do the free!!!
  rpwm.w=pwm.w;
  reversePWM(&rpwm,&pwm);

  tbs = twoBitSpecNew(fileSeq);
  tbf = twoBitOpen(tbs->fileName);

  while ((wordCount = lineFileChop(lf, row)) != 0){
    assert(wordCount>=3);
    chr_str = row[0];
    left = lineFileNeedNum(lf, row, 1)-0;
    right = lineFileNeedNum(lf, row, 2)+1;

    seqSize =  twoBitSeqSize(tbf,chr_str);
    if (right <= seqSize){
      // Expand window...
      left=max(0,left-bedWindow);
      right=min(seqSize,right+bedWindow);
      // Extract sequence
      seq = twoBitReadSeqFrag(tbf, chr_str, left, right);      
      scanPwmOneShortSeqMakeHist(seq,&pwm,&rpwm,hist);
      dnaSeqFree(&seq);
    }
  }
  twoBitSpecFree(&tbs);
  //carefulClose(&outFile);
  twoBitClose(&tbf);
  for (i=0; i<=maxBinCount; ++i){
    binStartR = i*binSizeR + minValR;
    fprintf(stdout,"%3d\t%g\t%g\t%g\n", i, binStartR, binStartR+binSizeR, hist[i]);
  }
}

/* ******************************************************************************** */
/* ******************************************************************************** */

void scanPwmMakeHist(char *fileMotif, char *fileSeq)
/* scanPwmMakeHist - Scan PWM and make a histogram of the scores. */
{
  struct dnaLoad *dl = dnaLoadOpen(fileSeq);
  struct dnaSeq *seq; 
  int i;
  double binStartR;
  double *hist = NULL;

  struct pssm pwm;
  FILE *pwmfd;
  
  initialise_pssm(&pwm,fileMotif,addPseudoCounts,useJaspar);
  if(!sameString(outpwm,"0")){
    pwmfd=mustOpen(outpwm,"w");
    fprintMatrix(pwmfd,&pwm);
    carefulClose(&pwmfd);
  }
  printMatrix(&pwm);
  convertPSSMToLogs(&pwm);  
  printMatrix(&pwm);

  AllocArray(hist, maxBinCount);
 

  //motifList = fixMotifs(dnaMotifLoadAll(motifFile));

  while ((seq = dnaLoadNext(dl)) != NULL)
    {
      verbose(2, "#\tprocessing: %s\n", seq->name);
      scanPwmOneSeqMakeHist(seq,&pwm,hist);   
    }

  for (i=0; i<=maxBinCount; ++i){
    binStartR = i*binSizeR + minValR;
    fprintf(stdout,"%3d\t%g\t%g\t%g\n", i, binStartR, binStartR+binSizeR, hist[i]);
  }
}



int main(int argc, char *argv[])
/* Process command line. */
{
optionInit(&argc, argv, options);
if (argc != 3)
    usage();

 dnaUtilOpen();

//motif = optionVal("motif", NULL);
//chr = optionVal("chr", NULL);
//strand = optionVal("strand", NULL);
 useJaspar = optionExists("j");
 useSnpRobust = optionExists("snp");
 maskRep = optionExists("mask");

 outpwm = optionVal("outpwm", "0");
 bedFile = optionVal("bedFile", "0");

 thresh= optionDouble("t",thresh);
 base= optionDouble("base",base);
 addPseudoCounts= optionDouble("p",addPseudoCounts);
 ompNumThreads= optionInt("omp",ompNumThreads);

 minValR = optionDouble("t",minValR);
 maxBinCount = optionInt("maxBinCount", maxBinCount);
 binSizeR = optionDouble("t",binSizeR);

 bedWindow = optionInt("bedWindow",bedWindow);


 if (binSizeR <= 0.0)
   errAbort("invalid binSize, must be greater than zero: %g\n", binSizeR);

 omp_set_num_threads(ompNumThreads); //Adjust this to env variable.... or option

 if(!sameString(bedFile,"0"))
   scanPwmBedMakeHist(argv[1],argv[2],bedFile);
 else
   scanPwmMakeHist(argv[1],argv[2]);

 return 0;
}

