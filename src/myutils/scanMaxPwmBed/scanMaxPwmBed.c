/* scanMaxPwmBed - Get maximum PWM score on a bed region. */
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
  "scanMaxPwmBed - Get maximum PWM score on a bed region\n"
  "usage:\n"
  "   scanMaxPwmBed motif.pwm genome.2bit regions.bed output.file\n"
  "options:\n"
  "   -j - use Jaspar *.pfm format (the default is Dan's TRANSFAC matrix)\n"
  "   -snp - discount most negative base in the pwm aligment if unkown snp hits the score \n"
  "   -window=X - expand the bed region by X base pairs (default 0) \n"
  //  "   -t=6.91 - use threshold cut-off \n"
  "   -base=2.71828183  - use log in that base (default natural log) \n"
  "   -p=0.0 - add p pseudo-counts to frequencies (0.0  automatic) \n"
  "            p=0  If no count is 0 then p=0 if not p=min(count)/2 \n"
  "   -omp=1 - multi-core implementation (default 1 == one thread)"
  "   -mask - remove repeat masked bases (default included) \n"
  );
}

static struct optionSpec options[] = {
   {"j", OPTION_BOOLEAN},
   {"mask", OPTION_BOOLEAN},
   {"snp", OPTION_BOOLEAN},
   //   {"t", OPTION_DOUBLE},
   {"p", OPTION_DOUBLE},
   {"base", OPTION_DOUBLE},
   {"omp" , OPTION_INT},
   {"window", OPTION_INT},
   {NULL, 0},
};


double addPseudoCounts = 0.0;
// double thresh = 6.91 ;
double base = 2.71828183 ;
boolean useJaspar = TRUE;
boolean useSnpRobust = FALSE;
boolean maskRep = FALSE;
int ompNumThreads= 1;
int window = 0;


/* ******************************************************************************** */
/* ******************************************************************************************** */

int scanMaxPwmOneSeq(struct dnaSeq *seq,struct pssm *pwm,struct pssm *rpwm, double *pMaxScore,char *pMaxStrand,int *piMaxScore)
{
  int i,j,Start,End,count;
  char nonATGCbase;
  double score;
  double MaxScore=-1E100;
  char MaxStrand='.';
  int iMaxScore=-1;
  
  count=0;
  
  Start=0;
  End=seq->size-pwm->w+1;
  //if((End-Start+1) < pwm->w)
  //  return 0; //Region too small to fit the motif..
  //#pragma omp parallel for private(nonATGCbase,j,score) reduction(+:count)
  for(i=Start;i<=End;i++){
    nonATGCbase=0;
    for(j=i;j<=(i+pwm->w-1);j++)
      if(ATGCbase(seq->dna[j],maskRep)){
	nonATGCbase++;
	//shift i to skip??
      }
    if(nonATGCbase==0){
      score=compare_subseq_to_pssm(seq->dna+i,pwm,useSnpRobust)/log(base); //Forward strand
      if(score >= MaxScore){
	MaxScore=score;
	iMaxScore=i;
	MaxStrand='+';
      }
      score=compare_subseq_to_pssm(seq->dna+i,rpwm,useSnpRobust)/log(base);  //Reverse strand
      if(score > MaxScore){
	MaxScore=score;
	iMaxScore=i;
	MaxStrand='-';
      }
      count++;
    }
  }
  //	fprintf(out,"%s\t%d\t%d\t-\t%f\n", seq->name,i,i+pwm->w-1,scoreF);
  // fprintf(out,"%f\t%d\t%c\n",MaxScore,iMaxScore,MaxStrand);      
  //return count;
  *pMaxScore=MaxScore;
  *pMaxStrand=MaxStrand;
  *piMaxScore=iMaxScore;
  return iMaxScore;
}

/* ******************************************************************************** */
/* ******************************************************************************** */

void scanMaxPwmBed(char *fileMotif, char *fileSeq, char *fileBed, char *fileOut)
/* scanMaxPwmBed - Get maximum PWM score on a bed region. */
{
  struct lineFile *lf = lineFileOpen(fileBed, TRUE);
  FILE *outFile = mustOpen(fileOut, "w");
  struct twoBitFile *tbf;
  struct twoBitSpec *tbs;
  struct dnaSeq *seq; 

  char *row[10];
  int wordCount;
  char *chr_str;
  int left,right;
  bits32 seqSize;
  // int j;
  
  struct pssm pwm;
  struct pssm rpwm;


  initialise_pssm(&pwm,fileMotif,addPseudoCounts,useJaspar);

  allocateMemoryToMatrix(&rpwm); // Where I do the free!!!
  rpwm.w=pwm.w;
  reversePWM(&rpwm,&pwm);

  tbs = twoBitSpecNew(fileSeq);
  tbf = twoBitOpen(tbs->fileName);

  while ((wordCount = lineFileChop(lf, row)) != 0){
    double MaxScore=-1E100;
    char MaxStrand='.';
    int iMaxScore=-1;

    assert(wordCount>=3);
    chr_str = row[0];
    left = lineFileNeedNum(lf, row, 1)-0;
    right = lineFileNeedNum(lf, row, 2)+1;
    //cStrand = row[5][0];

    //fprintf(outFile,"%s\t%d\t%d\t",chr_str,left,right);
    //for(j=0;j<wordCount;j++)
    //  fprintf(outFile,"%s\t",row[j]);

    seqSize =  twoBitSeqSize(tbf,chr_str);
    if (right > seqSize){// Not valid region...
      fprintf(outFile,"NA");
      fprintf(outFile,"\n");      
      /* errAbort("twoBitReadSeqFrag in %s end (%d) >= seqSize (%d)", name, fragEnd, seqSize); */
    }else{ 
      // Expand window...
      left=max(0,left-window);
      right=min(seqSize,right+window);
      // Extract sequence
      seq = twoBitReadSeqFrag(tbf, chr_str, left, right);      
      //if (noMask) Already taken care of inside scan???
      //toUpperN(seq->dna, seq->size);
      
      scanMaxPwmOneSeq(seq,&pwm,&rpwm,&MaxScore,&MaxStrand,&iMaxScore);
      if (iMaxScore < 0 ){// Not valid region...
	//fprintf(outFile,"NA");
	//fprintf(outFile,"\n");      
	fprintf(outFile,"%s\t%d\t%d\t%s\tNA\tNA\tNA\n",chr_str,left+iMaxScore,left+iMaxScore+pwm.w-1,row[3]);      
      }else{
	fprintf(outFile,"%s\t%d\t%d\t%s\t%f\t%c\t%d\n",chr_str,left+iMaxScore,left+iMaxScore+pwm.w-1,row[3],MaxScore,MaxStrand,iMaxScore-window);      
      }

      
      dnaSeqFree(&seq);
    }
  }
  twoBitSpecFree(&tbs);
  carefulClose(&outFile);
  twoBitClose(&tbf);
}

/* ******************************************************************************** */
/* ******************************************************************************** */

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

  //thresh= optionDouble("t",thresh);
  base= optionDouble("base",base);
  addPseudoCounts= optionDouble("p",addPseudoCounts);
  ompNumThreads= optionInt("omp",ompNumThreads);
  window = optionInt("window",window);

  // omp_set_num_threads(ompNumThreads); //Adjust this to env variable.... or option
  scanMaxPwmBed(argv[1],argv[2],argv[3],argv[4]);
  return 0;
}

