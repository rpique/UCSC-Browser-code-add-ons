/* scanPwm - scans a sequence for motif matches to a PWM. */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "dnautil.h"
#include "dnaseq.h"
#include "dnaLoad.h"

#include "pwm.h"

#include <omp.h>


static char const rcsid[] = "$Id: newProg.c,v 1.28 2009/07/02 17:14:39 angie Exp $";

void usage()
/* Explain usage and exit. */
{
errAbort(
  "scanPwm - scans a sequence for motif matches to a PWM\n"
  "usage:\n"
  "   scanPwm motif.pwm sequence\n"
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
   {NULL, 0},
};

double addPseudoCounts = 0.00;
double thresh = 6.91 ;
double base = 2.71828183 ;
boolean useJaspar = TRUE;
boolean useSnpRobust =FALSE;
boolean maskRep = FALSE;
int ompNumThreads= 1;

int scanPwmOneSeq(struct dnaSeq *seq,struct pssm *pwm)
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
#pragma omp parallel for private(nonATGCbase,j,score) reduction(+:count)
  for(i=Start;i<=End;i++){
    nonATGCbase=0;
    for(j=i;j<=(i+pwm->w-1);j++)
      if(ATGCbase(seq->dna[j],maskRep)){
	nonATGCbase++;
	//shift i to skip??
      }
    if(nonATGCbase==0){
      score=compare_subseq_to_pssm(seq->dna+i,pwm,useSnpRobust)/log(base); //Forward strand
      if(score >= thresh){
	fprintf(stdout,"%s\t%d\t%d\t+\t%f\n", seq->name,i,i+pwm->w-1,score);
	count++;
      }
      score=compare_subseq_to_pssm(seq->dna+i,rpwm,useSnpRobust)/log(base);  //Reverse strand
      if(score >= thresh){
	fprintf(stdout,"%s\t%d\t%d\t-\t%f\n", seq->name,i,i+pwm->w-1,score);
	count++;
      }   	
    }  
  }

  return count;
}




/* ******************************************************************************** */

void scanPwm(char *fileMotif, char *fileSeq)
/* scanPwm - scans a sequence for motif matches to a PWM. */
{
  struct dnaLoad *dl = dnaLoadOpen(fileSeq);
  struct dnaSeq *seq; 

  struct pssm pwm;

  initialise_pssm(&pwm,fileMotif,addPseudoCounts,useJaspar);
 
  //motifList = fixMotifs(dnaMotifLoadAll(motifFile));

  while ((seq = dnaLoadNext(dl)) != NULL)
    {
      verbose(2, "#\tprocessing: %s\n", seq->name);
      scanPwmOneSeq(seq,&pwm);   
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

 thresh= optionDouble("t",thresh);
 base= optionDouble("base",base);
 addPseudoCounts= optionDouble("p",addPseudoCounts);
 ompNumThreads= optionInt("omp",ompNumThreads);


 omp_set_num_threads(ompNumThreads); //Adjust this to env variable.... or option



// wigOutput = optionExists("wigOutput");

 scanPwm(argv[1],argv[2]);
return 0;
}
