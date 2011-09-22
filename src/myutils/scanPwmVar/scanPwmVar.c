/* scanPwmVar - Scan binding sites hit by SNPs. */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"

#include "dnautil.h"
#include "dnaseq.h"
#include "dnaLoad.h"

#include "pwm.h"

#include <omp.h>


static char const rcsid[] = "$Id: newProg.c,v 1.30 2010/03/24 21:18:33 hiram Exp $";

void usage()
/* Explain usage and exit. */
{
errAbort(
  "scanPwmVar - Scan binding sites hit by SNPs\n"
  "usage:\n"
  "   scanPwmVar motif.pwm sequence.2bit snpfile \n"
  "where:\n"
  "   sequence is a .fa , .nib or .2bit file or a file which is a list of sequence files.\n"
  "options:\n"
  "   -j - use Jaspar *.pfm format (the default is Dan's TRANSFAC matrix)\n"
  "   -snp - discount most negative base in the pwm aligment if unkown snp hits the score \n"
  "   -t=6.91 - use threshold cut-off \n"
  "   -base=2.71828183  - use log in that base (default natural log) \n"
  "   -p=0.0 - add p pseudo-counts to frequencies (-1 default automatic) \n"
  "            If no count is 0 then go ahed if not p=min(count)/2       \n"
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


double addPseudoCounts = 0.0;
double thresh = 6.91 ;
double base = 2.71828183 ;
boolean useJaspar = TRUE;
boolean useSnpRobust =FALSE;
boolean maskRep = FALSE;
int ompNumThreads= 1;

/* ******************************************************************************** */
/* ******************************************************************************** */
/* ******************************************************************************** */

#pragma pack(push)  /* push current alignment to stack */
#pragma pack(1)     /* set alignment to 1 byte boundary */
typedef struct{
  unsigned char chr; // index to the sequence list, 1-254, 0 (reserved), 255 (Multiple ocurring times)
  int pos;  //Position + for the plus strand, negative numbers for the negative strand. 
            //if chr=255, pos is the index to the array holding all occurrences...
}chrpos_t;

typedef struct{
  int pos;
  char ref;
  char alt;
}snp_t;

typedef struct{
  int pos;          // Position on the genome of the firt base of the indel. 0-based...
  short int reflen; // Length of the reference sequence
  short int altlen; // Length of the alternate sequence
  char *ref;        // Reference sequence
  char *alt;        // Alternate sequence to replace.  // We will use the convention that the first base of the Reference should stay.
}indel_t;

#pragma pack(pop)   /* restore original alignment from stack */


#include "khash.h"

//KHASH_MAP_INIT_INT(32, int)

KHASH_MAP_INIT_INT(hCount_t, int)            // Hash number of times read seen. 
KHASH_MAP_INIT_INT(hashPos_t, chrpos_t)      // 
KHASH_MAP_INIT_STR(hashChr_t, unsigned char)

#include "kvec.h"

typedef kvec_t(chrpos_t)repvec_t;  //Repeats of a kmer
typedef kvec_t(repvec_t)replist_t; //List of all kmer repeat vectors
typedef kvec_t(snp_t)snpvec_t;     //Vector of snps in a chromosome
typedef kvec_t(snpvec_t)snplist_t;    //List of snps vectors.  	
typedef kvec_t(indel_t)indelvec_t;      //Vector of indels in a chromosome
typedef kvec_t(indelvec_t)indellist_t;    //List of indel vectors.  	


//int loadSnpFile(char *snpFile, char ***chromNames, unsigned **chromSizes)
int loadSnpFile(char *snpFile, khash_t(hashChr_t) *hChr,  snplist_t *snplistp, indellist_t *indellistp)
{
  struct lineFile *lf = lineFileOpen(snpFile, TRUE);
  int k,j;

  unsigned int pos;

  int chrnum=kh_size(hChr)+1;  // I think I need to fix this... on the main program. It is a problem of not adding the empty chromosome?...

  char *row[10];
  int wordCount;
  char *chr_str;
  char *ref;
  char *alt;
  short int reflen;
  short int altlen;

  int snpCount=0;
  int indelCount=0;

  verbose(1,"# Reading SNP file: %s (1-based coordinates)\n",snpFile);

  //Initilalization of the snplist
  kv_init_size(snpvec_t,*snplistp,chrnum); 
  for(k=0;k<chrnum;k++)
    kv_init(kv_A(*snplistp,k));
  //Initilalization of the indellist
  kv_init_size(indelvec_t,*indellistp,chrnum); 
  for(k=0;k<chrnum;k++)
    kv_init(kv_A(*indellistp,k));
  
  // Read line by line 
  while ((wordCount = lineFileChop(lf, row)) != 0){
    //fscanf(bedF,"%s\t%d\t%d\t%c\%[^\n]\n",chr_str,&left,&right,&cStrand,buff);
    //  if(fscanf(fSnp,"%s\t%d\t%c\t%c%*[^\n]\n",cbuff,&snpaux.pos,&snpaux.ref,&snpaux.alt)!=4 || feof(fSnp))

    /* if(!((wordCount>=4)&&(wordCount<=8))){ */
    /*   for(k=0;k<wordCount;k++) */
    /* 	verbose(1,"%s",row[k]); */
    /*   verbose(1,"\n"); */
    /* }else{ */
    assert((wordCount>=4)&&(wordCount<=8));
    chr_str = row[0];
    pos = lineFileNeedNum(lf, row, 1)-1;

    //Determine if chr is valid 
    k=0;
    k = kh_get(hashChr_t, hChr , chr_str);
    assert(kh_exist(hChr,k));
    //assert(strcmp(kh_key(hChr,k),chr_str)==0);
    k = kh_val(hChr, k);

    // Determine if pos is valid ??? (left for later???)

    ref=row[2];    
    reflen=strlen(ref);
    //Loop if wordcount >4
    for(j=3;j<wordCount;j++){
      alt=row[j];
      altlen=strlen(alt);
      
      if((reflen>1)||(altlen>1)){ // INDEL
	indel_t aux;
	aux.pos=pos;
	aux.reflen=reflen;
	aux.altlen=altlen;
	aux.ref=cloneString(ref);
	aux.alt=cloneString(alt);
	kv_push(indel_t,kv_A(*indellistp,k),aux);
	indelCount++;
      }// INDEL
      else{ // SNP
	snp_t snpaux;
	snpaux.pos=pos;
	snpaux.ref=ref[0];
	snpaux.alt=alt[0];
	kv_push(snp_t,kv_A(*snplistp,k),snpaux);
	snpCount++;
      }// SNP
      //    }
    }// Alt alleles.... 

  }// while linefile

    verbose(1,"# There are %d snps and %d indels in the file\n",snpCount,indelCount);
  
    //  for(k=1;k<chrnum;k++)
    //   verbose(2,"#%d) %s %d Snps\n",k,chrNames[k],(int)kv_size(kv_A(snplist,k)));

    //close lf file...
    lineFileClose(&lf);
    return (snpCount+indelCount);
}



/* ******************************************************************************************** */

/*
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
  //#pragma omp parallel for private(nonATGCbase,j,score) reduction(+:count)
  for(i=Start;i<=End;i++){
    nonATGCbase=0;
    for(j=i;j<=(i+pwm->w-1);j++)
      if(ATGCbase(seq->dna[j])){
	nonATGCbase++;
	//shift i to skip??
      }
    if(nonATGCbase==0){
      score=compare_subseq_to_pssm(seq->dna+i,pwm); //Forward strand
      if(score >= thresh){
	fprintf(stdout,"%s\t%d\t%d\t+\t%f\n", seq->name,i,i+pwm->w-1,score);
	count++;
      }
      score=compare_subseq_to_pssm(seq->dna+i,rpwm);  //Reverse strand
      if(score >= thresh){
	fprintf(stdout,"%s\t%d\t%d\t-\t%f\n", seq->name,i,i+pwm->w-1,score);
	count++;
      }   	
    }  
  }

  return count;
}

*/

/* ******************************************************************************************** */


//    scankmersnp(seq,h,&replist,kh_val(hChr , k), kv_A(snplist,k));   
int scanPwmSeqSNP(struct dnaSeq *seq, struct pssm *pwm, snpvec_t *snpvecp){
  int i,j,n,k,Start,End,count;
  char nonATGCbase;
  double score;
  double MaxRefScore;
  double MaxAltScore;
  char MaxRefStrand;
  char MaxAltStrand;
  char oldval;
  int iMaxRef;
  int iMaxAlt;
  struct pssm revp,*rpwm;

  snp_t snpaux;

  rpwm=&revp;
  allocateMemoryToMatrix(rpwm);
  rpwm->w=pwm->w;
  reversePWM(rpwm,pwm);
  
  count=-1;

  for(j=0; j< kv_size(*snpvecp); j++){    
    snpaux=kv_A(*snpvecp,j);
    //assert(((snpaux.ref) & 0xDF) == ((seq->dna[snpaux.pos]) & 0xDF)); //Asserting Reference matches SNP
    assert(snpaux.pos < seq->size); //SNP out of limits

    Start = snpaux.pos - pwm->w + 1;
    if(Start<0)
      continue;
    End = snpaux.pos;// + pwm->w; //+seq->size;//100;//seq->size-pwm->w;
    if(End > (seq->size - 1))
      continue;
    if((End-Start+1) < pwm->w)
      continue; //Region too small to fit the kmer..

    MaxRefScore=-1E100;
    MaxAltScore=-1E100;
    iMaxRef=Start-1;
    iMaxAlt=Start-1;
    MaxRefStrand='.';
    MaxAltStrand='.';

    //
    oldval = seq->dna[snpaux.pos];
    seq->dna[snpaux.pos] = snpaux.ref; // This restores the to the previous state... 

    //#pragma omp parallel for private(nonATGCbase,j,score) reduction(+:count)
    for(i=Start,n=0;i<=End;i++,n++){
      nonATGCbase=0;
      for(k=i;k<=(i+pwm->w-1);k++)
	if(ATGCbase(seq->dna[k],maskRep)){
	  nonATGCbase++;
	  //shift i to skip??
	}
      if(nonATGCbase==0){
	score=compare_subseq_to_pssm(seq->dna+i,pwm,useSnpRobust)/log(base); //Forward strand
	if(score >= MaxRefScore){
	  MaxRefScore=score;
	  iMaxRef=i;
	  MaxRefStrand='+';
	}
	score=compare_subseq_to_pssm(seq->dna+i,rpwm,useSnpRobust)/log(base);  //Reverse strand
	if(score > MaxRefScore){
	  MaxRefScore=score;
	  iMaxRef=i;
	  MaxRefStrand='-';
	}
      }  
    }    
    
    seq->dna[snpaux.pos]=snpaux.alt; // Generates the alternate version...

    for(i=Start,n=0;i<=End;i++,n++){
      nonATGCbase=0;
      for(k=i;k<=(i+pwm->w-1);k++)
	if(ATGCbase(seq->dna[k],maskRep)){
	  nonATGCbase++;
	  //shift i to skip??
	}
      if(nonATGCbase==0){
	score=compare_subseq_to_pssm(seq->dna+i,pwm,useSnpRobust)/log(base); //Forward strand
	if(score >= MaxAltScore){
	  MaxAltScore=score;
	  iMaxAlt=i;
	  MaxAltStrand='+';
	}
	score=compare_subseq_to_pssm(seq->dna+i,rpwm,useSnpRobust)/log(base);  //Reverse strand
	if(score > MaxAltScore){
	  MaxAltScore=score;
	  iMaxAlt=i;
	  MaxAltStrand='-';
	}
      }  
    }    
    
    //    if( ((MaxRefScore>=thresh) || (MaxAltScore>=thresh)) &&  ((MaxRefScore>=0) || (MaxAltScore>=0)) ){ 
    if( ((MaxRefScore>=thresh) || (MaxAltScore>=thresh)) ){ 
      n=(iMaxRef - Start);
      k=(iMaxAlt - Start);
      if(MaxRefStrand=='+')
	n = pwm->w - 1 - n;
      if(MaxAltStrand=='+')
	k = pwm->w - 1 - k;      
      fprintf(stdout,
	      "%s\t%d\t%d\tR"
	      "\t%f\t%c\t%d"
	      "\t%f\t%c\t%d"
	      "\t%d\t%c\t%c\n"
	      , seq->name,iMaxRef,iMaxRef+pwm->w-1
	      ,MaxRefScore,MaxRefStrand,n
	      ,MaxAltScore,MaxAltStrand,k
	      ,snpaux.pos,snpaux.ref,snpaux.alt);
      count++;
    }
    
    seq->dna[snpaux.pos]=snpaux.ref; // This restores the to the previous state... 
    seq->dna[snpaux.pos]=oldval; // This restores the to the previous state... 

    //    lastpos=snpaux.pos;
    
  }
    return 1;
}

int scanPwmSeqInDel(struct dnaSeq *seq, struct pssm *pwm, indelvec_t *indelvecp){
  int i,j,n,k,Start,End,count;
  char nonATGCbase;
  double score;
  double MaxRefScore;
  double MaxAltScore;
  char MaxRefStrand;
  char MaxAltStrand;
  //  char oldval;
  int iMaxRef;
  int iMaxAlt;
  struct pssm revp,*rpwm;

  char buff[1000];

  //snp_t snpaux;
  indel_t indel;

  rpwm=&revp;
  allocateMemoryToMatrix(rpwm);
  rpwm->w=pwm->w;
  reversePWM(rpwm,pwm);
  
  count=-1;

  for(j=0; j< kv_size(*indelvecp); j++){    
    indel=kv_A(*indelvecp,j);

    //assert(((snpaux.ref) & 0xDF) == ((seq->dna[snpaux.pos]) & 0xDF)); //Asserting Reference matches SNP
    //assert(indelaux.pos < seq->size); //SNP out of limits
    // Assert Indel matches de reference....
    if(strncasecmp(&seq->dna[indel.pos],indel.ref,indel.reflen)!=0){
      verbose(1,"## %s:%d,%s\n",seq->name,indel.pos,(indel.ref));
      errAbort("Indel not matching the reference\n");
      //or use continue and issue a warning
    }

    Start = indel.pos - pwm->w + 1;
    if(Start<0)
      continue;
    End = indel.pos + indel.reflen -1;// + pwm->w; //+seq->size;//100;//seq->size-pwm->w;
    if(End > (seq->size - 1))
      continue;
    if((End-Start+1) < pwm->w)
      continue; //Region too small to fit the kmer..

    MaxRefScore=-1E100;
    MaxAltScore=-1E100;
    iMaxRef=Start-1;
    iMaxAlt=Start-1;
    MaxRefStrand='.';
    MaxAltStrand='.';


    // AS before scan the reference sequence!
    //#pragma omp parallel for private(nonATGCbase,j,score) reduction(+:count)
    for(i=Start;i<=End;i++){
      nonATGCbase=0;
      for(k=i;k<=(i+pwm->w-1);k++)
	if(ATGCbase(seq->dna[k],maskRep)){
	  nonATGCbase++;
	  //shift i to skip??  MAYBE I WANT TO ISSUE A NA IN THIS CASE
	}
      if(nonATGCbase==0){
	score=compare_subseq_to_pssm(seq->dna+i,pwm,useSnpRobust)/log(base); //Forward strand
	if(score >= MaxRefScore){
	  MaxRefScore=score;
	  iMaxRef=i;
	  MaxRefStrand='+';
	}
	score=compare_subseq_to_pssm(seq->dna+i,rpwm,useSnpRobust)/log(base);  //Reverse strand
	if(score > MaxRefScore){
	  MaxRefScore=score;
	  iMaxRef=i;
	  MaxRefStrand='-';
	}
      }  
    }    

    // Extract DNA sequence from reference...
    i=0;
    for(k=Start;k<indel.pos;k++) // First section...
      buff[i++]=seq->dna[k];
    //buff[i++]='_';
    for(k=0;k<indel.altlen;k++)
      buff[i++]=indel.alt[k];
    //buff[i++]='_';
    for(k=indel.pos+indel.reflen;k<=End;k++)
      buff[i++]=seq->dna[k];
    k=i;


    //#pragma omp parallel for private(nonATGCbase,j,score) reduction(+:count)
    for(i=0;i<k;i++){
      nonATGCbase=0;     
      for(n=i;n<=(i+pwm->w-1);n++)
	if(ATGCbase(buff[n],maskRep)){
	  nonATGCbase++;
	  //shift i to skip??
	}
      if(nonATGCbase==0){
	score=compare_subseq_to_pssm(buff+i,pwm,useSnpRobust)/log(base); //Forward strand
	if(score >= MaxAltScore){
	  MaxAltScore=score;
	  iMaxAlt=Start+i;
	  MaxAltStrand='+';
	}
	score=compare_subseq_to_pssm(buff+i,rpwm,useSnpRobust)/log(base);  //Reverse strand
	if(score > MaxAltScore){
	  MaxAltScore=score;
	  iMaxAlt=Start+i;
	  MaxAltStrand='-';
	}
      }  
    }    
    
    //    if( ((MaxRefScore>=thresh) || (MaxAltScore>=thresh)) &&  ((MaxRefScore>=0) || (MaxAltScore>=0)) ){ 
    if( ((MaxRefScore>=thresh) || (MaxAltScore>=thresh)) ){ 
      n=(iMaxRef - Start);
      k=(iMaxAlt - Start);
      if(MaxRefStrand=='+')
	n = pwm->w - 1 - n;
      if(MaxAltStrand=='+')
	k = pwm->w - 1 - k;      
      fprintf(stdout,
	      "%s\t%d\t%d\tR"
	      "\t%f\t%c\t%d"
	      "\t%f\t%c\t%d"
	      "\t%d\t%s\t%s\n"
	      , seq->name,iMaxRef,iMaxRef+pwm->w-1
	      ,MaxRefScore,MaxRefStrand,n
	      ,MaxAltScore,MaxAltStrand,k
	      ,indel.pos,indel.ref,indel.alt);
      count++;
    }
    
    //    lastpos=snpaux.pos;
    
  }
    return 1;
}




/*
  //DNA buff[kmerSize+1];
  //unsigned long int ForKmer2=0,RevKmer2=0;
  unsigned long int ForKmer=0,RevKmer=0;
  
  unsigned int mask=0xFFFFFFFF;  //0x000FFFFF
  char Base;
  unsigned long int b;

  int lastpos=-1000; // to check snps ordered... 
*/


/* ******************************************************************************************** */


/* ******************************************************************************** */
/* ******************************************************************************** */
/* ******************************************************************************** */


/*
void scanPwm(char *fileMotif, char *fileSeq)
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
*/


void scanPwmVar(char *fileMotif, char *fileSeq, char *fileSNPs)
/* scanPwmVar - Scan binding sites hit by SNPs. */
{
  struct dnaLoad *dl = dnaLoadOpen(fileSeq); // change to 2bit file alone?. 


  struct dnaSeq *seq; 
  int j;
  khiter_t k;
  //  size_t vsize;
  //  size_t cumsize=0;


  int hret;
  char *chrNames[256];
  //clock_t t;

  unsigned char chrnum=1;

  //khash_t(hashPos_t) *h = kh_init(hashPos_t);

  khash_t(hashChr_t) *hChr = kh_init(hashChr_t);  

  //repvec_t *repvecp;  
  replist_t replist;
  kv_init(replist);  
  
  //snpvec_t *snpvecp;  
  indellist_t indellist;
  snplist_t snplist;
  // snp_t snpaux;

  struct pssm pwm;
  initialise_pssm(&pwm,fileMotif,addPseudoCounts,useJaspar);

  while ((seq = dnaLoadNext(dl)) != NULL){
    k = kh_put(hashChr_t, hChr , cloneString(seq->name), &hret);
    if(hret){
      //assert chrnum<255
      kh_val(hChr, k) = chrnum;
      chrNames[chrnum]=cloneString(seq->name);
      chrnum++;
    }
    else{
      errAbort("Duplicated chromosome \n");
    }
    verbose(1, "#\tprocessing: %s %d %s key:%s %d\n", seq->name, (int)kh_val(hChr, k), chrNames[kh_val(hChr, k)],kh_key(hChr, k), hret);
    // //     scankmer(seq,h,&replist,kh_val(hChr , k));   
    //      scanPwmOneSeq(seq,&pwm);   
    
    dnaSeqFree(&seq);
    //if(chrnum>2)break;
  } 
  //  verbose(1,"Finished hashing reference genome\n");
  //  verbose(1,"khash: n_buckets=%d, size=%d, n_occupied=%d, upper_bound=%d \n",kh_n_buckets(h),kh_size(h),((h)->n_occupied),((h)->upper_bound));  
  dnaLoadClose(&dl);

  /* --------------------------------------------------------------------- */

  loadSnpFile(fileSNPs, hChr, &snplist, &indellist);

  verbose(1,"Now SNPs\n");
  for(k=1;k<chrnum;k++)
      verbose(2,"#  %d) %s %d Snps + %d Indels\n",k,chrNames[k],(int)kv_size(kv_A(snplist,k)),(int)kv_size(kv_A(indellist,k)));
  
    /*   
	 k=49;
	 for(j=0;j<10;j++){
      snp_t aux;
      aux=kv_A(kv_A(snplist,k),j);
      verbose(2,"%d %c %c\n",aux.pos,aux.ref,aux.alt);
      }*/
      
    /* --------------------------------------------------------------------- */
  
  verbose(1,"Scanning pwm using the SNP file \n");  
  dl = dnaLoadOpen(fileSeq);
  while ((seq = dnaLoadNext(dl)) != NULL){
    k = kh_get(hashChr_t, hChr , seq->name);
    j = kh_val(hChr, k);
    verbose(1,"#%d) Re-scanning %d SNPs in %s (1 single alt.)\n",
	    j, (int)kv_size(kv_A(snplist,j)), chrNames[j]);

    //	scankmersnp(seq, h, &replist, j, &kv_A(snplist,j));   
    scanPwmSeqSNP(seq,&pwm, &kv_A(snplist,j));

    verbose(1,"#%d) Re-scanning %d INDELs in %s (1 single alt.)\n",
	    j, (int)kv_size(kv_A(indellist,j)), chrNames[j]);

    scanPwmSeqInDel(seq,&pwm, &kv_A(indellist,j));

    /*
    verbose(1,"#%d) Re-hashing %d INDELs in %s (1 single alt.)\n",
	    j, (int)kv_size(kv_A(indellist,j)), chrNames[j]);
    // scankmerindel(seq, h, &replist, j, &kv_A(indellist,j));
    */
    /*
      if(hashAllSnpOpt){
      verbose(1,"#%d) Re-hashing %d SNPs in %s (complete alt.), idx=%d\n", 
      j, (int)kv_size(kv_A(snplist,j)), chrNames[j], idx);
      scankmersnp2(seq, h, &replist, j, &kv_A(snplist,j));   
      }
    */
    dnaSeqFree(&seq);
    //if(chrnum>5)break;
  }
  verbose(1,"Finished scanning alternative genomes\n");
  //    verbose(2,"khash: n_buckets=%d, size=%d, n_occupied=%d, upper_bound=%d \n",kh_n_buckets(h),kh_size(h),((h)->n_occupied),((h)->upper_bound));

  dnaLoadClose(&dl);

  
  
  /* --------------------------------------------------------------------- */  
  
}




int main(int argc, char *argv[])
/* Process command line. */
{
optionInit(&argc, argv, options);
if (argc != 4)
    usage();

 dnaUtilOpen();

 useJaspar = optionExists("j");
 useSnpRobust = optionExists("snp");
 maskRep = optionExists("mask");

 thresh= optionDouble("t",thresh);
 base= optionDouble("base",base);
 addPseudoCounts= optionDouble("p",addPseudoCounts);
 ompNumThreads= optionInt("omp",ompNumThreads);


 // omp_set_num_threads(ompNumThreads); //Adjust this to env variable.... or option


// wigOutput = optionExists("wigOutput");

 scanPwmVar(argv[1],argv[2],argv[3]);
return 0;
}

