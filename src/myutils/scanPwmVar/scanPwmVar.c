/* scanPwmVar - Scan binding sites hit by SNPs. */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"

#include "dnautil.h"
#include "dnaseq.h"
#include "dnaLoad.h"

#include <omp.h>


static char const rcsid[] = "$Id: newProg.c,v 1.30 2010/03/24 21:18:33 hiram Exp $";

void usage()
/* Explain usage and exit. */
{
errAbort(
  "scanPwmVar - Scan binding sites hit by SNPs\n"
  "usage:\n"
  "   scanPwmVar XXX\n"
  "usage:\n"
  "   scanPwm motif.pwm sequence\n"
  "where:\n"
  "   sequence is a .fa , .nib or .2bit file or a file which is a list of sequence files.\n"
  "options:\n"
  "   -j - use Jaspar *.pfm format (the default is Dan's TRANSFAC matrix)\n"
  "   -snp - discount most negative base in the pwm aligment if unkown snp hits the score \n"
  "   -t=6.91 - use threshold cut-off \n"
  "   -base=2.71828183  - use log in that base (default natural log) \n"
  "   -p=0.0 - add p pseudo-counts to frequencies (0.01 default to aboid -INF score) \n"
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


double addPseudoCounts = 0.01;
double thresh = 6.91 ;
double base = 2.71828183 ;
boolean useJaspar = TRUE;
boolean useSnpRobust =FALSE;
boolean maskRep = FALSE;
int ompNumThreads= 1;

#define maxPSSM 200
#define maxStr 500

struct pssm
{
  void *next;
  char *name;
  int w;
  double **matrix; /* */ 
};


void readInJasparPwm(struct pssm *pm, char *fileName);
void readInTransfacPwm(struct pssm *pm, char *fileName);
void allocateMemoryToMatrix (struct pssm *pwm);
void convertPSSMToLogs (struct pssm *pwm);
void reversePWM(struct pssm *op,const struct pssm *p);
void initialise_pssm(struct pssm *pm,char *fileName);
void printMatrix (struct pssm *pm);  
//void read_seq_fasta(struct parameters *p, struct sequence *s);
//void print_seq_fasta(FILE *f, struct sequence s);
int ATGCbase(char b);
double compare_subseq_to_pssm (char *seq, struct pssm *p);

//int get_pwm_scores(struct seq *seq, struct pssm *p);

void initialise_pssm(struct pssm *pwm,char *fileMotif) 
{
  int n,b,NumZeros=0;
  double FirstRowCounts=0.0,RowCounts=0.0;

  allocateMemoryToMatrix(pwm);
  if (useJaspar)
    readInJasparPwm(pwm,fileMotif);
  else
    readInTransfacPwm(pwm,fileMotif); 

  //Check if PWM is consistent and asses if probabilities or counts are provided. 
  for(n=0;n<pwm->w;n++){
    RowCounts=0.0;
    for(b=0;b<4;b++)
      RowCounts=RowCounts+pwm->matrix[n][b];
    for(b=0;b<4;b++)
      if(pwm->matrix[n][b]==0.0)
	NumZeros++;
    fprintf(stderr,"# Row %d: %g Counts\n",n,RowCounts);
    if(n==0)
      FirstRowCounts=RowCounts;
    //  else
      //assert(fabs(FirstRowCounts-RowCounts)<0.025);
      // if(fabs(FirstRowCounts-RowCounts)<0.025)
      //fprintf("#Row counts are not consistent across\n");
  }
  //fprintf(stderr, "Sum of first row %d\n",FirstRowCounts); 
  //Add pseudo_counts -- only possible if counts are provided...
  //Matrix is converted to probabilities to ...
  if(NumZeros>0)
    if(FirstRowCounts<1.1)
      addPseudoCounts=FirstRowCounts/100; 
  if(NumZeros==0)
    addPseudoCounts=0.0;  //Not necessary to add pseudo counts...

  // printMatrix(pwm);  
  fprintf(stderr,"#Adding pseudocounts %f\n",addPseudoCounts);
  for(n=0;n<pwm->w;n++){
    RowCounts=0.0;
    for(b=0;b<4;b++)
      pwm->matrix[n][b]=(pwm->matrix[n][b]+addPseudoCounts);
    for(b=0;b<4;b++)
      RowCounts=RowCounts+pwm->matrix[n][b];
    for(b=0;b<4;b++)
      pwm->matrix[n][b]=(pwm->matrix[n][b])/RowCounts;
  }
  printMatrix(pwm);  
  convertPSSMToLogs(pwm);
  printMatrix(pwm); 
}

void readInJasparPwm(struct pssm *pm, char *pwmFilePath)
{
  int row=0,col=0,num_rows=0;
  double c;
  FILE *inf;
  //  fprintf(stdout,"readInJasparPWM:\n");

  inf = mustOpen(pwmFilePath,"r");

  while(!feof(inf))
    {
    fscanf(inf,"%lf",&c);
    //   fprintf(stdout,"%d,",c);
    pm->matrix[row][col]=(double)c;    
    // fprintf(stdout,"%d,%f;",c,pm->matrix[row][col]);
    row++;
    if((c=getc(inf))=='\n')
      {      
        if(col==0)
	  num_rows=row;
	else
	  assert(num_rows==row); 
	
	row=0;
	col++;
	// fprintf(stdout,"--%d \n",col);
	if(col>3)
	  break;
      }
    else
      ungetc(c,inf);
    }

  fclose(inf);
  pm->w = num_rows;
  //p->pwmLength = num_rows;
  printMatrix(pm);  
}

void readInTransfacPwm(struct pssm *pm, char *pwmFilePath)
{
  int s,f=0,row=0,commentLine=0,nComLines=0,startPSSM=0;
  char str_buf,firstBit[maxStr];
  FILE *inf;

  inf = mustOpen(pwmFilePath,"r");

  while(!feof(inf)) {
    if(row >= maxPSSM) {
      fprintf(stderr,"PSSM too long (max is %d bp)\n",maxPSSM);
      exit(1);
    }
    if(startPSSM == 0) {
      s = fscanf(inf,"%c",&str_buf);
/*       printf("%c",str_buf); */
      if(s == EOF) break;
      if(commentLine) {
	if(str_buf == '\n') {
	  commentLine = 0;
	}
	continue;
      }
      if(str_buf == '#' || str_buf == 'A') { /* comment or header line  */
	commentLine++;
	nComLines++;
	continue;
      }
      startPSSM++;
    } 
    if(startPSSM > 0) {
      if(startPSSM==1) { /* just found the PSSM */
	firstBit[f] = str_buf; /* horrible hack to make sure we get the first digit on the first line */
	while((s = fscanf(inf,"%c",&str_buf))) {
	  if(str_buf == '\t' || str_buf == ' ') {
	    break;
	  }
	  f++;
	  firstBit[f] = str_buf;
	}
	firstBit[f+1] = '\0';
	pm->matrix[row][0] = atof(firstBit);
	s = fscanf(inf,"%lf %lf %lf",&pm->matrix[row][1],&pm->matrix[row][2],&pm->matrix[row][3]);
	startPSSM = 100;
      } else {
	s = fscanf(inf,"%lf %lf %lf %lf",&pm->matrix[row][0],&pm->matrix[row][1],&pm->matrix[row][2],&pm->matrix[row][3]);
	if(s == EOF) break;
      }
      row++;
    }
  }
  pm->w = row;
  //  p->pwmLength = row;
  fclose(inf);
}

void allocateMemoryToMatrix (struct pssm *p) {
  int i;
  p->matrix = malloc(maxPSSM*sizeof(double *));
 if(p->matrix == NULL) {
    fprintf(stderr,"Insufficient memory for p->matrix\n");
    exit(1);
  }
  for(i = 0;i < maxPSSM;i++) {
    p->matrix[i] = malloc(4*sizeof(double));
    if(p->matrix[i] == NULL) {
      fprintf(stderr,"Insufficient memory for p->matrix element %d\n",i);
      exit(1);
    }
    p->matrix[i][0] = 0.0;
  }
}

void convertPSSMToLogs (struct pssm *p) {
  int i,j;
  for(i=0;i<p->w;i++) {
    for(j=0;j<4;j++) {
      if(p->matrix[i][j]==0.0) { 
	fprintf(stderr,"Need to add pseudocounts to zero columns\n");
	exit(1);
      }
      p->matrix[i][j] = log(p->matrix[i][j]);
    }
  }
}

void reversePWM(struct pssm *op,const struct pssm *p){
  int i,I;
  double rowtmp[4];

  I=p->w;
  for(i=0;i<(I/2);i++) {
    rowtmp[0]=p->matrix[i][3];
    rowtmp[1]=p->matrix[i][2];
    rowtmp[2]=p->matrix[i][1];
    rowtmp[3]=p->matrix[i][0];
    op->matrix[i][0]=p->matrix[I-1-i][3];
    op->matrix[i][1]=p->matrix[I-1-i][2];
    op->matrix[i][2]=p->matrix[I-1-i][1];
    op->matrix[i][3]=p->matrix[I-1-i][0];
    op->matrix[I-1-i][0]=rowtmp[0];
    op->matrix[I-1-i][1]=rowtmp[1];
    op->matrix[I-1-i][2]=rowtmp[2];
    op->matrix[I-1-i][3]=rowtmp[3];
  }
  if((I%2)==1){ //if number of rows is odd we need to do this. 
    rowtmp[0]=p->matrix[i][3];
    rowtmp[1]=p->matrix[i][2];
    rowtmp[2]=p->matrix[i][1];
    rowtmp[3]=p->matrix[i][0];
    op->matrix[i][0]=rowtmp[0];
    op->matrix[i][1]=rowtmp[1];
    op->matrix[i][2]=rowtmp[2];
    op->matrix[i][3]=rowtmp[3];
  }    
}

void printMatrix (struct pssm *pm) {  
  int i;
  fprintf(stderr,"# N\tA\tC\tG\tT\n");
  for(i = 0;i < pm->w;i++) {
    fprintf(stderr,"# %d\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n",i,pm->matrix[i][0],pm->matrix[i][1],pm->matrix[i][2],pm->matrix[i][3]);
  }
}

/// Returns 1 if non ATGC letter...
int ATGCbase(char b) {
  if(!maskRep){
    if(b=='A' || b=='a') {
      return 0;
    } else if(b=='C' || b=='c' ) {
      return 0;
    } else if(b=='G' || b=='g' ) {
      return 0;
    } else if(b=='T' || b=='t' ) {
      return 0;
    }
    return 1;
  }
  else{
    if(b=='A') {
      return 0;
    } else if(b=='C') {
      return 0;
    } else if(b=='G') {
      return 0;
    } else if(b=='T') {
      return 0;
    }
    return 1; 
  }
}

/// Compares subseq to PWM
double compare_subseq_to_pssm (char *seq, struct pssm *p) {
  int i;
  char b;
  double score=0,aux;
  double vmin=0;
  for(i=0;i<p->w;i++) {
    aux=0;
    b=seq[i];
    if(b == 'A'  || b=='a') {
      aux = p->matrix[i][0];
    } else if(b == 'C'  || b=='c') {
      aux = p->matrix[i][1];
    } else if(b == 'G'  || b=='g') {
      aux = p->matrix[i][2];
    } else if(b == 'T'  || b=='t') {
      aux = p->matrix[i][3];
    } else {
      aux = -100;
      /*       fprintf(stderr,"Skipping ambiguous character %c\n",seq[i]); */
    }
    if(aux<vmin)
      vmin=aux;
    score += aux;
  }
  /*   printf("%lf %lf\n",score,p.w*log(0.25)); */
  

  if(useSnpRobust)
    score = score - vmin - p->w*log(0.25);
  else
    score = score - p->w*log(0.25);

  return score/log(base);  // 
}



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
	if(ATGCbase(seq->dna[k])){
	  nonATGCbase++;
	  //shift i to skip??
	}
      if(nonATGCbase==0){
	score=compare_subseq_to_pssm(seq->dna+i,pwm); //Forward strand
	if(score >= MaxRefScore){
	  MaxRefScore=score;
	  iMaxRef=i;
	  MaxRefStrand='+';
	}
	score=compare_subseq_to_pssm(seq->dna+i,rpwm);  //Reverse strand
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
	if(ATGCbase(seq->dna[k])){
	  nonATGCbase++;
	  //shift i to skip??
	}
      if(nonATGCbase==0){
	score=compare_subseq_to_pssm(seq->dna+i,pwm); //Forward strand
	if(score >= MaxAltScore){
	  MaxAltScore=score;
	  iMaxAlt=i;
	  MaxAltStrand='+';
	}
	score=compare_subseq_to_pssm(seq->dna+i,rpwm);  //Reverse strand
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



void scanPwm(char *fileMotif, char *fileSeq)
/* scanPwm - scans a sequence for motif matches to a PWM. */
{
  struct dnaLoad *dl = dnaLoadOpen(fileSeq);
  struct dnaSeq *seq; 

  struct pssm pwm;

  initialise_pssm(&pwm,fileMotif);
  
  //motifList = fixMotifs(dnaMotifLoadAll(motifFile));

  while ((seq = dnaLoadNext(dl)) != NULL)
    {
      verbose(2, "#\tprocessing: %s\n", seq->name);
      scanPwmOneSeq(seq,&pwm);   
    }
}



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
  initialise_pssm(&pwm,fileMotif);

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


 // omp_set_num_threads(ompNumThreads); //Adjust this to env variable.... or option



// wigOutput = optionExists("wigOutput");

 scanPwmVar(argv[1],argv[2],argv[3]);
return 0;
}

