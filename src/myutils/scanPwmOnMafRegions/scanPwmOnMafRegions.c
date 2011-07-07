/* scanPwmOnMafRegions - scan a Pwm in the multiple aligments of a maf on the regions given by a bed file. */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"

//#include "jksql.h"

#include "sqlNum.h"
#include "maf.h"
#include "bed.h"

#include "dnautil.h"
#include "dnaseq.h"
#include "dnaLoad.h"

#include <omp.h>


static char const rcsid[] = "$Id: newProg.c,v 1.28 2009/07/02 17:14:39 angie Exp $";

void usage()
/* Explain usage and exit. */
{
errAbort(
  "scanPwmOnMafRegions - scan a Pwm in the multiple aligments of a maf on the regions given by a bed file\n"
  "usage:\n"
  "   scanPwmOnMafRegions motif.pwm regions.bed species.txt out.tab multiz.maf(s)\n"
  "   species.txt is single column file with the list of species to extract from the maf files\n"
  "example:\n"
  "   scanPwmOnMafRegions -base=2 -p=0.05 -window=30 M00256.dat M00256.bed.gz species.txt M00256.maf.scores chr*.maf\n"
  "options:\n"
  "   -j - use Jaspar *.pfm format (the default is Dan's TRANSFAC matrix)\n"
  "   -base=2.71828183  - use log in that base (default natural log) \n"
  "   -window=0 - extending the region on the aligment to look for the motif\n"
  "   -p=0.0 - add p pseudo-counts to frequencies (0.01 default to aboid -INF score) \n"
  //  "   -omp=1 - multi-core implementation (default 1 == one thread)" 
  //  "   -mask - remove repeat masked bases \n" // look for another thing.
  // Add an option for output format.
  "   -xxx=XXX\n"
  );
}

static struct optionSpec options[] = {
   {"j", OPTION_BOOLEAN},
   {"mask", OPTION_BOOLEAN},
   {"p", OPTION_DOUBLE},
   {"base", OPTION_DOUBLE},
   {"window", OPTION_INT},
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
int window=0;

/***************************************************************/
/* TODO: Put this in a separate library identical in scanPwm.c */
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
void freeMatrix (struct pssm *p);
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

void freeMatrix (struct pssm *p) {
  int i;
  for(i = 0;i < maxPSSM;i++)
    free(p->matrix[i]);
  free(p->matrix);
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

double scanPwmOneSeqTopHit(struct dnaSeq *seq,struct pssm *pwm, int *iMaxRet, char *sMaxRet,double *maxScoreRet)
{
  int i,j,Start,End;
  char nonATGCbase;
  int iMax;
  char sMax;
  double score,maxScore;
  struct pssm revp,*rpwm;
  
  rpwm=&revp;
  allocateMemoryToMatrix(rpwm);
  rpwm->w=pwm->w;
  reversePWM(rpwm,pwm);

  maxScore=-1000;
  sMax='.';
  iMax=-1;
  *iMaxRet=iMax;
  *sMaxRet=sMax;
  *maxScoreRet=maxScore;

  
  Start=0;
  End=seq->size-pwm->w;
  
  if(seq->size < pwm->w)
    return -2000; //maxScore; //Region too small to fit the motif..

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
      if( score > maxScore ){
	maxScore=score;
	sMax='+';
	iMax=i;
      }
      score=compare_subseq_to_pssm(seq->dna+i,rpwm);  //Reverse strand
      if(score > maxScore){
	maxScore=score;
	sMax='-';
	iMax=i;
      }
    }  
  }

  *iMaxRet=iMax;
  *sMaxRet=sMax;
  *maxScoreRet=maxScore;

  //Need to free rpwm
  freeMatrix(rpwm);

  return maxScore;
}




/* ******************************************************************************** */
/* ******************************************************************************** */
/* ******************************************************************************** */
/* ******************************************************************************** */
/* MAF stuff  */

//boolean outDir = FALSE;
boolean keepInitialGaps = FALSE; // TRUE;
//char *dir = NULL;
//char *scoring = NULL;

struct hash *loadRegions(char *file)
/* load regions into a hash of lists by chrom */
{
  struct bed *bed = NULL, *bedList = NULL, *nextBed = NULL, *temp = NULL;
  struct hash *regionHash = newHash(6);
  struct bed *regions;

  regions = bedLoadNAll(file,3); //3 or 4 // number of fields in the bed
  /* order by chrom, start */
  slSort(&regions, bedCmp);
  verbose(2, "found %d regions\n", slCount(regions));
  bedList = regions;
  for (bed = regions; bed != NULL; bed = nextBed)
    {
      //Window extension! we do it in another point.
      // bed->chromStart = bed->chromStart - window;
      //bed->chromEnd = bed->chromEnd + window + 1; //

      verbose(3, "region %s:%d-%d\n", bed->chrom, bed->chromStart+1, bed->chromEnd);
      nextBed = bed->next;
      if ((bed->next == NULL) || (differentString(bed->chrom,bed->next->chrom)))
	{
	  temp = bed->next;
	  bed->next = NULL;
	  hashAdd(regionHash, bed->chrom, bedList);
	  verbose(2, "just added %d regions on %s\n", slCount(bedList), bed->chrom);
	  bedList = temp;
	}
    }
  return regionHash;
}

char *chromFromSrc(char *src)
/* get chrom name from <db>.<chrom> */
{
  char *p;
  if ((p = strchr(src, '.')) == NULL)
    errAbort("Can't find chrom in MAF component src: %s\n", src);
  return ++p;
}


char *addSuffixFreeHead(char *head, char *suffix)
/* Return a needMem'd string containing "headsuffix" and free head. Should be free'd
 when finished. */
{
  char *ret = NULL;
  int size = strlen(head) + strlen(suffix) +1;
  ret = needMem(sizeof(char)*size);
  snprintf(ret, size, "%s%s", head, suffix);
  freeMem(head);
  //fprintf(stdout,"%s,Buff: %s\n",shortname,ret);
  return ret;
}


double scanPwmOneSeqTopHitMyWrapper(char *buff, struct pssm *pwm)
{
  struct dnaSeq *myDna;
  double pwmScore=-1000;
  //fprintf(stdout,"%s ",buff);
  //dnaMixedCaseFilter(buff,dnatmp);
  myDna=newDnaSeq(buff,dnaFilteredSize(buff),"hg18");
  //pwmScore=scanPwmOneSeqTopHit(myDna,pwm);
  //fprintf(f,"%f\n",pwmScore);
  freeDnaSeq(&myDna);
  return pwmScore;
}


void extractMafs(char *file, FILE *f, struct hash *regionHash, struct hash *speciesHash, struct pssm *pwm)
/* extract MAFs in a file from regions specified in hash */
{
  char *chrom = NULL;
  struct bed *bed = NULL;
  struct mafFile *mf = mafOpen(file);
  struct mafAli *maf = NULL;
  struct mafComp *mc;
  struct mafComp *comp;

  struct dnaSeq *myDna;
  //  char path[256];
  //struct dnaSeq *myDna;
  // DNA *buff;
  // DNA *oldbuff;

  int mafNumber=0,i;
  int bedOverlap=0;
  int bedStart=0;
  int bedEnd=0;
  int numSpecies=speciesHash->elCount;
 
  DNA *dnabuff[numSpecies];  
  double pwmScore[numSpecies];
  char Rstatus[numSpecies];
  char Lstatus[numSpecies];

  int iMax;
  char sMax;
  double pwmScore2;


  for (i=0;i<numSpecies;i++)
    pwmScore[i]=-1000;

  verbose(1, "extracting from %s\n", file);
  maf = mafNext(mf);
  while (maf)
    {
      mc = maf->components;
      if (!chrom || differentString(chrom, chromFromSrc(mc->src)))
        chrom = cloneString(chromFromSrc(mc->src));         /* new chrom */
      bed = (struct bed *)hashFindVal(regionHash, chrom);
      if (!bed)
        {
	  /* no regions on this chrom -- skip to next chrom in maf */
	  do
            mafAliFree(&maf);
	  while (((maf = mafNext(mf)) != NULL) && sameString(chromFromSrc(maf->components->src), chrom));
	  continue;  // start over with this maf
        }

      /* Check if this region is overlapping with previous bed */
      if(bedOverlap)
	bedStart=bedEnd+1;
      else
	bedStart=bed->chromStart-window;
      
      /* Check if next region is overlapping with current bed */
      if( (bed->next) && sameString(bed->chrom,bed->next->chrom) && 
	  ((bed->chromEnd+1+window) > (bed->next->chromStart-window))){
	bedEnd=((bed->chromEnd+1 + bed->next->chromStart)/2); //at the midpoint
	bedOverlap=1;
      }
      else{
	bedEnd=bed->chromEnd+window+1;
	bedOverlap=0;
      }	
      
      verbose(2, "region: %s:%d-%d\n", 
	      bed->chrom, bedStart+1, bedEnd);
      /* skip mafs before region, stopping if chrom changes */
      while (maf && (mc = maf->components) && sameString(chrom, chromFromSrc(mc->src)) &&
	     (mc->start + mc->size) <= bedStart)
        {
	  mafAliFree(&maf);
	  maf = mafNext(mf);
        }
      

      mafNumber=0;
      //      buff=cloneString("");

      //if(bedOverlap==0){
	for (i=0;i<numSpecies;i++)
	  dnabuff[i]=cloneString(""); //Initialize dnabuff
	for (i=0;i<numSpecies;i++)
	  Rstatus[i]='-';
	for (i=0;i<numSpecies;i++)
	  Lstatus[i]='-';
	//}

      //for (i=0;i<numSpecies;i++)
      //	pwmScore[i]=-1000;

      /* extract all mafs and pieces of mafs in region */
      /* RPR, I'm working with small regions, so the region probably falls in a single maf or two if boundary */
      /* I think i should concatenate sequences first and then scan*/
      while (maf && (mc = maf->components) && sameString(chrom, chromFromSrc(mc->src)) &&
	     (bedStart < mc->start + mc->size && bedEnd > mc->start))
        {
	  int mafStart = mc->start;
	  int mafEnd = mc->start + mc->size;
	  struct mafAli *full = maf;
	  if (mafStart < bedStart || mafEnd > bedEnd)
            {
	      full = maf;
	      maf = mafSubsetE(full, mc->src, bedStart, bedEnd, keepInitialGaps);
	      mc = maf->components;
            }
	  verbose(2, "   %s:%d-%d\n", chrom, mc->start+1, mc->start + mc->size);
	  /* Code added by DG and Roger */ 
	  mafNumber++;
	  // f=stdout; //
	  for (comp = maf->components; comp != NULL; comp = comp->next) {
	    if ((comp->size == 0) && (comp->leftStatus)) {
	      /* 	    fprintf(f, "%s Empty\n",comp->src); */
	      /* Try to find why it is empty */
	      //fprintf(f, " Empty %s LeftStatus%c RightStatus%c \n",comp->src,comp->leftStatus,comp->rightStatus);
	      //fprintf(f,"%s:%d-%d\t",bed->chrom, bedStart+1, bedEnd);
	      //fprintf(f,"%s\tNA\n",comp->src);
	      char *shortname=cloneString(comp->src);
	      chopSuffix(shortname);
	      i=hashIntValDefault(speciesHash,shortname,-1);
	      if(i>-1){
		Rstatus[i]=comp->rightStatus;
		Lstatus[i]=comp->leftStatus;
	      }
	      freeMem(shortname);
	    } else {
	      DNA *dnatmp=cloneString(comp->text);
	      char *shortname=cloneString(comp->src);
	      //stripChar(dnatmp,'-');
	      //stripChar(dnatmp, '.');
	      chopSuffix(shortname);
	      //dnaMixedCaseFilter(dnatmp,dnatmp);     
	      i=hashIntValDefault(speciesHash,shortname,-1);
	      if(i>-1)
		dnabuff[i]=addSuffixFreeHead(dnabuff[i],dnatmp);
	      if(i>-1){
		Rstatus[i]=comp->rightStatus;
		Lstatus[i]=comp->leftStatus;
	      }
	      if(i==0){
		Rstatus[i]='H';
		Lstatus[i]='H';
	      }


	      // fprintf(f,"%d,%s:%d-%d\t",mafNumber,bed->chrom, bedStart+1, bedEnd);
	      // fprintf(f,"%s-%s\t%s\t%f\n",shortname,comp->src,dnatmp,1.0);
	      freeMem(dnatmp);
	      freeMem(shortname);
	    }
	  }	  
	  //mafWrite(f, maf);
	  struct mafAli *nextMaf = (mafEnd > bedEnd+1)
            ? mafSubset(full, mc->src, bedEnd+1, mafEnd) : mafNext(mf); //bedEnd+1
	  if (maf != full)
            mafAliFree(&maf);
	  mafAliFree(&full);
	  maf = nextMaf;
        }
      //fprintf(f,"HG %d,%d,%s:%d-%d",bedOverlap,mafNumber,bed->chrom, bedStart+1, bedEnd);
      //fprintf(f,"%s:%d-%d",bed->chrom, bedStart+1, bedEnd);
      fprintf(f,"#%s\t%d\t%d\n",bed->chrom, bed->chromStart, bed->chromEnd);
      // fprintf(f,"#%s:%d-%d\n",bed->chrom, bedStart+window, bed->chromEnd-window-1);
      //fprintf(f,"%s,%s ",b1,b2);
      //      if(mafNumber>0)
      if(bedOverlap==1){  
	fprintf(stderr,"#OVERLAP %s\t%d\t%d|\t%d\t%d\n",bed->chrom, bed->chromStart, bed->chromEnd, bedStart, bedEnd);
      }

      //if(bedOverlap==0)
      for (i=0;i<numSpecies;i++){ // I can put something to verify if NA
	//pwmScore[i]=scanPwmOneSeqTopHitMyWrapper(dnabuff[i],pwm,&iMax,&sMax); // this frees dnabuff!!
	if(strlen(dnabuff[i])>0){
	  //fprintf(f,"%s\n",dnabuff[i]);
	  stripChar(dnabuff[i],'-');
	  stripChar(dnabuff[i],'.');
	  dnaMixedCaseFilter(dnabuff[i],dnabuff[i]);
	  if(strlen(dnabuff[i])>0){
	    myDna=newDnaSeq(dnabuff[i],dnaFilteredSize(dnabuff[i]),"hg18");
	    pwmScore[i]=scanPwmOneSeqTopHit(myDna,pwm,&iMax,&sMax,&pwmScore2);
	    DNA *dnatmp=cloneStringZ(dnabuff[i]+iMax,pwm->w);
	    toUpperN(dnatmp,pwm->w);
	    if(sMax=='-')
	      reverseComplement(dnatmp,pwm->w);
	    if(pwmScore[i]>-900) //Control that iMax is within the window...
	      fprintf(f,"%d\t%c\t%c\t%d\t%c\t%f\t%s\n",i,Lstatus[i],Rstatus[i],
		        iMax-(bed->chromStart-bedStart),sMax,pwmScore[i],dnatmp);
	    else 
	      fprintf(f,"%d\t%c\t%c\tNA\tNA\tNA\tNA\n",i,Lstatus[i],Rstatus[i]);
	    /*	    if(bedOverlap==1)  
	      fprintf(stderr,"##%d)\t%c\t%c\t%d\t%c\t%f\t%s\n",
		      i,Lstatus[i],Rstatus[i],iMax-window,sMax,pwmScore[i],dnabuff[i]);
	    */
	    freeDnaSeq(&myDna);
	    freeMem(dnatmp);
	  }
	  else{
	    //	      fprintf(f,"%d\t%c\t%c\tMiss\tMiss\tMiss\tMiss\n",i,Lstatus[i],Rstatus[i]);
	    fprintf(f,"%d\t%c\t%c\tNA\tNA\tNA\tNA\n",i,Lstatus[i],Rstatus[i]);
	  }
	}
	else{
	  fprintf(f,"%d\t%c\t%c\tNA\tNA\tNA\tNA\n",i,Lstatus[i],Rstatus[i]);
	}
      }
      

      //      for (i=0;i<numSpecies;i++)
      //	fprintf(f,"\t%f",pwmScore[i]);
      //      fprintf(f,"\n");
      
      //dnaMixedCaseFilter(buff,dnatmp);
      //myDna=newDnaSeq(buff,dnaFilteredSize(buff),"hg18");
      //if(mafNumber>0)
      //	pwmScore=scanPwmOneSeqTopHit(myDna,pwm);
      //fprintf(f,"%f\n",pwmScore);
      //freeDnaSeq(&myDna);
      //freeMem(buff); already freed in myDna
      // Reset to previous maf if regions can overlap... 
      // but I cannot do that, if skeeped....! :( ... 
      // Practical solution: copy previous pwmScore or just eliminate Overlaps...
      
      /* get next region from the bed */
      hashRemove(regionHash, bed->chrom);
      
      //bedOverlap=0;
      if (bed->next){
	/* Check if next region is overlapping with current bed */
	//	if( sameString(bed->chrom,bed->next->chrom) && (bedEnd > bed->next->chromStart))
	//  bedOverlap=1;
        hashAdd(regionHash, bed->chrom, bed->next);
      }
      /*
      if(bedOverlap==0){
	for (i=0;i<numSpecies;i++)
	  freeMem(dnabuff[i]); //Free dnabuff
      }
      else{
	for (i=0;i<numSpecies;i++){
	  //copy the last window 
	  int lenOverlap = (bedEnd - bed->next->chromStart);
	  int lenStr = strlen(dnabuff[i]);
	  if(lenStr>lenOverlap){
	    char *newst=cloneStringZ(dnabuff[i]+(lenStr-lenOverlap),lenOverlap);
	    verbose(1,"##Overlap!%d) %s:%d-%d, %d, %d, %s-->",i,bed->chrom, 
		    bedStart+window, bedEnd-window-1,
		    lenStr,lenOverlap,dnabuff[i]);
	    freeMem(dnabuff[i]);
	    dnabuff[i]=newst;
	    verbose(1,"%s\n",dnabuff[i]);
	  }
	}
      }
      */
    }
  mafFileFree(&mf);
}

/* ******************************************************************************** */
/* ******************************************************************************** */

struct hash *hashWordsInFilePos(char *fileName, int hashSize)
/* Create a hash of space delimited words in file. */
{
struct hash *hash = newHash(hashSize);
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *line, *word;
int pos=0;

while (lineFileNext(lf, &line, NULL))
    {
    while ((word = nextWord(&line)) != NULL)
      hashAddInt(hash, word,pos++);
      //hashAdd(hash, word, NULL);
    }
lineFileClose(&lf);
return hash;
}



void scanPwmOnMafRegions(char *regionFile, char *pwmFile, char *speciesFile, char *out, int mafCount, char *mafFiles[])
/* scanPwmOnMafRegions - scan a Pwm in the multiple aligments of a maf on the regions given by a bed file. */
{
  int i = 0;
  struct hash *bedHash = NULL;

  struct hash *speciesHash = NULL; 
  struct hashEl *hel;
  struct hashCookie cookie;

  FILE *f = NULL;

  struct pssm pwm;

  initialise_pssm(&pwm,pwmFile);

  
  verbose(1, "Extracting from %d files to %s\n", mafCount, out);
  bedHash = loadRegions(regionFile);
  
  f = mustOpen(out, "w");
  verbose(3, "creating %s\n", out);

  verbose(1, "Reading species file %s\n",speciesFile);
  speciesHash = hashWordsInFilePos(speciesFile, 0);
  verbose(1, " Number of species = %d\n",speciesHash->elCount);
  /* Loop through hash */
  //hashIntVal(struct hash *hash, char *name);
  //hashIntValDefault(struct hash *hash, char *name, int defaultInt);
  cookie = hashFirst(speciesHash);
  while ((hel = hashNext(&cookie)) != NULL)
    verbose(2," %ld) %s \n", (long int)(hel->val),hel->name);

  
  for (i = 0; i < mafCount; i++) //Chromosome by chromosome.
    extractMafs(mafFiles[i], f, bedHash,speciesHash,&pwm);
 
  carefulClose(&f);
}


int main(int argc, char *argv[])
/* Process command line. */
{
  optionInit(&argc, argv, options);
  //* keepInitialGaps = optionExists("keepInitialGaps");

 useJaspar = optionExists("j");
 useSnpRobust = optionExists("snp");
 maskRep = optionExists("mask");

 thresh= optionDouble("t",thresh);
 base= optionDouble("base",base);
 addPseudoCounts= optionDouble("p",addPseudoCounts);
 window= optionInt("window",window);

 // ompNumThreads= optionInt("omp",ompNumThreads);
 // omp_set_num_threads(ompNumThreads); //Adjust this to env variable.... or option

  if (argc < 6)
    usage();
  scanPwmOnMafRegions(argv[1],argv[2],argv[3], argv[4], argc-5, &argv[5]);
  return 0;
}


