#include "common.h"
#include "pwm.h"


//void initialise_pssm(struct pssm *pwm,char *fileMotif)
void initialise_pssm(struct pssm *pwm,char *fileMotif,double addPseudoCounts, boolean useJaspar) 
{
  int n,b,NumZeros=0;
  double FirstRowCounts=0.0,RowCounts=0.0,MinNonZeroCount=1.0;

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
    for(b=0;b<4;b++){
      if(pwm->matrix[n][b]==0.0){
	NumZeros++;
      }else{	
	if(pwm->matrix[n][b]<MinNonZeroCount)
	  MinNonZeroCount=pwm->matrix[n][b];
      }
    }
    verbose(2,"# Row %d: %g Counts\n",n,RowCounts);
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

  if(addPseudoCounts<=0){
    if(NumZeros>0)
      if(FirstRowCounts<1.1)
	addPseudoCounts=MinNonZeroCount/2; //FirstRowCounts/100; 
    if(NumZeros==0)
      addPseudoCounts=0.0;  //Not necessary to add pseudo counts...
  }

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
  if(verboseLevel()==2){
    printMatrix(pwm);  
    convertPSSMToLogs(pwm);
    printMatrix(pwm); 
  }
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
  //printMatrix(pm);  
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
int ATGCbase(char b, boolean maskRep) {
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
double compare_subseq_to_pssm (char *seq, struct pssm *p, boolean useSnpRobust) {
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

  return score; ///log(base);  // 
}
