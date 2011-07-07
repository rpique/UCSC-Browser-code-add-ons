
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <sys/mman.h> 


// The following is set for hg18 
char *chromosomes[]={"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
  "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
  "chr21","chr22","chrX","chrY"};
int chromSizes[]={247249719,242951149,199501827,191273063,180857866,170899992,
                  158821424,146274826,140273252,135374737,134452384,132349534,
                  114142980,106368585,100338915,88827254 ,78774742 ,76117153 ,
                  63811651 ,62435964 ,46944323 ,49691432 ,154913754,57772954 };
#define NUMCHR 24

//#define GENOME_DATA_FOLDER "/Users/rpique/local/motif/ucscData/chromosomes/" 

//#define SUFFIX_SEQ ".bin"
//#define SUFFIX_SCORES ".scores"

#define maxStr 500

#define TEMP_FILE_PAT  "/tmp/rpr.tmp.XXXXXX"

typedef unsigned short int myCutSiteIntType; //2 byte
//typedef unsigned char myCutSiteIntType;

int main(int argc,char *argv[]){

  myCutSiteIntType *CutSites[2][NUMCHR];
  myCutSiteIntType *mmapall;

  int iChr,j,Window;
  int MaxCutSites,strand;
  int aux,left,right;
  int fd;
  int count;
  int skipLine=0;
  long int cF,cR;
  long int offset;
  double offsetCheck;
  char cStrand;
  char buff[maxStr];
  char chr_str[maxStr];
  FILE *inF,*outF;

  //  fprintf(stderr,"# argv[1]=%s\n",argv[1]);

  if(argc !=3){
    fprintf(stderr,"%s Reads.bin Window < SegmentsFile.bed > CutSites.txt\n",argv[0]);
    fprintf(stderr,"SegmentsFile is a tab delimited file with at least three columns: \n");
    fprintf(stderr," Chr1\t1000100\t1000150\t[Other Columns]\n");
    fprintf(stderr," .... \n"); 
    fprintf(stderr," ChrY\t1000100\t1000150\t[Other Columns]\n");
    fprintf(stderr,"Window is the size in bases around the center of the Segments, i.e., \n (start+end)/2 +/- window, in which the cutsites are retrieved.\n");
    fprintf(stderr,"The output in stdout is in the following form: \n");
    //    fprintf(stderr," Chr1\t1000050\t1000250\t[Other Columns]\tNumber of reads with cutsite in 1000050,1000051,... 1000250\n");
    fprintf(stderr,"Instance1: Forward reads, Reverse strand reads.");
    fprintf(stderr," .... \n\n");    
    return 0;
  }

  //fprintf(stderr,"# Creating binary file for memory mapping...\n");
  aux=0;
  offset=0;
  outF=fopen(argv[1],"r+b");
  assert(outF!=NULL);

  for(strand=0;strand<2;strand++)
    for(iChr=0;iChr<NUMCHR;iChr++){
      offset=offset+chromSizes[iChr]*sizeof(myCutSiteIntType);
      offsetCheck=offsetCheck+chromSizes[iChr]*sizeof(myCutSiteIntType);
    }


  fd=fileno(outF);
  fprintf(stderr,"# Mapping %s with %ld=%g bytes in memory\n",argv[1],offset,offsetCheck);
  
  mmapall=(myCutSiteIntType *)mmap(0,offset,PROT_READ,MAP_PRIVATE,fd,0);   
  assert(mmapall!=NULL);

  //mmapall=(myCutSiteIntType *)malloc(offset);
  // //read(fd,mmapall,offset);
  // cF=fread(mmapall,1,offset,outF);
  // fprintf(stderr,"# Bytes read %ld\n",cF);

  offset=0;
  for(strand=0;strand<2;strand++)
    for(iChr=0;iChr<NUMCHR;iChr++){
      CutSites[strand][iChr]=&(mmapall[offset]);
      //CutSites[strand][iChr]=(mmapall+offset*sizeof(myCutSiteIntType));
      //CutSites[strand][iChr]=(myCutSiteIntType *)mmap(NULL,chromSizes[iChr]*sizeof(myCutSiteIntType),
      //					      PROT_READ,MAP_PRIVATE,fd,offset*sizeof(myCutSiteIntType));
      offset=offset+chromSizes[iChr];
      //assert(CutSites[strand][iChr]!=NULL);
      //fprintf(stderr,"# Chri%d Str%d Address %p, %ld\n",iChr,strand,CutSites[strand][iChr],CutSites[strand][iChr]-CutSites[0][0]);
    }
  fprintf(stderr,"# offset %ld \n",offset);

  /*
  cR=0;cF=0;
  for(iChr=0;iChr<NUMCHR;iChr++)
    for(j=0;j<chromSizes[iChr];j++){
      cF+=CutSites[1][iChr][j];
      cR+=CutSites[0][iChr][j];
    }
  fprintf(stderr,"# Forw reads %ld, Revr read %ld\n",cF,cR);
  */   

  count=1;
  fprintf(stderr,"# Opening segment files %s\n","stdin"); //argv[2]);
  Window=atoi(argv[2]);
  fprintf(stderr,"# Window=%d\n",Window);
  inF=stdin; //fopen(argv[2],"r");
  //fprintf(stderr,"# Storing segment files %s\n",argv[3]);
  //outF=fopen(argv[3],"wb");

  cF=0;cR=0;
  while(!feof(inF)){
    fscanf(inF,"%s\t%d\t%d\t%c\%[^\n]\n",chr_str,&left,&right,&cStrand,buff);
    for(iChr=0;iChr<NUMCHR;iChr++)
      if(strcmp(chromosomes[iChr],chr_str)==0)
	break;
    if(iChr<NUMCHR){
      //fprintf(stderr,"#%s\t%d\t%d\t%c\n",chromosomes[iChr],left,right,cStrand);
      //aux=(left+right)/2;
      //assert(aux<=chromSizes[iChr]);
      
      left=left-Window;
      right=right+Window;
     
      if(left>=0 && right<chromSizes[iChr]){
	//fprintf(stdout,"%s\t%d\t%d\t%s",chromosomes[iChr],left,right,buff);
	//fwrite(&(CutSites[iChr][left]),sizeof(int),right-left+1,outF);
	if(cStrand=='+'){
	  for(j=left;j<=right;j++)
	    fprintf(stdout,"\t%d",CutSites[1][iChr][j]);
	  for(j=left;j<=right;j++)
	    fprintf(stdout,"\t%d",CutSites[0][iChr][j]);
	  fprintf(stdout,"\n");
	  for(j=left;j<=right;j++)
            cF+=CutSites[1][iChr][j];
          for(j=left;j<=right;j++)
            cR+=CutSites[0][iChr][j];

	}
	else if(cStrand=='-'){
	  for(j=right;j>=left;j--)
	    fprintf(stdout,"\t%d",CutSites[0][iChr][j]);
	  for(j=right;j>=left;j--)
	    fprintf(stdout,"\t%d",CutSites[1][iChr][j]);
	  fprintf(stdout,"\n");
	  for(j=left;j<=right;j++)
            cF+=CutSites[1][iChr][j];
          for(j=left;j<=right;j++)
            cR+=CutSites[0][iChr][j];
	}
	else{
  	 fprintf(stderr,"# Unrecognized strand parameter!!!\n");	
	 skipLine=1;
	}  
      }
      else{
	fprintf(stderr,"# Skipping segment off-limits\n");
	skipLine=1;
      }
      
      if((count%100)==0)
	fprintf(stderr,"# Processed lines: %d\r",count);
    }
    else{
      fprintf(stderr,"# Unreq %s:%d-%d in line %d !!!\n",chr_str,left,right,count);
      skipLine=1;
    }
    
    if(skipLine==1){
      /* for(j=left;j<=right;j++) */
      /* 	fprintf(stdout,"\t%d",0); */
      /* for(j=left;j<=right;j++) */
      /* 	fprintf(stdout,"\t%d",0); */
      fprintf(stdout,"NA\n");
      skipLine=0;
    }


    count++;
  }
  fclose(inF);
  fprintf(stderr,"# Mapped %ld on F, %ld on R strands \n",cF,cR);
  fprintf(stderr,"# Succesfully processed %d segments\n",count);
  return 0;
}
