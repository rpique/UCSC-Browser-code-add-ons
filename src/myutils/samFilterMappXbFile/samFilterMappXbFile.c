/* samFilterMappXbFile - Filter for mappability sam aligned reads using and mappability Xb file. */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"

#include "khash.h"
#include "xbfiles.h"

//KHASH_MAP_INIT_STR(hashChr_t, unsigned char)

void usage()
/* Explain usage and exit. */
{
errAbort(
  "samFilterMappXbFile - Filter for mappability sam aligned reads using an Xb file and seting quality to 0\n"
  "usage:\n"
  "   samFilterMappXbFile mapp.x8b input.sam/bam output.sam/bam\n"
  " * supports stdin for input.sam and stdout for output.sam\n"
  " * input/output is uncompressed i.e., sam but can be piped to samtools \n"
  // "options:\n"
  // "   -readSize=20\n"
  );
}

static struct optionSpec options[] = {
   {NULL, 0},
};



void samFilterMappXbFile(char *xbFileName, char *samFileName, char *outFileName)
/* samFilterMappXbFile - Filter for mappability sam aligned reads using and mappability Xb file. */
{
  struct lineFile *lf = lineFileOpen(samFileName, TRUE);

  FILE *outF=mustOpen(outFileName,"w");

  xbList_t *xbl=xbLoadMmap(xbFileName);

  int skipLine=0;
  //long int cF,cR;

  int iChr,j;
  int left;//,right;
  int count=0;
  char cStrand;
  int flag;
  int mappState;
  int readLen;
  

  char *row[1000];  // I don't know what is an appropriate number!
  int wordCount;
  char *chr_str;
  

  khiter_t k;
  khiter_t khit;
  int hret;

  khash_t(hashChr_t) *hChr= kh_init(hashChr_t); 
  
  //Hash chromosome names
  for(iChr=0;iChr<xbl->count;iChr++){
    k = kh_put(hashChr_t, hChr , cloneString(xbl->names[iChr]), &hret);
    assert(hret==1);
    kh_val(hChr, k) = iChr;
  }
  
  // cF=0;cR=0;
  while ((wordCount = lineFileChop(lf, row)) != 0){
    if(row[0][0]!='@'){
    //if(wordCount>=10){
      chr_str = row[2];
      left = lineFileNeedNum(lf, row, 3)-1;
      //    right = lineFileNeedNum(lf, row, 2);
      flag = lineFileNeedNum(lf, row, 1);
      readLen = strlen(row[9]);
    
      cStrand='+';
      if(flag|0x0010)
	cStrand='-';
    
      //    for(j=0;j<wordCount;j++)
      //  fprintf(outF,"%s\t",row[j]);

      /* SHOULD I DO SOMETHING ABOUT STRAND?*/
      /* No; because I use the first position or the forward even if the match is in the reverse. */
      khit = kh_get(hashChr_t, hChr , chr_str);
      if(kh_exist(hChr,khit)){
	iChr=kh_val(hChr,khit);
	if((left>0) && (left<=xbl->sizes[iChr])){
	  mappState=0;
	  //fprintf(stdout,"|%d,%s,%d,%d,%s|",count,chr_str,left,readLen,row[4]);
	  for(j=left;j<=(left+readLen-20);j++){
	    if(xbl->vec[iChr].a[j]==1)
	      mappState++;
	    //fprintf(stdout,"%d,",xbl->vec[iChr].a[j]);
	  }
	  //fprintf(stdout,".%d\n",mappState);
	  // PRINT sam line
	  if(mappState<1)
	    sprintf(row[4],"%d",mappState);
	  //if(mappState!=1)
	  //  sprintf(row[4],"%s","X");
	  for(j=0;j<(wordCount-1);j++)
	    fprintf(outF,"%s\t",row[j]);
	  fprintf(outF,"%s\n",row[j]);
	}else{
	  verbose(2,"# Skipping segment off-limits\n");
	  //verbose(1,"# Unreq %s:%d-%d in line %d !!!\n",chr_str,left,right,count);
	  skipLine=1;
	}
      }else{
	verbose(2,"# Skipping segment off-limits\n");
	//verbose(1,"# Unreq %s:%d-%d in line %d !!!\n",chr_str,left,right,count);
	skipLine=1;
      }  
    }else{
      // This should help to pass through the header part of the bam/sam file
      for(j=0;j<(wordCount-1);j++)
	fprintf(outF,"%s\t",row[j]);
      fprintf(outF,"%s\n",row[j]);
    }
    if(verboseLevel()==2){
      if((count%100000)==0){
	verboseDot();
	if((count%8000000)==0)
	  verbose(2,"\n");
      }
    }
    if(skipLine==1){
      //fprintf(outF,"NA\n");
      skipLine=0;
    }    
    count++;
  }
  
  lineFileClose(&lf);

  carefulClose(&outF);
  //verbose(1,"# Mapped %ld on F, %ld on R strands \n",cF,cR);
  verbose(1,"# Succesfully processed %d segments\n",count);

  //free xbl?
}


int main(int argc, char *argv[])
/* Process command line. */
{
  optionInit(&argc, argv, options);
  if (argc != 4)
    usage();

  //  window = optionInt("window", window);
  //  readsize = optionInt("readsize", readsize);
  
  samFilterMappXbFile(argv[1],argv[2],argv[3]);
  return 0;
}



