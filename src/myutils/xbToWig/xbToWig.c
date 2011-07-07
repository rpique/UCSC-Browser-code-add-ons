/* xbToWig - Converts a xb file to a fixed step wig file. */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"

#include "xbfiles.h"

static char const rcsid[] = "$Id: newProg.c,v 1.30 2010/03/24 21:18:33 hiram Exp $";

void usage()
/* Explain usage and exit. */
{
  errAbort(
  "xbToWig - Converts a xb file to a fixed step wig file\n"
  "usage:\n"
  "  xbToWig in.xb out.F.wig out.R.wig\n"
  "  xbToWig -combined in.xb out.wig\n"
  "in.xb - is the xb binary file \n"
  "out.F.wig - wig file with the forward strand of the xb file \n"
  "out.R.wig - wig file with the reverse strand of the xb file \n"
  "options:\n"
  "   -smoothing = 0 (default) no smoothing, \n"
  "                1 -> (0.25 0.5 0.25)      \n" //, 2 - (1/9,2/9,3/9,2/9,1/9)"
  "   -combined -> if set, combine + and - strand into a single output\n"
  "   -verbose = verbosity lebel (1 default)\n"
  );
}

boolean combined=FALSE;
int smoothing=0;

static struct optionSpec options[] = {
   {"combined", OPTION_BOOLEAN},
   {"smoothing", OPTION_INT},
   {NULL, 0},
};

void xbToWig(char *xbFile, char *outFrFile, char *outRvFile)
/* xbToWig - Converts a xb file to a fixed step wig file. */
{
  xbList_t *xbl=xbLoadMmap(xbFile);
  FILE *oFf=mustOpen(outFrFile,"w");
  FILE *oRf=mustOpen(outRvFile,"w");
  int cR,cF;
  int iChr,j;
  
  //Write chromsizes ?
  
  cR=0;cF=0;
  for(iChr=1;iChr<xbl->count;iChr++){
    verbose(1,"# Processing %s\n",xbl->names[iChr]);
    fprintf(oFf,"variableStep chrom=%s\n",xbl->names[iChr]);
    fprintf(oRf,"variableStep chrom=%s\n",xbl->names[iChr]);
    for(j=0;j<xbl->sizes[iChr];j++){
      cF+=xbl->vec[(iChr<<1)].a[j];
      cR+=xbl->vec[(iChr<<1)+1].a[j];
      if((xbl->vec[(iChr<<1)].a[j])>0)
	fprintf(oFf,"%d\t%d\n",j+1,xbl->vec[(iChr<<1)].a[j]);
      if((xbl->vec[(iChr<<1)+1].a[j])>0)
	fprintf(oRf,"%d\t%d\n",j+1,xbl->vec[(iChr<<1)+1].a[j]);
    }
  }
  verbose(1,"# Forw reads %d, Revr read %d\n",cF,cR);
  
  fclose(oFf);
  fclose(oRf);
}

void xbToWigCombined(char *xbFile, char *outFile)
/* xbToWig - Converts a xb file to a fixed step wig file. */
{
  xbList_t *xbl=xbLoadMmap(xbFile);
  FILE *of=mustOpen(outFile,"w");
  int cR,cF;
  int iChr,j,k;
  int sum;
  float filter[30]; // global??
  float res=0.0;
  
  if(smoothing==1){
    filter[0]=0.25;
    filter[1]=0.5;
    filter[2]=0.25;
  }
  if(smoothing==2){
    filter[0]=1/(float)16;
    filter[1]=4/(float)16;
    filter[2]=6/(float)16;
    filter[3]=4/(float)16;
    filter[4]=1/(float)16;
  }
  if(smoothing==3){
    filter[0]=1/(float)64;
    filter[1]=6/(float)64;
    filter[2]=15/(float)64;
    filter[3]=20/(float)64;
    filter[4]=15/(float)64;
    filter[5]=6/(float)64;
    filter[6]=1/(float)64;
  }
  if(smoothing==4){
    filter[0]=1/(float)256;
    filter[1]=8/(float)256;
    filter[2]=28/(float)256;
    filter[3]=56/(float)256;
    filter[4]=70/(float)256;
    filter[5]=56/(float)256;
    filter[6]=28/(float)256;
    filter[7]=8/(float)256;
    filter[8]=1/(float)256;
  }
    
  cR=0;cF=0;
  for(iChr=1;iChr<xbl->count;iChr++){
    verbose(1,"# Processing %s smoothing=%d\n",xbl->names[iChr],smoothing);
    fprintf(of,"variableStep chrom=%s\n",xbl->names[iChr]);
    for(j=0;j<xbl->sizes[iChr];j++){
      cF+=xbl->vec[(iChr<<1)].a[j];
      cR+=xbl->vec[(iChr<<1)+1].a[j];
      if(smoothing==0){
	sum=(xbl->vec[(iChr<<1)].a[j]+xbl->vec[(iChr<<1)+1].a[j]);      
	if(sum>0)
	  fprintf(of,"%d\t%d\n",j+1,sum);
      }
      else{
	res=0.0;
	for(k=0;k<((smoothing*2)+1);++k)
	  res = res + filter[k]*(float)(xbl->vec[(iChr<<1)].a[j-k+smoothing] + xbl->vec[(iChr<<1)+1].a[j-k+smoothing]);
	if(res>0.0001)
	  fprintf(of,"%d\t%lf\n",j+1,res);
      }
    }
    verbose(1,"# Forw reads %d, Revr read %d\n",cF,cR);
    fflush(of);
  }
  verbose(1,"# Forw reads %d, Revr read %d\n",cF,cR);
  
  fclose(of);
}


int main(int argc, char *argv[])
/* Process command line. */
{
  optionInit(&argc, argv, options);
  combined = optionExists("combined");
  smoothing = optionInt("smoothing", smoothing);
  
  if(!combined){
    if (argc != 4)
      usage();
    xbToWig(argv[1],argv[2],argv[3]);
  }
  if(combined){
    if (argc != 3)
      usage();
    xbToWigCombined(argv[1],argv[2]);
  }

  return 0;
}
