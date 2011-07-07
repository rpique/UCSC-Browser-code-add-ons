/* bigWigExtractBedRegions - Extract data in big wig corresponding to segments defined in a Bed file. */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "sqlNum.h"
#include "udc.h"
#include "bigWig.h"

#include "dnaseq.h"
#include "twoBit.h"


static char const rcsid[] = "$Id: newProg.c,v 1.28 2009/07/02 17:14:39 angie Exp $";

void usage()
/* Explain usage and exit. */
{
errAbort(
  "bigWigExtractBedRegions - Extract data in big wig corresponding to segments defined in a Bed file\n"
  "usage:\n"
  "   bigWigExtractBedRegions input.bw bedregions.bed output.txt\n"
  "options:\n"
  "   -Window=0 Window around the bed\n"
  "   -udcDir=/dir/to/cache - place to put cache for remote bigBed/bigWigs\n"
  );
}

int Window = 0;	/* Window size from command line. */
char *summaryType = "mean";


static struct optionSpec options[] = {
   {"Window", OPTION_INT},
   {"udcDir", OPTION_STRING},
   {NULL, 0},
};


void bigWigSummary(struct bbiFile *bwf,FILE *outFile, char *chrom, int start, int end, 
		   char cStrand, int dataPoints)
/* bigWigSummary - Extract summary information from a bigWig file.. */
{

  /* Make up values array initialized to not-a-number. */
  //double nan0 = strtod("NaN", NULL);
  double summaryValues[dataPoints];
  int i;
  for (i=0; i<dataPoints; ++i)
    summaryValues[i] = 0; //nan0;
  
  bigWigSummaryArray(bwf, chrom, start, end, bbiSummaryTypeFromString(summaryType), 
		     dataPoints, summaryValues);

  if (cStrand == '+'){
    for (i=0; i<dataPoints; ++i){
      double val = summaryValues[i];
      if (i != 0)
	fprintf(outFile,"\t");
      if (isnan(val))
	fprintf(outFile,"NaN");
      else
	fprintf(outFile,"%g", val); 
    }
    fprintf(outFile,"\n");
  }
  if (cStrand == '-'){
    for (i=(dataPoints-1);i>=0; --i){
      double val = summaryValues[i];
      if (i != (dataPoints-1))
	fprintf(outFile,"\t");
      if (isnan(val))
	fprintf(outFile,"NaN");
      else
	fprintf(outFile,"%g", val); 
    }
    fprintf(outFile,"\n");
  }
//fprintf(stderr,"# No data in region %s:%d-%d in BigWig file\n", chrom, start, end);
}



void bigWigExtractBedRegions(char *bigWigFile,char *bedFileName, char *outName)
/* bigWigExtractBedRegions - Extract data in big wig corresponding to segments defined in a Bed file. */
{
  struct bbiFile *bwf = bigWigFileOpen(bigWigFile);
  FILE *inFile = mustOpen(bedFileName, "r");
  FILE *outFile = mustOpen(outName, "w");

  char chr_str[100];
  int left;
  int right;
  char cStrand;
  char buff[500];


  while(!feof(inFile)){
    fscanf(inFile,"%s\t%d\t%d\t%c\%[^\n]\n",chr_str,&left,&right,&cStrand,buff);
    left=(left-Window);
    right=(right+Window);
    
    bigWigSummary(bwf,outFile,chr_str,left,right,cStrand,right-left);
  }

  // bigWigSummary(argv[1], argv[2], sqlUnsigned(argv[3]), sqlUnsigned(argv[4]), sqlUnsigned(argv[5]));

  carefulClose(&outFile);
  carefulClose(&inFile);
  bigWigFileClose(&bwf);
}


int main(int argc, char *argv[])
/* Process command line. */
{
  optionInit(&argc, argv, options);
  if (argc != 4)
    usage();
  Window = optionInt("Window", Window);
  //  dnaUtilOpen();
  // summaryType = optionVal("type", summaryType);
  udcSetDefaultDir(optionVal("udcDir", udcDefaultDir()));

  bigWigExtractBedRegions(argv[1], argv[2],argv[3]);
  return 0;
}

