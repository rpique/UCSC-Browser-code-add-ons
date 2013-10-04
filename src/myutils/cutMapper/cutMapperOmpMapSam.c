/* cutMapperOmpMapSam - Maps reads in FastQ file using and indexed genome and makes output into jfd format. */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"

#include <omp.h>

#include "dnautil.h"

#include "cutMapperCore.h"

//#include "cutMapperResampler.h"

static char const rcsid[] = "$Id: newProg.c,v 1.30 2010/03/24 21:18:33 hiram Exp $";

void usage()
/* Explain usage and exit. */
{
  errAbort(
  "cutMapperOmpMapSam - Maps reads in FastQ file using and indexed genome and makes output into jfd format\n"
  "usage:\n"
  "   cutMapperOmpMapSam indexFolder file.fastQ mappedreads.txt\n"
  "options:\n"
  //  "   -sensitivity=folder --> Sensitivity folder containing the files, default not used?, \n"
  "                sensitivity files should be located as <folder/chr???.sensitivites.txt \n"
  "   -omp=1 - multi-core implementation (default 1 == one thread)"
  "   -verbose=1,2,3 (default 1) \n"
  );
}

//boolean applySensitivity=FALSE;
//char *sensitivityFolder;
//double **sens;

int ompNumThreads= 1;

boolean rmdup=FALSE;


static struct optionSpec options[] = {
  //  {"sensitivity", OPTION_STRING},
  {"omp" , OPTION_INT},
  {"rmdup", OPTION_BOOLEAN},
  {NULL, 0},
};


void reParseFastqFileAndMap(char *inFastq, khash_t(hashPos_t) **hReads, FILE *f, char **chromNames)
/* hashFastqFile - Parse file and map reads using the hash and repeat table */
{
  //  FILE *f=mustOpen(outFile,"w");
  //  FILE *fR=mustOpen(repFile,"w");
  
  struct lineFile *lf = lineFileOpen(inFastq, TRUE);
  char *seqName, *seq, *seqString, *qName, *qual, *seqQual;
  int line = 0;
  //int numKmers = 0;
  int numUni = 0;
  int numRep = 0;
  int numSeq = 0;
  int numUnm = 0;
  int numWithN=0;
  boolean startOfFile = TRUE;
  
  unsigned long int kmer;
  unsigned char Prefix=0;
  unsigned int Suffix;
  
  khiter_t khit;
  //int hret;
  int j;
  
  int rloc=0;

  j=0;

  verbose(2,"[%s %3d] file(%s)\n", __func__, __LINE__, inFastq);
  while ( lineFileNext(lf, &seqName, NULL)){
    ++line;
    if (startOfFile){
      if (*seqName == '#'){
	//free seqName...???
	continue;
      }else{
    	  startOfFile = FALSE;
      }
    }
    seqName=cloneString(seqName);
    wantNewLine(lf, inFastq, ++line, &seq, "fastq sequence line");
    seqString=cloneString(seq);
    kmer=packKmer(seq,kmerSize);
    wantNewLine(lf, inFastq, ++line, &qName, "fastq sequence name (quality line)");
    wantNewLine(lf, inFastq, ++line, &qual, "quality line");
    seqQual=cloneString(qual);
    //Cut the strings to 20-mers 
    seqQual[kmerSize]='\0';
    seqString[kmerSize]='\0';
    verbose(3,"[%s %3d] Mapping Kmer %lx: %s\t%s\n", __func__, __LINE__,kmer,seqName+1,seqString);
    if(kmer != -1){ // otherWise check how many N, and hash them in a different way... ? 4 times??
      Prefix=getPrefix(kmer); // ((ForKmer&kmermask)>>32);
      Suffix=getSuffix(kmer); // &0xFFFFFFFF;
      khit = kh_get(hashPos_t, hReads[Prefix] , Suffix);
      verbose(3,"[%s %3d] Mapping Kmer %lx: %x\t%x (%lx)\n", __func__, __LINE__,kmer,Prefix,Suffix,khit);
      if ( kh_exist(hReads[Prefix], khit)){
	//In the hash... OUTPUT READ
	verbose(3,"[%s %3d] Mapping Kmer %lx: %x\t%x (%lx) %d:%d\t%s\n", __func__, __LINE__,kmer,Prefix,Suffix,khit,(int)kh_value(hReads[Prefix], khit).chr,kh_value(hReads[Prefix], khit).pos,seqString);
	if(kh_value(hReads[Prefix], khit).chr==255){//Non uniquely mappable //255 Use 254 for debugging
	  numRep++;
	  verbose(3,"[%s %3d] Reapeating 255 Kmer %lx: %x\t%x (%lx)\n", __func__, __LINE__,kmer,Prefix,Suffix,khit);
	}
	else{ //Uniquely mappable
	  //unPackKmer((((unsigned long int)j)<<32) + Suffix, 20, kmerString);
	  numUni++;
	  rloc=kh_value(hReads[Prefix], khit).pos;
	  if(rloc<0){
	    rloc = rloc - kmerSize + (int)strlen(seqString);
	    reverseBytes(seqQual,strlen(seqString));
	    reverseComplement(seqString,strlen(seqString));
	  }
	  // assert(kh_value(hReads[Prefix], khit).chr < 95 );
	  // assert(kh_value(hReads[Prefix], khit).chr > 0  );

	  fprintf(f,"%s\t%d\t%s\t%d\t%d\t%dM\t*\t0\t0\t%s\t%s\tX0:i:1\n",
		  seqName+1, //1) Query template NAME
		  (rloc<0) ? 16 : 0, //2) bitwise FLAG
		  chromNames[kh_value(hReads[Prefix], khit).chr],    //3) Reference sequence NAME
		  (int) abs(rloc),           //4) 1-based leftmost mapping POSition
		  20,                //5) Quality
		  (int)strlen(seqString), //6) CIGAR string
		  //7) Ref Name of the mate/next segment 
	          //8) Position of the mate/next segment
		  //9) Observed Template LENgth
		  seqString, //10) Sequence
		  seqQual      //11) ASCII of Phred-scaled base QUALity+33
		  //12) Optional 
		  );

	  verbose(3,"[%s %3d] Kmer %lx: %s\t%d\t%d\t%d\t%s\n", __func__, __LINE__,kmer,seqName+1,(rloc<0) ? 16 : 0, (int) kh_value(hReads[Prefix], khit).chr,(int) abs(rloc),seqString);
	  
	  // fprintf(f,"%d\n",kh_value(hReads[Prefix], khit).pos);

	}// else unmappable
	//fprintf(f,"");
      	if(rmdup){ // This will remove the from the hash.
      	  kh_del(hashPos_t, hReads[Prefix], khit);
      	  //verbose(3,"[%s %3d] Remove Kmer %lx from hash\n", __func__, __LINE__,kmer);
	}
      }
      else{
	//UnMappable!!
	verbose(3,"[%s %3d] Unmappable2 Kmer %lx: %s\t%d\t%d\t%d\t%s\n", __func__, __LINE__,kmer,seqName+1,(rloc<0) ? 16 : 0, (int) kh_value(hReads[Prefix], khit).chr,(int) abs(rloc),seqString);
	numUnm++;
      }
    }
    else{
      //unMappable with Ns
      numWithN++;
    }
    ++numSeq;
    freeMem(seqName);
    freeMem(seqString);
    freeMem(seqQual);
  }

  verbose(1,"#Summary2:\t%s\t%d\t%d\t%d\t%d\t%d\t%f\t%f\n"
	  ,inFastq,numUni,numRep,numUnm,numWithN,numSeq
	  ,numUnm/(float)(numUni+numRep+numUnm)*100
	  ,(numUnm+numWithN)/(float)(numSeq)*100);  
  verbose(1,"#\tUnique %d, Repeats <10 %d reads\n",numUni,numRep); 
  verbose(1,"#\tprocessed %d lines  %d reads\n",line,numSeq);
  verbose(1,"# Total # reads with >0 Ns: %d\n",numWithN);
  
  //  verbose(1, "#\tTotal number of k-mers %d, %d \n",line,numKmers); 
  //verbose(1, "#\tTotal number of rep k-mers %d \n",numRepKmers); 

  //  FILE *f=mustOpen(outFile,"w");
  //  carefulClose(&f);
  
  //  kh_destroy(hashPos_t, h);
  
  //close lf file...
  lineFileClose(&lf);
}

void cutMapperOmpMapSam(char *inFastq, char *indexFolder, char *uniFile)
/* cutMapperOmpMapSam - Maps reads in FastQ file using and indexed genome and makes output into jfd format. */
{
  FILE *f;
  char cbuff[strlen(indexFolder)+256];
  //char kmerString[21];
  
  unsigned int j,k;
  //int i;
  int TotalMappedReads=0;
  int UniquelyMappedReads=0;
  int UnMappableReads=0;
  int LessThan10Reps=0;

  unsigned int rloc;
  
  clock_t t;
  clock_t t0=clock();
    
  khash_t(hCount_t) *hReads[256];
  
    
  khash_t(hashPos_t) *hMapped[256]; 

  khash_t(hashChr_t) *hChr; 
  char **chromNames;
  unsigned *chromSizes;

  /* --------------------------------------------------------------------- */

  verbose(1,"Loading Chromosome sizes\n");
  sprintf(cbuff,"%s/chromSizes.txt",indexFolder);
  hChr = loadChromSizes(cbuff, &chromNames, &chromSizes);
  for(j=0;j<kh_end(hChr);j++)
    if(kh_exist(hChr,j)){
      rloc=kh_val(hChr,j);
      verbose(2,"%d %s %s %d\n",rloc,kh_key(hChr,j),chromNames[rloc],chromSizes[rloc]);
    }
  verbose(1,"#Loaded %d chromosomes\n",kh_size(hChr));
  //Initializing hash maps, for mapped reads!
  for(j=0;j<256;j++) 
    hMapped[j]=kh_init(hashPos_t);   
  //Initializing repeat maps, for mapped reads
  
  //Initializing hash maps.
  for(j=0;j<256;j++) 
    hReads[j]=kh_init(hCount_t);   

  /* --------------------------------------------------------------------- */
  //  if(applySensitivity)
  //   sens=read_sens(chromNames,kh_size(hChr),sensitivityFolder);
  
  /* --------------------------------------------------------------------- */

  //Hashing input reads.
  t = clock();
  hashFastqFile(inFastq, hReads);
  verbose(1,"# Hashed fastq reads in %f seconds (%f total)\n",
	  (double)(clock() - t)/CLOCKS_PER_SEC, 
	  (double)(clock() - t0)/CLOCKS_PER_SEC); 
  
  /* --------------------------------------------------------------------- */
  
  //for(j=57;j<59;j++){ //256  
  // USE MPI or OPENMP HERE?????
#pragma omp parallel for private(t,j,k) reduction(+:LessThan10Reps,UniquelyMappedReads,UnMappableReads) schedule(dynamic,1) 
  for(j=0;j<256;j++){ //256
    int hret;
    char cbuff2[strlen(indexFolder)+256];

    char bluff[40];

    khint_t khit,khit2,indexSize;

    khash_t(hashPos_t) *hIndex;
    //    repvec_t *repvecp;  
    //    replist_t replist;

    //Load Hash! 
    sprintf(cbuff2,"%s/Hash%03d.bin",indexFolder,j);
    t = clock();
    hIndex = hashFileOpen(cbuff2);
    verbose(1,"# Opened hash for index %d in %f seconds (%f total)\n",j,
	    (double)(clock() - t)/CLOCKS_PER_SEC,
	    (double)(clock() - t0)/CLOCKS_PER_SEC);    

    /*
    for(k = kh_begin(hIndex); k < kh_end(hIndex); ++k){
      if(kh_exist(hIndex, k)){ // I found it!
	unsigned int Suffix = kh_key(hIndex, k);
	if(kh_value(hIndex, k).chr<255 & (Suffix==0x8ee01093)){
	  unPackKmer((((unsigned long int)j)<<32) + (unsigned long int)Suffix, 20, bluff);
	  verbose(3,"[%s %3d] %x\t%x (%u) %d:%d\t %s\n", __func__, __LINE__,j,Suffix,k,
		  (int)kh_value(hIndex, k).chr, kh_value(hIndex, k).pos,bluff);
	  assert(kh_value(hIndex, k).chr>0);
	  assert(kh_value(hIndex, k).chr<96);
	}
      }
    }
    */
    
    
    /* --------------------------------------------------------------------- */

    t = clock();
    //Traverse the count hash
    indexSize=kh_size(hIndex);
    for (k = kh_begin(hReads[j]); k < kh_end(hReads[j]); ++k){
      if (kh_exist(hReads[j], k)){
	unsigned int Suffix=kh_key(hReads[j], k);
	unsigned int num=kh_val(hReads[j], k);
	// Look up the key, on the Suffix hash table hIndex
	khit = kh_get(hashPos_t, hIndex , Suffix);
	if(kh_exist(hIndex, khit)){// && (kh_end(hIndex) < khit)){ // I found it!
	  TotalMappedReads+=num;
	  if(kh_value(hIndex, khit).chr==255){//Non uniquely mappable
	    LessThan10Reps+=num;
	  }
	  else{ // Uniquely mappable
	    UniquelyMappedReads+=num;
	    khit2 = kh_put(hashPos_t, hMapped[j], Suffix, &hret);
	    assert(hret==1);
	    kh_val(hMapped[j],khit2) = kh_val(hIndex, khit);
	    //unPackKmer((((unsigned long int)j)<<32) + (unsigned long int)Suffix, 20, bluff);
	    //verbose(3,"[%s %3d] %x\t%x (%u,%u) %d:%d\t %s\n", __func__, __LINE__,j,Suffix,khit,khit2,
	    //	    (int)kh_value(hMapped[j], khit2).chr,kh_value(hMapped[j], khit2).pos,bluff);
	    assert(kh_value(hIndex, khit).chr==kh_val(hMapped[j],khit2).chr);
	    assert(kh_value(hIndex, khit).pos==kh_val(hMapped[j],khit2).pos);
	    assert(kh_value(hMapped[j], khit2).chr < kh_size(hChr));
	    assert(kh_value(hMapped[j], khit2).chr > 0);
	  }
	}
	else{
	  UnMappableReads+=num;
	}
      }
    }/* index k for*/

    /* --------------------------------------------------------------------- */

    
    verbose(1, "# Mapped index %d: Cum. total mapped %d,"
	    " Cum Total unique %d, Cum Total <10rep %d,"
	    " Cum Total Unmappable %d (%f %%)\n"
	    ,j,TotalMappedReads,UniquelyMappedReads
	    ,LessThan10Reps,UnMappableReads
	    ,UnMappableReads/(float)(UnMappableReads+TotalMappedReads)*100);
    verbose(1,"#\tIndex %d: Mapped %d kmers in %f seconds (%f total)\n",
	    j,indexSize,
	    (double)(clock() - t)/CLOCKS_PER_SEC,
	    (double)(clock() - t0)/CLOCKS_PER_SEC);
    //DESTROYING OBJECTS 

    kh_destroy(hashPos_t,hIndex);    
  } /* Index j for */   
  verbose(1,"# Total time: %f seconds = %f minutes\n",
	  (double)(clock() - t0)/CLOCKS_PER_SEC,
	  (double)(clock() - t0)/CLOCKS_PER_SEC/60);

  
  for(j=0;j<256;j++) 
    kh_destroy(hCount_t,hReads[j]);    

  /* --------------------------------------------------------------------- */

  //Now we can reparse the input file and output for each kmer the right location

  f=mustOpen(uniFile,"w");
  for(j=0;j<kh_end(hChr);j++)
    if(kh_exist(hChr,j)){
      rloc=kh_val(hChr,j);
      if(chromSizes[rloc]>0){
	fprintf(f,"@SQ\tSN:%s\tLN:%d\n",chromNames[rloc],chromSizes[rloc]);
      }
    }
  reParseFastqFileAndMap(inFastq ,hMapped,f,chromNames);
  carefulClose(&f);


  /* --------------------------------------------------------------------- */

  verbose(1,
	  "# Total mapped %d\n"
	  "# Total unique %d (%f %%)\n"
	  "# Total <=10 times %d (%f %%)\n"
	  "# Total Unmappable %d (%f %%)\n"
	  ,TotalMappedReads,
	  UniquelyMappedReads,UniquelyMappedReads/(float)TotalMappedReads*100,
	  LessThan10Reps,LessThan10Reps/(float)TotalMappedReads*100,
	  UnMappableReads,UnMappableReads/(float)(UnMappableReads+TotalMappedReads)*100);
  verbose(1,
	  "#Summary:%s\t%d\t%d\t%d\t%d\t%f\n"
	  ,inFastq
	  ,UniquelyMappedReads
	  ,TotalMappedReads-UniquelyMappedReads
	  ,UnMappableReads
	  ,UniquelyMappedReads+UnMappableReads
	  ,UnMappableReads/(float)(UnMappableReads+TotalMappedReads)*100);

  verbose(1,"# Total time: %f seconds = %f minutes\n",
	  (double)(clock() - t0)/CLOCKS_PER_SEC,
	  (double)(clock() - t0)/CLOCKS_PER_SEC/60);
}

/* ******************************************************************************************** */
/* ******************************************************************************************** */

int main(int argc, char *argv[])
/* Process command line. */
{
  optionInit(&argc, argv, options);

  rmdup = optionExists("rmdup");
  //  applySensitivity= optionExists("sensitivity");
  //  if(applySensitivity)
  //   sensitivityFolder=optionVal("sensitivity","");
  //  verbose(1,"%s \n",sensitivityFolder);

  ompNumThreads= optionInt("omp",ompNumThreads);

  omp_set_num_threads(ompNumThreads); //Adjust this to env variable.... or option

  if (argc != 4)
    usage();

  cutMapperOmpMapSam(argv[2],argv[1],argv[3]);

  return 0;
}

