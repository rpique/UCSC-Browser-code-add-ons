/* cutMapperFastQtoJfd - Maps reads in FastQ file using and indexed genome and makes output into jfd format. */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"

#include "cutMapperCore.h"


static char const rcsid[] = "$Id: newProg.c,v 1.30 2010/03/24 21:18:33 hiram Exp $";

void usage()
/* Explain usage and exit. */
{
errAbort(
  "cutMapperFastQtoJfd - Maps reads in FastQ file using and indexed genome and makes output into jfd format\n"
  "usage:\n"
  "   cutMapperFastQtoJfd indexFolder file.fastQ unique.ujdf repeats.rjfd \n"
  "options:\n"
  "   -xxx=XXX\n"
  );
}

static struct optionSpec options[] = {
   {NULL, 0},
};

/*
khash_t(hashChr_t) *chromSizes(char *chromFile, char ***chromNames, unsigned **chromSizes)
{
  FILE *f=mustOpen(chromFile,"r");
  int chrnum=1;
  char buff[256];
  unsigned size;
  khash_t(hashChr_t) *hChr= kh_init(hashChr_t); 
  
  khiter_t k;
  int hret;

  kvec_t(unsigned) kvChromSizes;
  kvec_t(char *)   kvChromNames;
  
  kv_init(kvChromSizes);
  kv_init(kvChromNames);

  kv_push(char *,kvChromNames,NULL);
  kv_push(unsigned,kvChromSizes,0);
    
  while(!feof(f)){
    chrnum++;
    assert(chrnum<256);
    fscanf(f,"%s\t%u\t", buff, &size);
    verbose(1,"#%s\t%u\n", buff, size);
    k = kh_put(hashChr_t, hChr , cloneString(buff), &hret);
    kh_val(hChr, k) = chrnum;
    kv_push(char *,kvChromNames,(char *)kh_key(hChr, k));
    kv_push(unsigned,kvChromSizes,size);
  }
  kv_trim(char *,kvChromNames);
  kv_trim(unsigned,kvChromSizes);
  
  *chromNames=kvChromNames.a;
  *chromSizes=kvChromSizes.a;

  return hChr;
}
*/

void reParseFastqFileAndMap(char *inFastq, khash_t(hashPos_t) **hReads, replist_t *replistp, char *uniFile, char *repFile, char **chromNames)
/* hashFastqFile - Parse file and map reads using the hash and repeat table */
{
  FILE *fU=mustOpen(uniFile,"w");
  FILE *fR=mustOpen(repFile,"w");
  
  struct lineFile *lf = lineFileOpen(inFastq, TRUE);
  char *seqName, *seq, *qName, *qual;
  int line = 0;
  //int numKmers = 0;
  int numUni = 0;
  int numRep = 0;
  int numSeq = 0;
  int numWithN=0;
  boolean startOfFile = TRUE;
  
  unsigned long int kmer;
  unsigned char Prefix=0;
  unsigned int Suffix;
  
  khiter_t khit;
  //int hret;
  int j;
  
  repvec_t *repvecp;
  unsigned int rloc;
  
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
    kmer=packKmer(seq,kmerSize);
    wantNewLine(lf, inFastq, ++line, &qName, "fastq sequence name (quality line)");
    wantNewLine(lf, inFastq, ++line, &qual, "quality line");
    if(kmer != -1){ // otherWise check how many N, and hash them in a different way... ? 4 times??
      Prefix=getPrefix(kmer); // ((ForKmer&kmermask)>>32);
      Suffix=getSuffix(kmer); // &0xFFFFFFFF;
      khit = kh_get(hashPos_t, hReads[Prefix] , Suffix);
      if ( kh_exist(hReads[Prefix], khit)){
	//In the hash... OUTPUT READ 	
	if(kh_value(hReads[Prefix], khit).chr==255){//Non uniquely mappable //255 Use 254 for debugging
	  //unPackKmer((((unsigned long int)j)<<32) + Suffix, 20, kmerString);
	  //Go to repeat table
	  rloc=kh_val(hReads[Prefix], khit).pos;
	  //fprintf(stderr,"#%d) %d,%x,%s,%d, %d \n",line,(int)Prefix,Suffix,seq,khit,rloc);
	  assert(rloc<=kv_size(*replistp));
	  repvecp=&kv_A(*replistp,rloc);
	  //fprintf(stdout,"repeats %d", (int) kv_size(*repvecp));
	  //Traverse repeats:
	  if(kv_size(*repvecp)<=10){
	    ++numRep;
	    fprintf(fR,"%s",seqName);	  
	    for(j=0;j<kv_size(*repvecp);j++){
	      fprintf(fR,"\t%s\t%d\t%c",
		      chromNames[kv_A(*repvecp,j).chr],
		      abs(kv_A(*repvecp,j).pos),
		      kv_A(*repvecp,j).pos<0 ? '-' : '+');
	    }
	    fprintf(fR,"\n");
	  }
	}
	else{ //Uniquely mappable
	  //unPackKmer((((unsigned long int)j)<<32) + Suffix, 20, kmerString);
	  numUni++;
	  //assert(kh_value(hReads[Prefix], khit).chr==14);
	  fprintf(fU,"%s\t%s\t%d\t%c\n",seqName,
		  chromNames[kh_value(hReads[Prefix], khit).chr],
		  abs(kh_value(hReads[Prefix], khit).pos),
		  kh_value(hReads[Prefix], khit).pos<0 ? '-' : '+');
	  
	  // fprintf(fU,"%d\n",kh_value(hReads[Prefix], khit).pos);

	}// else unmappable	  
	//fprintf(fU,"");
      }
      else{
	//UnMappable!!
      }
    }
    else{
      //unMappable with Ns
      numWithN++;
    }
    ++numSeq;
    freeMem(seqName);
  }
  
  verbose(1, "#\tUnique %d, Repeats <10 %d reads\n",numUni,numRep);
 
  verbose(1, "#\tprocessed %d lines  %d reads\n",line,numSeq);
  verbose(1, "#\tTotal number of reads with >0 N: %d\n",numWithN);
  
  //  verbose(1, "#\tTotal number of k-mers %d, %d \n",line,numKmers); //WHAT THE HELL IS GOING ON!!!
  //verbose(1, "#\tTotal number of rep k-mers %d \n",numRepKmers); 
  
  carefulClose(&fU);
  carefulClose(&fR);
  
  //  kh_destroy(hashPos_t, h);
  
  //close lf file...
  lineFileClose(&lf);
}

void cutMapperFastQtoJfd(char *inFastq, char *indexFolder, char *uniFile, char *repFile)
/* cutMapperFastQtoJfd - Maps reads in FastQ file using and indexed genome and makes output into jfd format. */
{
  //  FILE *f;
  char cbuff[strlen(indexFolder)+256];
  //char kmerString[21];
  
  int j,k;
  //int i;
  int TotalMappedReads=0;
  int UniquelyMappedReads=0;
  int UnMappableReads=0;
  int LessThan10Reps=0;
  khint_t khit,khit2,indexSize;  

  int hret;
  unsigned int rloc;
  
  clock_t t;
  clock_t t0=clock();
    
  khash_t(hCount_t) *hReads[256];
  
  khash_t(hashPos_t) *hIndex;  
  repvec_t *repvecp;  
  replist_t replist;
    
  khash_t(hashPos_t) *hMapped[256]; 
  replist_t replistMapped; 
  repvec_t *repvec2p;  

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
  //Initializing hash maps, for mapped reads!
  for(j=0;j<256;j++) 
    hMapped[j]=kh_init(hashPos_t);   
  //Initializing repeat maps, for mapped reads
  kv_init(replistMapped);
  
  //Initializing hash maps.
  for(j=0;j<256;j++) 
    hReads[j]=kh_init(hCount_t);   
  
  /* --------------------------------------------------------------------- */

  //Hashing input reads.
  t = clock();
  hashFastqFile(inFastq, hReads);
  verbose(1,"# Hashed fastq reads in %f seconds (%f total)\n",
	  (double)(clock() - t)/CLOCKS_PER_SEC, 
	  (double)(clock() - t0)/CLOCKS_PER_SEC); 
  
  /* --------------------------------------------------------------------- */
  
  //for(j=0;j<256;j++){ //256
  for(j=0;j<256;j++){ //256
    //Load repeats!!
    sprintf(cbuff,"%s/repeats%03d.bin",indexFolder,j);
    t = clock();
    repeatFileOpen(cbuff,&replist);
    verbose(1,"# Opened repeats in %f seconds (%f total)\n",
	    (double)(clock() - t)/CLOCKS_PER_SEC,
	    (double)(clock() - t0)/CLOCKS_PER_SEC);
    
    //Load Hash! 
    sprintf(cbuff,"%s/Hash%03d.bin",indexFolder,j);
    t = clock();
    hIndex = hashFileOpen(cbuff);
    verbose(1,"# Opened hash in %f seconds (%f total)\n",
	    (double)(clock() - t)/CLOCKS_PER_SEC,
	    (double)(clock() - t0)/CLOCKS_PER_SEC);
    
    //Now we can map the reads of prefix j, traversing the hReads[j] hash
    // Repeats can be resolved now, if we have a sensitivity index available, or later
    // making a second pass, or by storing the repeats in a separate hash. 
    // If we store the repeats, we could probably do it multiple times, if it fits in memory.
    
    /* --------------------------------------------------------------------- */

    t = clock();
    //Traverse the count hash
    indexSize=kh_size(hIndex);
    for (k = kh_begin(hReads[j]); k < kh_end(hReads[j]); ++k){
      if (kh_exist(hReads[j], k)){
	unsigned int Suffix=kh_key(hReads[j], k);
	unsigned int num=kh_val(hReads[j], k);
	//Look up the key, on the Suffix hash table hIndex
	khit = kh_get(hashPos_t, hIndex , Suffix);
	if(kh_exist(hIndex, khit)){ // I found it!
	  TotalMappedReads+=num;
	  if(kh_value(hIndex, khit).chr==255){//Non uniquely mappable
	    //print locations... if less than XXX locations...`
	    //unPackKmer((((unsigned long int)j)<<32) + Suffix, 20, kmerString);
	    //fprintf(stdout,"%s %d ",kmerString,num);
	    //Go to repeat table
	    repvecp=&kv_A(replist,kh_val(hIndex, khit).pos);
	    //fprintf(stdout,"repeats %d", (int) kv_size(*repvecp));
	    //Traverse repeats:
	    if(kv_size(*repvecp)<=10){
	      LessThan10Reps+=num;
	      khit2 = kh_put(hashPos_t, hMapped[j], Suffix, &hret);
	      assert(hret==1);
	      kv_push_empty(repvec_t,replistMapped);
	      rloc = kv_size(replistMapped)-1;
	      //fprintf(stderr,"#%d,%x,%d %d ",j,Suffix,khit2,rloc);
	      repvec2p = &kv_A(replistMapped,rloc);
	      kv_init(*repvec2p);
	      kv_copy(chrpos_t ,*repvec2p, *repvecp);
	      kh_val(hMapped[j],khit2).chr=255;
	      kh_val(hMapped[j],khit2).pos=rloc;	     	      
	      //fprintf(stderr,"= %d\n",kh_val(hMapped[j],khit2).pos);
	    }
	    /*
	    for(i=0;i<kv_size(*repvecp);i++){
	      fprintf(stdout," (%d,%d)",kv_A(*repvecp,i).chr,kv_A(*repvecp,i).pos);
	      }
	    fprintf(stdout,"\n");
	    */
	  }
	  else{ //Uniquely mappable
	    UniquelyMappedReads+=num;
	    /* 
	    unPackKmer((((unsigned long int)j)<<32) + Suffix, 20, kmerString);
	    fprintf(stdout,"%s %d at %d %d\n",kmerString,num,
		    kh_value(hIndex, khit).chr,
		    kh_value(hIndex, khit).pos);
	    
	    fprintf(stdout,"%d\n",kh_value(hIndex, khit).pos);
	    */
	    khit2 = kh_put(hashPos_t, hMapped[j], Suffix, &hret);
	    assert(hret==1);
	    kh_val(hMapped[j],khit2) = kh_val(hIndex, khit);	    
	  }
	}
	else{
	  UnMappableReads+=num;
	  //unPackKmer((((unsigned long int)j)<<32) + Suffix, 20, kmerString);
	  //fprintf(stdout,"%s %d Not mapped\n",kmerString,num);
	}
      }
    }/* index k for*/
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
    for (k = 0; k < kv_size(replist); ++k)
      kv_destroy(kv_A(replist,k));
    kv_destroy(replist);
    kh_destroy(hashPos_t,hIndex);    
  } /* Index j for */   
  verbose(1,"# Total time: %f seconds = %f minutes\n",
	  (double)(clock() - t0)/CLOCKS_PER_SEC,
	  (double)(clock() - t0)/CLOCKS_PER_SEC/60);

  
  for(j=0;j<256;j++) 
    kh_destroy(hCount_t,hReads[j]);    

  /* --------------------------------------------------------------------- */

  //Now we can reparse the input file and output for each kmer the right location
  
  reParseFastqFileAndMap(inFastq ,hMapped,&replistMapped,uniFile,repFile,chromNames);

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
  if (argc != 5)
    usage();
  cutMapperFastQtoJfd(argv[2],argv[1],argv[3],argv[4]);
  return 0;
}

