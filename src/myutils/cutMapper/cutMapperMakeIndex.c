/* cutMapperMakeIndex - Makes and index for an Specialized 20bp short read mapper for DNase-seq. */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"

//#include "twoBit.h"
#include "dnautil.h"
#include "dnaseq.h"
#include "dnaLoad.h"

#include <time.h>

#include "cutMapperCore.h"

static char const rcsid[] = "$Id: newProg.c,v 1.30 2010/03/24 21:18:33 hiram Exp $";

unsigned int idx = 0xFF; /* 4-mer index (0-255)*/
boolean hashOneSnpOpt = 0;
boolean hashAllSnpOpt = 0;


void usage()
/* Explain usage and exit. */
{
  errAbort(
	   "cutMapperMakeIndex - Makes and index for an Specialized 20bp short read mapper for DNase-seq\n"
	   "usage:\n"
	   "   cutMapperMakeIndex input.2bit idxfolder \n"
	   "options:\n"
	   "   -idx= index to generate (0-255)\n"
	   "   -one  hash the alt. genome, switching one SNP by one\n"
	   "   -all  hash the alt. genome, switching all SNPs once\n"
	   "   -xxx=XXX\n" 
	   );
}



static struct optionSpec options[] = {
  {"idx", OPTION_INT},
  {"one", OPTION_BOOLEAN},
  {"all", OPTION_BOOLEAN},  
  {NULL, 0},
};

/* ******************************************************************************************** */
/* ******************************************************************************************** */

void hashKmer(unsigned long int Kmer,
	      khash_t(hashPos_t) *h, 
	      replist_t *replistp,
	      unsigned char chr,
	      int pos){
  unsigned char Prefix;
  unsigned int Suffix;
  unsigned int rloc;
  
  khiter_t khit;
  int hret;
  
  repvec_t *repvecp;
  chrpos_t auxPos;
  
  auxPos.chr=chr;
  auxPos.pos=pos;
  
  /* Hash reverse strand */
  // function(Kmer,Strand,h,replistp,chr)
  Prefix=getPrefix(Kmer); // ((ForKmer&kmermask)>>32);
  Suffix=getSuffix(Kmer); // &0xFFFFFFFF;
  if(Prefix==idx){
    khit = kh_put(hashPos_t, h , Suffix, &hret);
    if (!hret){
      //Already in the hash
      if(kh_value(h, khit).chr==255){
	//kh_value(h, khit).pos++;
	repvecp=&kv_A(*replistp, kh_val(h, khit).pos);
	if(kv_size(*repvecp) < 128){
	  kv_push(chrpos_t,*repvecp,auxPos);
	}// if 128 or more we could change flag to 0 
	//  and just count repeats, and free memory. no longer need positions
      }
      else{
	//Check if the chr, and pos are different, otherwise it is not a problem!
	// This is important when scanning the reference genome...
	// kh_value(h , khit).chr = 255;
	//kh_value(h , khit).pos = 2;
	kv_push_empty(repvec_t,*replistp);
	rloc = kv_size(*replistp)-1;
	repvecp = &kv_A(*replistp,rloc);
	kv_init(*repvecp);
	//kv_init(kv_A(*replistp,rloc));
	kv_push(chrpos_t,*repvecp,kh_value(h, khit));
	kv_push(chrpos_t, *repvecp, auxPos);
	kh_value(h, khit).chr = 255;
	kh_value(h, khit).pos = rloc;	    
      }
    }
    else{
      kh_value(h , khit).pos = pos; //-(i-kmerSize+1);
      kh_value(h , khit).chr = chr;
    }	
  } 
}

/* ******************************************************************************************** */
/* ******************************************************************************************** */

//This function modifies the previous one to only add repeats.
void hashKmerCheckRepeatPos(unsigned long int Kmer,
	      khash_t(hashPos_t) *h, 
	      replist_t *replistp,
	      unsigned char chr,
	      int pos){
  unsigned char Prefix;
  unsigned int Suffix;
  unsigned int rloc;
  
  khiter_t khit;
  int hret;
  int k;
  
  repvec_t *repvecp;
  chrpos_t auxPos;
  
  auxPos.chr=chr;
  auxPos.pos=pos;
  
  // function(Kmer,Strand,h,replistp,chr)
  Prefix=getPrefix(Kmer); // ((ForKmer&kmermask)>>32);
  Suffix=getSuffix(Kmer); // &0xFFFFFFFF;
  if(Prefix==idx){
    khit = kh_put(hashPos_t, h , Suffix, &hret);
    if (!hret){
      //Already in the hash
      if(kh_value(h, khit).chr==255){
	//kh_value(h, khit).pos++;
	repvecp=&kv_A(*replistp, kh_val(h, khit).pos);
	if(kv_size(*repvecp) < 128){
	  //Need to check if in the list.
	  int isNotIn=1;
	  for(k=0; k<kv_size(*repvecp); k++)
	    if((kv_A(*repvecp,k).chr==auxPos.chr) && 
	       (kv_A(*repvecp,k).pos==auxPos.pos)){
	      isNotIn=0;
	      break;
	    } 
	  if(isNotIn) // It is entering an repeat from the alternate genome.
	    kv_push(chrpos_t,*repvecp,auxPos);
	}// if 128 or more we could change flag to 0 
	//  and just count repeats, and free memory. no longer need positions
	// But this just makes handling the vectors more complicated....
      }
      else{
	//Check if the chr, and pos are different, otherwise it is not a problem!
	// This is important when scanning the non-reference genome...
	//  = 255;
	//kh_value(h , khit).pos = 2;
	if((kh_value(h , khit).chr==auxPos.chr) && 
	   (kh_value(h , khit).pos==auxPos.pos)){
	  //Do nothing... 
	}else{
	  kv_push_empty(repvec_t,*replistp);
	  rloc = kv_size(*replistp)-1;
	  repvecp = &kv_A(*replistp,rloc);
	  kv_init(*repvecp);
	  //kv_init(kv_A(*replistp,rloc));
	  kv_push(chrpos_t,*repvecp,kh_value(h, khit));
	  kv_push(chrpos_t, *repvecp, auxPos);
	  kh_value(h, khit).chr = 255;
	  kh_value(h, khit).pos = rloc;	 
	}   
      }
    }
    else{
      kh_value(h , khit).pos = pos; //-(i-kmerSize+1);
      kh_value(h , khit).chr = chr;
    }	
  } 
}

/* ******************************************************************************************** */
/* ******************************************************************************************** */

int scankmer(struct dnaSeq *seq,khash_t(hashPos_t) *h, replist_t *replistp,unsigned char chr){
  int i,Start,End;

  //DNA buff[kmerSize+1];
  //unsigned long int ForKmer2=0,RevKmer2=0;
  unsigned long int ForKmer=0,RevKmer=0;
  
  unsigned int mask=0xFFFFFFFF;  //0x000FFFFF
  char Base;
  unsigned long int b;

  Start=0;
  End=seq->size;//100;//seq->size-kmerSize;
  if((End-Start+1) < kmerSize){
    return 0; //Region too small to fit the kmer..
  }
  
  //#pragma omp parallel for private(nonATGCbase,j,score) reduction(+:count)
  for(i=Start;i<=End;i++){
    Base=seq->dna[i];

    Base=Base&0xDF; // 11011111 sets letter to uppercase
    mask<<=1;       // This bitMask monitors how many letters are valid 
    if(Base=='A') {
      b=0;
    }else if(Base=='C'){
      b=1;
    }else if(Base=='G'){
      b=2;
    }else if(Base=='T'){
      b=3;
    }else{
      b=0;
      mask|=1;
    }
    ForKmer = (ForKmer<<2) + b;
    RevKmer = (RevKmer>>2) + ((((~b)&0x03))<<38);
    
    if(!(mask&0x00FFFFF)){
      //Debugging stuff. 
      //unPackKmer(ForKmer&KmerMask,kmerSize,buff);
      //verbose(1,"#1)%d %s\t%lx\t",i-kmerSize+1,buff,ForKmer&KmerMask);      
      //unPackKmer(RevKmer,kmerSize&KmerMask,buff);
      //verbose(1,"%s\t%lx\n",buff,RevKmer&KmerMask);
      /* Hash forward strand */
      hashKmerCheckRepeatPos(ForKmer,h,replistp,chr,(i + 1 - kmerSize + 1));  //add 1 for 1-based coords. 
      /* Hash reverse strand */
      hashKmerCheckRepeatPos(RevKmer,h,replistp,chr,-(i + 1 - kmerSize + 1));
    }
  } 
  return 1;
}

/* I CAN'T USE 0 BASED COORDINATES IF I USE THE SIGN FOR STRAND. +0 and -0 are the same thing! */

/* ******************************************************************************************** */
/* ******************************************************************************************** */

//    scankmersnp(seq,h,&replist,kh_val(hChr , k), kv_A(snplist,k));   
int scankmersnp(struct dnaSeq *seq,khash_t(hashPos_t) *h, replist_t *replistp,unsigned char chr,snpvec_t *snpvecp){
  int i,j,Start,End;

  //DNA buff[kmerSize+1];
  //unsigned long int ForKmer2=0,RevKmer2=0;
  unsigned long int ForKmer=0,RevKmer=0;
  
  unsigned int mask=0xFFFFFFFF;  //0x000FFFFF
  char Base;
  unsigned long int b;
  
  int lastpos=-1000; // to check snps ordered... 

  snp_t snpaux;

  // handle if lastpos < pos-19 ?? double switch??
  
  for(j=0; j< kv_size(*snpvecp); j++){
    snpaux=kv_A(*snpvecp,j);
    //  assert(lastpos<snpaux.pos); //Asserting in order...  ACTUALLY NOT NECESSARY 
    //  verbose(2,"%d,%d, %c %c\n",(int)chr,snpaux.pos,(seq->dna[snpaux.pos]),(snpaux.ref));

    assert(((snpaux.ref) & 0xDF) == ((seq->dna[snpaux.pos]) & 0xDF)); //Asserting Reference matches SNP
    assert(snpaux.pos < seq->size); //SNP out of limits

    // handle if lastpos < pos-19 ?? double switch??
    Start = snpaux.pos - kmerSize + 1;
    if(Start<0)
      continue;
    End = snpaux.pos + kmerSize; //+seq->size;//100;//seq->size-kmerSize;
    if(End > (seq->size - 1))
      continue;
    if((End-Start+1) < kmerSize)
      continue; //Region too small to fit the kmer..
    
    seq->dna[snpaux.pos]=snpaux.alt;

    //#pragma omp parallel for private(nonATGCbase,j,score) reduction(+:count)
    for(i=Start;i<=End;i++){
      Base=seq->dna[i];
      
      Base=Base&0xDF; // 11011111 sets letter to uppercase
      mask<<=1;       // This bitMask monitors how many letters are valid 
      if(Base=='A') {
	b=0;
      }else if(Base=='C'){
	b=1;
      }else if(Base=='G'){
	b=2;
      }else if(Base=='T'){
	b=3;
      }else{
	b=0;
	mask|=1;
      }
      ForKmer = (ForKmer<<2) + b;
      RevKmer = (RevKmer>>2) + ((((~b)&0x03))<<38);
    
      if(!(mask&0x00FFFFF)){
	//Debugging stuff. 
	//unPackKmer(ForKmer&KmerMask,kmerSize,buff);
	//verbose(1,"#1)%d %s\t%lx\t",i-kmerSize+1,buff,ForKmer&KmerMask);      
	//unPackKmer(RevKmer,kmerSize&KmerMask,buff);
	//verbose(1,"%s\t%lx\n",buff,RevKmer&KmerMask);
	/* Hash forward strand */
	hashKmerCheckRepeatPos(ForKmer,h,replistp,chr,(i + 1 - kmerSize + 1));
	/* Hash reverse strand */
	hashKmerCheckRepeatPos(RevKmer,h,replistp,chr,-(i + 1 - kmerSize + 1));
      }
    } 

    // handle if lastpos < pos-19 ?? double switch??
    
    seq->dna[snpaux.pos]=snpaux.ref; // This restores the to the previous state... 
    lastpos=snpaux.pos;
    
  }
    return 1;
}

/* ******************************************************************************************** */

int scankmerindel(struct dnaSeq *seq,khash_t(hashPos_t) *h, replist_t *replistp,unsigned char chr,indelvec_t *indelvecp){
  int i,j,k,Start,End,pos;

  //DNA buff[kmerSize+1];
  //unsigned long int ForKmer2=0,RevKmer2=0;
  unsigned long int ForKmer=0,RevKmer=0;
  
  unsigned int mask=0xFFFFFFFF;  //0x000FFFFF
  char Base;

  char buff[1000];
  unsigned long int b;
  
  indel_t aux;

  // handle if lastpos < pos-19 ?? double switch??
  
  for(j=0; j< kv_size(*indelvecp); j++){
    aux=kv_A(*indelvecp,j);
    //verbose(2,"%d,%d, %c %c\n",(int)chr,snpaux.pos,(seq->dna[snpaux.pos]),(snpaux.ref));

    // Assert Indel matches de reference....
    if(strncasecmp(&seq->dna[aux.pos],aux.ref,aux.reflen)!=0){
      verbose(1,"## %d) %d,%s\n",(int)chr,aux.pos,(aux.ref));
      errAbort("Indel not matching the reference\n");
    }

    //Check or skip
    assert(aux.pos < seq->size); //SNP out of limits
    assert(aux.pos >= 0); //SNP out of limits


    // handle if lastpos < pos-19 ?? double switch??
    Start = aux.pos - kmerSize + 1;
    if(Start<0)
      continue;
    End = aux.pos + kmerSize + aux.reflen -1; // I need to add the deletion size... //+seq->size;//100;//seq->size-kmerSize;
    if(End > (seq->size - 1))
      continue;
    if((End-Start+1) < kmerSize)
      continue; //Region too small to fit the kmer..
    
    //verbose(1,"## %d) %d,%s --> %s\n",(int)chr,aux.pos,(aux.ref),aux.alt);
    /*
    for(k=Start;k<End;k++)
      verbose(1,"%c",seq->dna[k]);
    verbose(1,"\n");
    */

    // Extract DNA sequence from reference...
    i=0;
    for(k=Start;k<aux.pos;k++) // First section...
      buff[i++]=seq->dna[k];
    //buff[i++]='_';
    for(k=0;k<aux.altlen;k++)
      buff[i++]=aux.alt[k];
    //buff[i++]='_';
    for(k=aux.pos+aux.reflen;k<=End;k++)
      buff[i++]=seq->dna[k];
    k=i;

    /*
    for(i=0;i<k;i++)
      verbose(1,"%c",buff[i]);
    verbose(1,"\n");
    */

    // Edit the DNA sequence in the buffer...
    //  seq->dna[snpaux.pos]=snpaux.alt;

    pos=Start+0+1-kmerSize+1;
    //#pragma omp parallel for private(nonATGCbase,j,score) reduction(+:count)
    for(i=0;i<k;i++){
      Base=buff[i];
      
      Base=Base&0xDF; // 11011111 sets letter to uppercase
      mask<<=1;       // This bitMask monitors how many letters are valid 
      if(Base=='A') {
	b=0;
      }else if(Base=='C'){
	b=1;
      }else if(Base=='G'){
	b=2;
      }else if(Base=='T'){
	b=3;
      }else{
	b=0;
	mask|=1;
      }
      ForKmer = (ForKmer<<2) + b;
      RevKmer = (RevKmer>>2) + ((((~b)&0x03))<<38);
    
      if(!(mask&0x00FFFFF)){
	//Debugging stuff. 
	//unPackKmer(ForKmer&KmerMask,kmerSize,buff);
	//verbose(1,"#1)%d %s\t%lx\t",i-kmerSize+1,buff,ForKmer&KmerMask);      
	//unPackKmer(RevKmer,kmerSize&KmerMask,buff);
	//verbose(1,"%s\t%lx\n",buff,RevKmer&KmerMask);
	if(i<kmerSize)
	  pos=Start+i+1-kmerSize+1;
	if((aux.reflen>1)&&(i==kmerSize))
	  pos=pos+(aux.reflen-1);
	if((aux.altlen>1)&&(i>=kmerSize)&&(i<kmerSize+aux.altlen-1))
	  pos=pos;
	if(i>=kmerSize+aux.altlen-1)
	  pos++;
	//verbose(1,"%d,",pos-(Start+1-kmerSize+1));	
	/* Hash forward strand */
	hashKmerCheckRepeatPos(ForKmer,h,replistp,chr,pos);
	/* Hash reverse strand */
	hashKmerCheckRepeatPos(RevKmer,h,replistp,chr,-pos);
      }
    } 
    //verbose(1,"\n");
  }
    return 1;
}

/* ******************************************************************************************** */
/* ******************************************************************************************** */

//   scankmersnp(seq,h,&replist,kh_val(hChr , k), kv_A(snplist,k));   
//   This hashes the complete alternate genome. (If two SNPs < 20bp  it will hash the kmers with the two alternate)
//   This function also changes the sequence!, so we have to be careful. that we use it last!
int scankmersnp2(struct dnaSeq *seq,khash_t(hashPos_t) *h, replist_t *replistp,unsigned char chr,snpvec_t *snpvecp){
  int i,j,Start,End;
  
  //DNA buff[kmerSize+1];
  //unsigned long int ForKmer2=0,RevKmer2=0;
  unsigned long int ForKmer=0,RevKmer=0;
  
  unsigned int mask=0xFFFFFFFF;  //0x000FFFFF
  char Base;
  unsigned long int b;
  
  int lastpos=-1000; // to check snps ordered... // we don't need to check in this case 

  snp_t snpaux;
  
  // Changing all the bases to the alternate version
  for(j=0; j< kv_size(*snpvecp); j++){
    snpaux=kv_A(*snpvecp,j);
    if(lastpos==snpaux.pos)
      continue;  // 
    //assert(lastpos<snpaux.pos); //Asserting in order...
    //verbose(2,"%d,%d, %c %c\n",(int)chr,snpaux.pos,(seq->dna[snpaux.pos]),(snpaux.ref));
    assert(((snpaux.ref) & 0xDF) == ((seq->dna[snpaux.pos]) & 0xDF)); //Asserting Reference matches SNP
    assert(snpaux.pos < seq->size);    //SNP out of limits
    seq->dna[snpaux.pos]=snpaux.alt;   // Changing to the alternate version.
    lastpos=snpaux.pos;
  }

  // Hashing the new alternate version
  for(j=0; j< kv_size(*snpvecp); j++){
    snpaux=kv_A(*snpvecp,j);

    Start = max(snpaux.pos - kmerSize + 1, 0);
    End = min(snpaux.pos + kmerSize, seq->size); //+seq->size;//100;//seq->size-kmerSize;
    if((End-Start+1) < kmerSize)
      continue; //Region too small to fit the kmer..
    
    //#pragma omp parallel for private(nonATGCbase,j,score) reduction(+:count)
    for(i=Start;i<=End;i++){
      Base=seq->dna[i];      
      Base=Base&0xDF; // 11011111 sets letter to uppercase
      mask<<=1;       // This bitMask monitors how many letters are valid 
      if(Base=='A') {
	b=0;
      }else if(Base=='C'){
	b=1;
      }else if(Base=='G'){
	b=2;
      }else if(Base=='T'){
	b=3;
      }else{
	b=0;
	mask|=1;
      }
      ForKmer = (ForKmer<<2) + b;
      RevKmer = (RevKmer>>2) + ((((~b)&0x03))<<38);
    
      if(!(mask&0x00FFFFF)){
	//Debugging stuff. 
	//unPackKmer(ForKmer&KmerMask,kmerSize,buff);
	//verbose(1,"#1)%d %s\t%lx\t",i-kmerSize+1,buff,ForKmer&KmerMask);      
	//unPackKmer(RevKmer,kmerSize&KmerMask,buff);
	//verbose(1,"%s\t%lx\n",buff,RevKmer&KmerMask);
	/* Hash forward strand */
	hashKmerCheckRepeatPos(ForKmer,h,replistp,chr,(i + 1 - kmerSize + 1));
	/* Hash reverse strand */
	hashKmerCheckRepeatPos(RevKmer,h,replistp,chr,-(i + 1 - kmerSize + 1));
      }
    }     
  }

  //Restore reference genome, here???

    return 1;
}


/* ******************************************************************************************** */
/* ******************************************************************************************** */

/* ******************************************************************************************** */
/* ******************************************************************************************** */

void cutMapperMakeIndex(char *fileSeq, char *fileSNPs,char *outFolder)
/* cutMapperMakeIndex - Makes and index for an Specialized 20bp short read mapper for DNase-seq. */
{
  //FILE *fSnp = mustOpen(fileSNPs, "r");
  FILE *f; //For dumping index
  char cbuff[strlen(outFolder)+256];
  struct dnaLoad *dl = dnaLoadOpen(fileSeq); // change to 2bit file alone?. 
  //  struct twoBitFile *tbf;
  //  struct twoBitSpec *tbs;
  
  struct dnaSeq *seq; 
  int j;
  khiter_t k;
  size_t vsize;
  size_t cumsize=0;
  
  int hret;
  char *chrNames[256];
  clock_t t;
  
  //khash_t(32) *h = kh_init(32);  
  //DNA buff[kmerSize+1];
  unsigned char chrnum=1;

  khash_t(hashPos_t) *h = kh_init(hashPos_t);  

  khash_t(hashChr_t) *hChr = kh_init(hashChr_t);  

  repvec_t *repvecp;  
  replist_t replist;
  kv_init(replist);  
  
  //snpvec_t *snpvecp;  
  indellist_t indellist;
  snplist_t snplist;
  // snp_t snpaux;

  //Opening input sequence, now in 2bit format!
  // Why?, because I will like to go to each chromosome and hash with SNPs...
  // also the entire sequence is in memory I think and make this easier. 
  // It will also ensure check file consistency.
  //  tbs = twoBitSpecNew(fileSeq);
  //  tbf = twoBitOpen(tbs->fileName);


  // Create output folder if not existing...
  sprintf(cbuff,"mkdir -p %s",outFolder);
  assert(system(cbuff)==0);
  
  //  processSeqsFromBinFile(tbf, binName, outFile);  

  /* --------------------------------------------------------------------- */

  while ((seq = dnaLoadNext(dl)) != NULL)
    {
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
      scankmer(seq,h,&replist,kh_val(hChr , k));   
      dnaSeqFree(&seq);
      //if(chrnum>2)break;
    } 
  verbose(1,"Finished hashing reference genome\n");
  verbose(1,"khash: n_buckets=%d, size=%d, n_occupied=%d, upper_bound=%d \n",kh_n_buckets(h),kh_size(h),((h)->n_occupied),((h)->upper_bound));  
  dnaLoadClose(&dl);

  /* --------------------------------------------------------------------- */
  
  if(hashOneSnpOpt || hashAllSnpOpt){
    
    loadSnpFile(fileSNPs, hChr, &snplist, &indellist);
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
  
    verbose(1,"Hashing alterative kmers using a SNP file \n");  
    dl = dnaLoadOpen(fileSeq);
    while ((seq = dnaLoadNext(dl)) != NULL){
      k = kh_get(hashChr_t, hChr , seq->name);
      j = kh_val(hChr, k);
      if(hashOneSnpOpt){
	verbose(1,"#%d) Re-hashing %d SNPs in %s (1 single alt.), idx=%d\n",
		j, (int)kv_size(kv_A(snplist,j)), chrNames[j], idx);
	scankmersnp(seq, h, &replist, j, &kv_A(snplist,j));   
	verbose(1,"#%d) Re-hashing %d INDELs in %s (1 single alt.), idx=%d\n",
		j, (int)kv_size(kv_A(indellist,j)), chrNames[j], idx);
	scankmerindel(seq, h, &replist, j, &kv_A(indellist,j));
      }
      if(hashAllSnpOpt){
	verbose(1,"#%d) Re-hashing %d SNPs in %s (complete alt.), idx=%d\n", 
		j, (int)kv_size(kv_A(snplist,j)), chrNames[j], idx);
	scankmersnp2(seq, h, &replist, j, &kv_A(snplist,j));   
      }
      dnaSeqFree(&seq);
      //if(chrnum>5)break;
    }
    verbose(1,"Finished hashing alternative genomes\n");
    verbose(2,"khash: n_buckets=%d, size=%d, n_occupied=%d, upper_bound=%d \n",kh_n_buckets(h),kh_size(h),((h)->n_occupied),((h)->upper_bound));

    dnaLoadClose(&dl);
  
    // Destroy SNPlist
    //for(k=1;k<chrnum;k++){
    //  Destroy SnpLists.....
    //  verbose(2,"#%d) %s %d Snps\n",k,chrNames[k],(int)kv_size(kv_A(snplist,k)));
    //}
 



  }
  
  /* --------------------------------------------------------------------- */  

  //Dumping repeats into file
  sprintf(cbuff,"%s/repeats%03d.bin",outFolder,idx);
  verbose(1,"Dumping repeats into file %s\n",cbuff);
  verbose(1," size int= %ld  chrpos_t= %ld\n",sizeof(int),sizeof(chrpos_t));
  f=mustOpen(cbuff,"wb");
  verbose(1," Number of repeated kmers %d\n",(unsigned int)kv_size(replist));
  fwrite(&kv_size(replist),sizeof(size_t),1,f);
  for (k = 0; k < kv_size(replist); ++k){
    repvecp=&kv_A(replist,k);
    vsize=kv_size(*repvecp);
    fwrite(&vsize,sizeof(size_t),1,f);
    fwrite(repvecp->a,sizeof(chrpos_t),vsize,f);
    kv_destroy(*repvecp);
    cumsize+=vsize;
  }
  carefulClose(&f);
  kv_destroy(replist);
  verbose(1," Repeats written, %d repeats. \n", (unsigned)cumsize);

  //Dump hash table into file
  sprintf(cbuff,"%s/Hash%03d.bin",outFolder,idx);
  verbose(1,"Dumping hash table into file %s\n",cbuff);
  verbose(1," size int= %ld  chrpos_t= %ld\n",sizeof(int),sizeof(chrpos_t));
  f=mustOpen(cbuff,"wb");
  fwrite(&kh_size(h),sizeof(khint_t),1,f);
  for (k = kh_begin(h); k < kh_end(h); ++k)
    if (kh_exist(h, k)){
      fwrite(&kh_key(h, k),sizeof(uint32_t),1,f);
      fwrite(&(kh_value(h, k)),sizeof(chrpos_t),1,f);
      //fwrite(&(kh_value(h, k).chr),sizeof(char),1,f);
      //fwrite(&(kh_value(h, k).pos),sizeof(int),1,f);
    }
  carefulClose(&f);
  kh_destroy(hashPos_t, h);

  /* --------------------------------------------------------------------- */

  //Reload repeats!!
  sprintf(cbuff,"%s/repeats%03d.bin",outFolder,idx);
  t = clock();
  repeatFileOpen(cbuff,&replist);
  verbose(1," Loaded in %f seconds\n",(double)(clock() - t)/CLOCKS_PER_SEC);
  
  //Rehash! 
  sprintf(cbuff,"%s/Hash%03d.bin",outFolder,idx);
  t = clock();
  h = hashFileOpen(cbuff);
  verbose(1," Loaded in %f seconds\n",(double)(clock() - t)/CLOCKS_PER_SEC);
  verbose(1,"khash: n_buckets=%d, size=%d, n_occupied=%d, upper_bound=%d \n",
	  kh_n_buckets(h),kh_size(h),((h)->n_occupied),((h)->upper_bound));

  /* --------------------------------------------------------------------- */

  
  //Debug hash by traverse. again to  see if files match! 
  /*
  f=mustOpen("test1.txt","w");
  for (k = kh_begin(h); k < kh_end(h); ++k)
    if (kh_exist(h, k)){
      if (kh_value(h, k).chr==255){
	unPackKmer((((unsigned long int) idx )<<32) + kh_key(h, k), 20, buff);
	repvecp=&kv_A(replist,kh_val(h, k).pos);
	fprintf(f,"%x\t%s\t%d\tP:%d\tR:%d\n", kh_key(h, k), buff , 255, kh_val(h, k).pos, (int) kv_size(*repvecp));
	//Traverse repeats:
	for(j=0;j<kv_size(*repvecp);j++){
	  fprintf(f,"  %d\t%s\t%d\n",j,chrNames[kv_A(*repvecp,j).chr],kv_A(*repvecp,j).pos);
	  //fprintf("  %d\n",j);
	}
      }else{
	unPackKmer((((unsigned long int) idx )<<32) + kh_key(h, k), 20, buff);
	fprintf(f,"%x\t%s\t%s\t%d\n", kh_key(h, k), buff ,chrNames[kh_val(h, k).chr], kh_val(h, k).pos );	
      }
    }
  carefulClose(&f);
  */

  /* --------------------------------------------------------------------- */
  
  //Destroy repeats
  kh_destroy(hashPos_t, h);
  kh_destroy(hashChr_t, hChr);
}


int main(int argc, char *argv[])
/* Process command line. */
{
  optionInit(&argc, argv, options);
  if (argc != 4)
    usage();
  hashOneSnpOpt = optionExists("one");
  hashAllSnpOpt = optionExists("all");
  idx=optionInt("idx", -1);

  if(idx<0 || idx>255)
    errAbort("Please, specify a valid index idx=0,...,255\n");

  dnaUtilOpen();
  cutMapperMakeIndex(argv[1],argv[2],argv[3]);
  return 0;
}

  // Discarded stuff in this function  (TO CLEAN UP) 


// Old SNP file opening!
    /*
    verbose(1,"# Reading SNP file\n");
    kv_init_size(snpvec_t,snplist,chrnum); // Make it into a function 
    // readSNPfile(snplist,indellist,snpFile)
    for(k=0;k<chrnum;k++)
      kv_init(kv_A(snplist,k));
    j=0;
    while(!feof(fSnp)){
      if(fscanf(fSnp,"%s\t%d\t%c\t%c%*[^\n]\n",cbuff,&snpaux.pos,&snpaux.ref,&snpaux.alt)!=4 || feof(fSnp))
	break;
      ++j;
      verbose(4,"#%d)%s\t%d\t%c\t%c\n",j,cbuff,snpaux.pos,snpaux.ref,snpaux.alt);
      k = kh_get(hashChr_t, hChr , cbuff);
      assert(kh_exist(hChr,k));
      k = kh_val(hChr, k);
      kv_push(snp_t,kv_A(snplist,k),snpaux);
    }
    verbose(2,"#%d)%s\t%d\t%c\t%c\n",j,cbuff,snpaux.pos,snpaux.ref,snpaux.alt);
    verbose(1,"# There are %d snps\n",j);
    //  for(k=1;k<chrnum;k++)
    //   verbose(2,"#%d) %s %d Snps\n",k,chrNames[k],(int)kv_size(kv_A(snplist,k)));
    carefulClose(&fSnp);
    */



/*
  //For debug // kv_size(replist)
  f=mustOpen("testA.txt","w");
  for (k = 0; k < kv_size(replist); ++k){ 
    repvecp=&kv_A(replist,k);
    vsize=kv_size(*repvecp);
    fprintf(f," Unit: %d Size: %d\n",(unsigned)k, (unsigned)vsize);
    for(j = 1; j < vsize; j++)
      fprintf(f,"  %d\t%s\t%d\n",j,chrNames[kv_A(*repvecp,j).chr],kv_A(*repvecp,j).pos);
  }
  carefulClose(&f);

*/
  /*
  //Debug hash by traverse. 
  f=mustOpen("test1.txt","w");
  for (k = kh_begin(h); k < kh_end(h); ++k)
    if (kh_exist(h, k)){
      if (kh_value(h, k).chr==255){
	unPackKmer((((unsigned long int) idx )<<32) + kh_key(h, k), 20, buff);
	repvecp=&kv_A(replist,kh_val(h, k).pos);
	fprintf(f,"%x\t%s\t%d\tP:%d\tR:%d\n", kh_key(h, k), buff , 255, kh_val(h, k).pos, (int) kv_size(*repvecp));
	//Traverse repeats:
	for(j=0;j<kv_size(*repvecp);j++){
	  fprintf(f,"  %d\t%s\t%d\n",j,chrNames[kv_A(*repvecp,j).chr],kv_A(*repvecp,j).pos);
	  //fprintf("  %d\n",j);
	}
      }else{
	unPackKmer((((unsigned long int) idx )<<32) + kh_key(h, k), 20, buff);
	fprintf(f,"%x\t%s\t%s\t%d\n", kh_key(h, k), buff ,chrNames[kh_val(h, k).chr], kh_val(h, k).pos );	
      }
    }
  carefulClose(&f);
  */

  /*
  //Debug hash by traverse. 
  for (k = kh_begin(h); k < kh_end(h); ++k)
  if (kh_exist(h, k)) 
  if (kh_value(h, k)>100){
    unPackKmer((((unsigned long int) idx )<<32) + kh_key(h, k), 20, buff);
	fprintf(stdout,"%x\t%s\t%d\n", kh_key(h, k), buff , kh_value(h, k));
      }
  */

      /*
      // function(Kmer,Strand,h,replistp,chr)
      Prefix=getPrefix(RevKmer); // ((ForKmer&kmermask)>>32);
      Suffix=getSuffix(RevKmer); //&0xFFFFFFFF;
      if(Prefix==idx){
	khit = kh_put(hashPos_t, h , Suffix, &hret);
	if (!hret){
	  //Already in the hash
	  if(kh_value(h, khit).chr==255){
	    //kh_value(h, khit).pos++;
	    repvecp=&kv_A(*replistp, kh_val(h, khit).pos);
	    if(kv_size(*repvecp) < 128){
	      auxPos.pos = -(i - kmerSize + 1);
	      kv_push(chrpos_t,*repvecp,auxPos);
	    }
	  }
	  else{
	    //Check if the chr, and pos are different, otherwise it is not a problem!
	    // This is important when scanning the reference genome...
	    // kh_value(h , khit).chr = 255;
	    //kh_value(h , khit).pos = 2;
	    kv_push_empty(repvec_t,*replistp);
	    rloc=kv_size(*replistp)-1;
	    repvecp=&kv_A(*replistp,rloc);
	    kv_init(*repvecp);
	    //kv_init(kv_A(*replistp,rloc));
	    kv_push(chrpos_t,*repvecp,kh_value(h, khit));
	    auxPos.pos = -(i-kmerSize+1);
	    kv_push(chrpos_t, *repvecp,auxPos);
	    kh_value(h, khit).chr = 255;
	    kh_value(h, khit).pos = rloc;	    
	  }
	}
	else{
	  kh_value(h , khit).pos = -(i-kmerSize+1);
	  kh_value(h , khit).chr = chr;
	}	
      }

      Prefix=getPrefix(ForKmer); // ((ForKmer&kmermask)>>32);
      Suffix=getSuffix(ForKmer); //&0xFFFFFFFF;
      if(Prefix==idx){
	//unPackKmer(Prefix,4,buff);
	//verbose(1,"#3a)%d %s\t%x\t",i,buff,(unsigned int)Prefix);      
	// unPackKmer(Suffix,16,buff);
	//verbose(1,"%s\t%x\n",buff,Suffix);	
	khit = kh_put(hashPos_t, h , Suffix, &hret);
	if (!hret){
	  //Already in the hash
	  if(kh_value(h , khit).chr==255){ //Already in the repeat list
	    //kh_value(h , khit).pos++; // For the moment do nothing 
	    // maybe limit to 1024 max_hits...
	    repvecp = &kv_A(*replistp, kh_val(h, khit).pos);
	    if(kv_size(*repvecp) < 128){
	      auxPos.pos = i - kmerSize + 1;
	      kv_push(chrpos_t, *repvecp, auxPos);
	    }
	  }
	  else{ //Entering to the repeat list
	    // typedef kvec_t(chrpos_t)repvec_t;  //Repeats of a kmer
	    // typedef kvec_t(repvec_t)replist_t; //List of all kmer repeat vectors
	    //repvecp=kv_pushp(repvec_t,*replistp);
	    kv_push_empty(repvec_t,*replistp);
	    rloc = kv_size(*replistp) - 1;
	    repvecp = &kv_A(*replistp,rloc);
	    kv_init(*repvecp);
	    //kv_init(kv_A(*replistp,rloc));
	    kv_push(chrpos_t,*repvecp,kh_value(h, khit));
	    auxPos.pos = i - kmerSize + 1;
	    kv_push(chrpos_t,*repvecp,auxPos);
	    kh_value(h, khit).chr = 255;
	    kh_value(h, khit).pos = rloc;	    
	  }
	}
	else{
	  kh_value(h , khit).pos = i-kmerSize+1;
	  kh_value(h , khit).chr = chr;
	}	
      }
      */

    /*
    if(i>=kmerSize-1)
      ForKmer2 = packKmer(seq->dna+i-kmerSize+1,kmerSize);
    else
      ForKmer2 = -1;
    if(ForKmer2 != -1){
      strncpy(buff,seq->dna+i-kmerSize+1,kmerSize);
      verbose(1,"#2)%d %s\t%lx\t",i-kmerSize+1,buff,ForKmer2);      
      reverseComplement(buff,kmerSize);
      RevKmer2 = packKmer(buff,kmerSize);
      verbose(1,"%s\t%lx\n",buff,RevKmer2);
    }
    */  

  /* Debugging stuff. 
	 unPackKmer(ForKmer,kmerSize,buff);
	 verbose(1,"#2)%d %s\t%lx\t",i,buff,ForKmer);      
	 unPackKmer(RevKmer,kmerSize,buff);
	 verbose(1,"%s\t%lx\n",buff,RevKmer);
      */
    /* 
    ForKmer2 = packKmer(seq->dna+i,kmerSize);
    if(ForKmer != -1){
      //buff=cloneStringZ(seq->dna+i,kmerSize);
      strncpy(buff,seq->dna+i,kmerSize);
      //verbose(1,"#1) %s\t%lx\t",buff,ForKmer);      
      reverseComplement(buff,kmerSize);
      RevKmer = packKmer(buff,kmerSize);
      //verbose(1,"%s\t%lx\n",buff,RevKmer);
      */
      /* Hash forward strand */
    /*
      Prefix=(ForKmer>>32)&0xFF;
      Suffix=(ForKmer)&0xFFFFFFFF;
      if(Prefix==idx){
	//unPackKmer(Prefix,4,buff);
	//verbose(1,"#3a)%d %s\t%x\t",i,buff,(unsigned int)Prefix);      
	unPackKmer(Suffix,16,buff);
	//verbose(1,"%s\t%x\n",buff,Suffix);
	
	khit = kh_put(hashPos_t, h, Suffix, &hret);
	if (!hret){
	  //Already in the hash
	  // kh_del(hashPos_t, h, k);
	  kh_value(h, khit) = kh_value(h, khit)+1;
	  //verbose(1,"#3b) Already seen %d times\n", kh_value(h, khit));
	  fprintf(stdout,"#3b)%d %x-%x,%s Already seen %d times\n",i,(unsigned int)Prefix,Suffix,buff, kh_value(h, khit));
	}
	else{
	  kh_val(h, khit) = 1;
	}	
      }
      }
    */





