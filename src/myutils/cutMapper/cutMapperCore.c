/* Putting this in a core library... cutMapperCore.c */

#include <sys/mman.h>

#include "cutMapperCore.h"

#include "common.h"
#include "linefile.h"


int encodeACGT(char c){
  c=c&0xDF; // 11011111 sets letter to uppercase
  switch(c){
  case 'A':
    return 0;
  case 'C':
    return 1;
  case 'G':
    return 2;
  case 'T':
    return 3;
  default:
    return -1;
  }
}

/*
int encodeACGT(char b){
  if(b=='A' || b=='a') {
    return 0;
  } else if(b=='C' || b=='c' ) {
    return 1;
  } else if(b=='G' || b=='g' ) {
    return 2;
  } else if(b=='T' || b=='t' ) {
    return 3;
  }
  return -1;
}
*/

unsigned long int packKmer(char *s,int k){
  int Base,j;
  unsigned long int kmer=0;
  //  assert((k<<1)<sizeof(unigned long int));
  for(j=0;j<k;j++){
    Base=encodeACGT(s[j]);
    if(Base==-1)
      return 0xFFFFFFFFFFFFFFFF; //
    kmer=(kmer<<2)+Base;
  }
  return kmer;
}

void unPackKmer(unsigned long int kmer,int k,char *s){
  int j;
  char ACGT[]={"ACGT"};
  for(j=k-1;j>=0;j--){
    s[j]=ACGT[(kmer & 0x03)];
    kmer=(kmer>>2);
  }
  s[k]='\0';
}

/*
// If I was to do this function as a define then it would be general purpose, as in khash.h
khash_t(hashPos_t) *hashFileOpen(char *fileName){
  chrpos_t auxPos;  
  uint32_t auxKey;
  khint_t size;
  khiter_t khit;
  int hret;
  FILE *f;
  int i=0;
  
  khash_t(hashPos_t) *h = kh_init(hashPos_t);  
  
  f=mustOpen(fileName,"rb");
  // We save the size of the hash at the beginning of the hash file.
  fread(&size,sizeof(khint_t),1,f);
  // We set hash size.. so it does not need to reincrease...
  kh_resize(hashPos_t, h, size);
  verbose(1," Reading hash table from %s with %d entries\n",fileName,size);

  //while(!feof(f)){
  while((size>0) & !feof(f)){
    i++;
    fread(&auxKey,sizeof(uint32_t),1,f);
    fread(&auxPos,sizeof(chrpos_t),1,f);
    khit = kh_put(hashPos_t, h , auxKey, &hret);
    if(hret!=1)errAbort("ERROR: %d\t%d\t%x\t%d\t%d\t%d\t%d\n",i,size,auxKey,(int) auxPos.chr,(int) auxPos.pos,hret,feof(f));
    //assert(hret==1); //It should be unique in the dump!
    kh_value(h , khit)=auxPos;
    size--;
  }
  
  assert(size==0);
  carefulClose(&f);
  verbose(1," Hash read completed succesfully.\n");

  return h;
}

*/
// If I was to do this function as a define then it would be general purpose, as in khash.h
// Read in one gulp, faster... ?
// Not much faster, only a little bit, maybe I can try mmap
khash_t(hashPos_t) *hashFileOpen(char *fileName){
  chrpos_t auxPos;  
  uint32_t auxKey;
  char *buff;
  //char *pmap;
  khint_t size;
  khiter_t khit;
  int hret;
  int tsize=(sizeof(uint32_t)+sizeof(chrpos_t));
  FILE *f=mustOpen(fileName,"r+b");
  //int fdi=fileno(f);
  int i=0;
  
  khash_t(hashPos_t) *h = kh_init(hashPos_t);  
  
  // We save the size of the hash at the beginning of the hash file.
  fread(&size,sizeof(khint_t),1,f);
  // We set hash size.. so it does not need to reincrease...
  kh_resize(hashPos_t, h, size);
  

  //What if I do it with an mmap ???
  
  buff = (char *)calloc(tsize, size);
  verbose(2,"# Reading hash table from %s with %d entries...\n", fileName, size);
  //Read in one big gulp... 
  assert(fread(buff,tsize,size,f)==size);
  carefulClose(&f);
  
  /*  
  rewind(f);
  verbose(2,"# Reading hash table from %s with %d entries using mmap...", fileName, size);
  pmap = (char *)mmap(0,((long long int)tsize*size+sizeof(khint_t)),PROT_READ,MAP_PRIVATE,fdi,0);
  buff=pmap+sizeof(khint_t);
  verbose(2,"... mmap OK\n");
  */
    
  //while(!feof(f))
  while(size>0){
    auxKey=*((uint32_t *)(buff+i*tsize));
    auxPos=*((chrpos_t *)(buff+i*tsize+sizeof(uint32_t)));
    khit = kh_put(hashPos_t, h , auxKey, &hret);
    //    verbose(2,"Key %u, Pos %d,%d, khit%d, hret%d\n",auxKey,(int)auxPos.chr
    if(hret!=1)errAbort("ERROR: %d\t%d\t%x\t%d\t%d\t%d\t%d\n",i,size,auxKey,(int) auxPos.chr,(int) auxPos.pos,hret,khit);
    //assert(hret==1); //It should be unique in the dump!
    kh_value(h , khit)=auxPos;
    size--;
    i++;
  }
  free(buff);
  /*
  //unmap
  munmap(pmap,tsize*size+sizeof(khint_t));
  carefulClose(&f);
  */


  verbose(3,"# ... completed succesfully.\n");
  verbose(2,"# khash: n_buckets=%d, size=%d, n_occupied=%d, upper_bound=%d \n",
	  kh_n_buckets(h),kh_size(h),((h)->n_occupied),((h)->upper_bound));

  return h;
}


replist_t *repeatFileOpen(char *fileName, replist_t *replistp)
{
  FILE *f;
  size_t lsize;
  size_t vsize;
  size_t cumsize=0;
  int k=0;
  repvec_t *repvecp;
  
  f=mustOpen(fileName,"rb");
  //  fwrite(&kv_size(replist),sizeof(size_t),1,f);
  fread(&lsize,sizeof(size_t),1,f);
  kv_init_size(repvec_t, *replistp, lsize);
  
  verbose(3,"# Reading repeats from %s with %d entries...\n",fileName, (unsigned)lsize);
  for (k = 0; k < lsize; ++k){
    fread(&vsize,sizeof(size_t),1,f);
    //kv_destroy(*repvecp); init new unit
    //// kv_push_empty(repvec_t,*replistp); // No because I alread created the unit!
    repvecp = &kv_A(*replistp,k);
    kv_init_size(chrpos_t, *repvecp, vsize);
    kv_resize(chrpos_t, *repvecp, vsize);
    assert(fread(repvecp->a,sizeof(chrpos_t),vsize,f)==vsize);
    cumsize+=vsize;
  }
  carefulClose(&f);
  verbose(3,"# ...completed succesfully. (%d repeats)\n", (unsigned)cumsize);
  
  return replistp;
}



/*
boolean checkSeq(char *file, int line, char *row, char *s, char *name)
// Return TRUE if string has non-zero length and contains only chars [ACGTNacgtn0-3]
// Othewise print warning that name column is empty and return FALSE
{
  verbose(3,"[%s %3d] inputLine=%d %s seq(%s) [%s]\n", __func__, __LINE__, line, name, s, row);
  int i;
  for ( i = 0; s[i] ; ++i){
    if (!dnaChars[(int)s[i]]){
      if (s==row)
	warn("Error [file=%s, line=%d]: invalid DNA chars in %s(%s)", file, line, name, s);
      else
	warn("Error [file=%s, line=%d]: invalid DNA chars in %s(%s) [%s]", file, line, name, s, row);
      return FALSE;
    }
  }
if (i == 0)
    {
    if(privateData)  // PrivateData means sequence should be empty
        return TRUE;
    if (s==row)
	warn("Error [file=%s, line=%d]: %s empty", file, line, name);
    else
	warn("Error [file=%s, line=%d]: %s empty in line [%s]", file, line, name, row);
    return FALSE;
    }
else if(privateData) { // PrivateData means sequence should be empty
    if (s==row)
        warn("Error [file=%s, line=%d]: %s is not empty but this should be private data", file, line, name);
    else
        warn("Error [file=%s, line=%d]: %s  is not empty but this should be private data in line [%s]", file, line, name, row);
    return FALSE;
    }
return TRUE;
}
*/

boolean wantNewLine(struct lineFile *lf, char *file, int line, char **row, char *msg)
{
  boolean res = lineFileNext(lf, row, NULL);
  if (!res)
    warn("Error [file=%s, line=%d]: %s not found", file, line, msg);
  return res;
}

boolean checkSeqName(char *file, int line, char *s, char firstChar, char *name)
// Return TRUE if string has non-zero length and contains only seqName[] chars
// Othewise print warning that seqName is empty and return FALSE
{
  //int i;
  if (s[0] == 0){
    warn("Error [file=%s, line=%d]: %s empty [%s]", file, line, name, s);
    return FALSE;
  }
  else if (s[0] != firstChar){
    warn("Error [file=%s, line=%d]: %s first char invalid (got '%c', wanted '%c') [%s]",
	 file, line, name, s[0], firstChar, s);
    return FALSE;
  }
  /*
  for ( i = 1; s[i] ; ++i){
    if (s[i] == ' ')
      break;
    if (!seqName[(int)s[i]]){
      warn("Error [file=%s, line=%d]: invalid %s chars in [%s]", file, line, name, s);
      return FALSE;
    }
    }*/
  return TRUE;
}

boolean checkQual(char *file, int line, char *s)
// Return TRUE if string has non-zero length and contains only qualChars[] chars
// Othewise print warning that quality is empty and return FALSE
{
  /*int i;
    int qual;
  for ( i = 0; s[i] ; ++i){
    s[i]=s[i]-64; //converting to quality sore [0..62]
    if ((s[i] < 0) || (s[i] > 62)){
      warn("Error [file=%s, line=%d]: invalid quality chars", file, line);
      return FALSE;
    }
  }
  if (i == 0){
    warn("Error [file=%s, line=%d]: quality empty [%s]", file, line, s);
    return FALSE;
    }*/
  return TRUE;
}


// fastq:
// @NINA_1_FC30G3VAAXX:5:1:110:908
// ATCGTCAGGTGGGATAATCCTTACCTTTTCCTCCTC
// +NINA_1_FC30G3VAAXX:5:1:110:908
// aa`]`a`XQ^VQQ^`aaaaaaa^[[ZG[aXUX[[[X

void hashFastqFile(char *inFastq, khash_t(hCount_t) **h)
/* hashFastqFile - Hash all the reads of the inptut file in fastq format int the hash tables.. */
{
  struct lineFile *lf = lineFileOpen(inFastq, TRUE);
  char *seqName, *seq, *qName, *qual;
  int line = 0;
  int numKmers = 0;
  int numRepKmers =0;
  int numSeq =0;
  int numWithN=0;
  int errs = 0;
  boolean startOfFile = TRUE;

  unsigned long int kmer;
  unsigned char Prefix=0;
  unsigned int Suffix;
  
  khiter_t khit;
  int hret;
  int j;
  
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
    if ( !(checkSeqName(inFastq, line, seqName, '@', "sequence name")
	   && (wantNewLine(lf, inFastq, ++line, &seq, "fastq sequence line"))
	   && (seq=cloneString(seq))
	   && (strlen(seq)>0) //checkSeq(inFastq, line, seq, seq, "sequence")
	   && (wantNewLine(lf, inFastq, ++line, &qName, "fastq sequence name (quality line)"))
	   && checkSeqName(inFastq, line, qName, '+', "quality name")
	   && (wantNewLine(lf, inFastq, ++line, &qual, "quality line"))
	   && checkQual(inFastq, line, qual)
	   )){
      //      if (printFailLines)
      printf("%s\n%s\n%s\n%s\n", seqName, seq, qName, qual);
      if (++errs >= 1)//maxErrors)
	errAbort("Aborting .. found %d errors\n", errs);
    }
    else{
      // OK now we are fine....
      // Maybe I should use quality to flag out some bases...
      kmer=packKmer(seq,kmerSize);
      if(kmer != -1){ // otherWise check how many N, and hash them in a different way... ? 4 times???
	Prefix=getPrefix(kmer); // ((ForKmer&kmermask)>>32);
	Suffix=getSuffix(kmer); // &0xFFFFFFFF;
	khit = kh_put(hCount_t, h[Prefix] , Suffix, &hret);
	if (!hret){
	  //Already in the hash
	  kh_value(h[Prefix], khit)++;
	  numRepKmers++;
	}
	else{
	  kh_value(h[Prefix], khit)=1;
	  numKmers++;
	}
      }
      else{  
	numWithN++;
      }
      ++numSeq;
      freeMem(seq);
    }
  }
  
  verbose(1, "#\tprocessed %d lines  %d reads\n",line,numSeq);
  verbose(1, "#\tTotal number of reads with >0 N: %d\n",numWithN);
  
  line=0;
  for(j=0;j<256;j++) 
    line+=kh_size(h[j]);  
  
  verbose(1, "#\tTotal number of k-mers %d, %d \n",line,numKmers); //WHAT THE HELL IS GOING ON!!!
  verbose(1, "#\tTotal number of rep k-mers %d \n",numRepKmers); 
  
  line=0; //Traverse the hashes
  for(j=0;j<256;j++) 
    for (khit = kh_begin(h[j]); khit < kh_end(h[j]); ++khit)
      if (kh_exist(h[j], khit))
	line+=kh_val(h[j], khit); 
  verbose(1, "#\tTotal number of reads potentially mappable %d \n",line); 


  // I WOULD LIKE TO PUT IT IN A DIFFERENT FUNCTION
  // I could check if the last two letters  contain the beggining of the adapter... 
  // 5.-P-TCGTATGCCGTCTTCTGCTTG
  // 3.-NNAGCATACGGCAGAAGACGAAC
  // TCGTATGCCGTCTTCTGCTTGA, CGTATGCCGTCTTCTGCTTGAA, *TCGTATGCCGTCTTCTGCTTG
  // Removing adapter sequence
  seq="TCGTATGCCGTCTTCTGCTT";
  kmer=packKmer(seq,kmerSize);
  Prefix=getPrefix(kmer); // ((ForKmer&kmermask)>>32);
  Suffix=getSuffix(kmer); // &0xFFFFFFFF;
  khit = kh_get(hCount_t, h[Prefix] , Suffix);
  if(kh_exist(h[Prefix], khit)){
    verbose(1, "#\tAdapter kmer %s seen %d times\n",seq,kh_val(h[Prefix],khit));
    kh_del(hCount_t, h[Prefix], khit);
  }
  seq="CGTATGCCGTCTTCTGCTTG";
  kmer=packKmer(seq,kmerSize);
  Prefix=getPrefix(kmer); // ((ForKmer&kmermask)>>32);
  Suffix=getSuffix(kmer); // &0xFFFFFFFF;
  khit = kh_get(hCount_t, h[Prefix] , Suffix);  
  if(kh_exist(h[Prefix], khit)){
    verbose(1, "#\tAdapter kmer %s seen %d times\n",seq,kh_val(h[Prefix],khit));
    kh_del(hCount_t, h[Prefix], khit);
  }
  seq="GTATGCCGTCTTCTGCTTGA";
  kmer=packKmer(seq,kmerSize);
  Prefix=getPrefix(kmer); // ((ForKmer&kmermask)>>32);
  Suffix=getSuffix(kmer); // &0xFFFFFFFF;
  khit = kh_get(hCount_t, h[Prefix] , Suffix);
  if(kh_exist(h[Prefix], khit)){
    verbose(1, "#\tAdapter kmer %s seen %d times\n",seq,kh_val(h[Prefix],khit));
    kh_del(hCount_t, h[Prefix], khit);
  }
  seq=strdup("ATCGTATGCCGTCTTCTGCT");
  kmer=packKmer(seq,kmerSize);
  Prefix=getPrefix(kmer); // ((ForKmer&kmermask)>>32);
  Suffix=getSuffix(kmer); // &0xFFFFFFFF;
  for(j=0;j<4;j++){
    Prefix=getPrefix(kmer); // ((ForKmer&kmermask)>>32);
    Prefix=Prefix+(unsigned char)(j<<6);
    //fprintf(stdout,"%d %d\n",(int)Prefix,(int)j);
    unPackKmer((((unsigned long int)Prefix)<<32) + (unsigned long int)Suffix, 20, seq);
    khit = kh_get(hCount_t, h[Prefix] , Suffix);
    if(kh_exist(h[Prefix], khit)){
      verbose(1, "#\tAdapter kmer %s seen %d times\n",seq,kh_val(h[Prefix],khit));
      kh_del(hCount_t, h[Prefix], khit);
    }
  }
  free(seq);
  seq=strdup("AATCGTATGCCGTCTTCTGC");
  kmer=packKmer(seq,kmerSize);
  Prefix=getPrefix(kmer); // ((ForKmer&kmermask)>>32);
  Suffix=getSuffix(kmer); // &0xFFFFFFFF;
  for(j=0;j<16;j++){
    Prefix=getPrefix(kmer); // ((ForKmer&kmermask)>>32);
    Prefix=Prefix+(unsigned char)(j<<4);
    unPackKmer((((unsigned long int)Prefix)<<32) + (unsigned long int)Suffix, 20, seq);
    khit = kh_get(hCount_t, h[Prefix] , Suffix);
    if(kh_exist(h[Prefix], khit)){
      verbose(1, "#\tAdapter kmer %s seen %d times\n",seq,kh_val(h[Prefix],khit));
      kh_del(hCount_t, h[Prefix], khit);
    }
  }
  free(seq);
  //Adapter K-mer =  TCGTATGCCGTCTTCTGCTT

  line=0; //Traverse the hashes
  for(j=0;j<256;j++) 
    for (khit = kh_begin(h[j]); khit < kh_end(h[j]); ++khit)
      if (kh_exist(h[j], khit))
	line+=kh_val(h[j], khit); 
  verbose(1, "#\tTotal number of reads potentially mappable %d \n",line); 

  //close lf file...
  lineFileClose(&lf);
}

khash_t(hashChr_t) *loadChromSizes(char *chromFile, char ***chromNames, unsigned **chromSizes)
{
  FILE *f=mustOpen(chromFile,"r");
  int chrnum=0;
  char buff[256];
  char empty[]="---";
  //char *buff2;
  unsigned size;
  khash_t(hashChr_t) *hChr= kh_init(hashChr_t); 
  
  khiter_t k;
  int hret;

  kvec_t(unsigned) kvChromSizes;
  kvec_t(char *)   kvChromNames;
  
  kv_init(kvChromSizes);
  kv_init(kvChromNames);

  //Pushing the empty crhomosome.
  //buff=cloneString(empty);
  k = kh_put(hashChr_t, hChr , cloneString(empty), &hret);
  kv_push(char *,kvChromNames, (char *)kh_key(hChr, k));
  kv_push(unsigned,kvChromSizes,0);
    
  while(!feof(f)){
    chrnum++;
    assert(chrnum<256);
    fscanf(f,"%s\t%u\t", buff, &size);
    verbose(2,"#%d) %s\t%u\n",chrnum, buff, size);
    //buff2=cloneString(buff);
    //kv_push(char *,kvChromNames, buff2 );
    k = kh_put(hashChr_t, hChr , cloneString(buff), &hret);
    assert(hret==1);
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


