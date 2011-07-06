#ifndef CUT_MAPPER_CORE_H
#define CUT_MAPPER_CORE_H

// I could avoid this to libraries including hear if I use void pointers..
#include "common.h"
#include "linefile.h"


#pragma pack(push)  /* push current alignment to stack */
#pragma pack(1)     /* set alignment to 1 byte boundary */
typedef struct{
  unsigned char chr; // index to the sequence list, 1-254, 0 (reserved), 255 (Multiple ocurring times)
  int pos;  //Position + for the plus strand, negative numbers for the negative strand. 
            //if chr=255, pos is the index to the array holding all occurrences...
}chrpos_t;

typedef struct{
  int pos;
  char ref;
  char alt;
}snp_t;

typedef struct{
  int pos;          // Position on the genome of the firt base of the indel. 0-based...
  short int reflen; // Length of the reference sequence
  short int altlen; // Length of the alternate sequence
  char *ref;        // Reference sequence
  char *alt;        // Alternate sequence to replace.  // We will use the convention that the first base of the Reference should stay.
}indel_t;

#pragma pack(pop)   /* restore original alignment from stack */


#include "khash.h"

//KHASH_MAP_INIT_INT(32, int)

KHASH_MAP_INIT_INT(hCount_t, int)            // Hash number of times read seen. 
KHASH_MAP_INIT_INT(hashPos_t, chrpos_t)      // 
KHASH_MAP_INIT_STR(hashChr_t, unsigned char)

#include "kvec.h"

typedef kvec_t(chrpos_t)repvec_t;  //Repeats of a kmer
typedef kvec_t(repvec_t)replist_t; //List of all kmer repeat vectors
typedef kvec_t(snp_t)snpvec_t;     //Vector of snps in a chromosome
typedef kvec_t(snpvec_t)snplist_t;    //List of snps vectors.  	
typedef kvec_t(indel_t)indelvec_t;      //Vector of indels in a chromosome
typedef kvec_t(indelvec_t)indellist_t;    //List of indel vectors.  	


		
#define kmerSize 20 
#define prefixSize 4
#define suffixSize 16
#define KmerMask  0x000000FFFFFFFFFF
#define getPrefix(Kmer) ((Kmer&0x000000FFFFFFFFFF)>>32)
#define getSuffix(Kmer) ((Kmer)&0xFFFFFFFF)


/* Function headers shared by *.c */

int encodeACGT(char b);

unsigned long int packKmer(char *s,int k);

void unPackKmer(unsigned long int kmer,int k,char *s);

khash_t(hashPos_t) *hashFileOpen(char *fileName);

replist_t *repeatFileOpen(char *fileName, replist_t *replistp);

boolean checkQual(char *file, int line, char *s);

boolean checkSeqName(char *file, int line, char *s, char firstChar, char *name);

boolean wantNewLine(struct lineFile *lf, char *file, int line, char **row, char *msg);
		   
void hashFastqFile(char *inFastq, khash_t(hCount_t) **h);

khash_t(hashChr_t) *loadChromSizes(char *chromFile, char ***chromNames, unsigned **chromSizes);

int loadSnpFile(char *snpFile, khash_t(hashChr_t) *hChr,  snplist_t *snplistp, indellist_t *indellistp);


//int repeatFileOpen(char *fileName, replist_t *replistp);



#endif /* CUT_MAPPER_CORE_H */

