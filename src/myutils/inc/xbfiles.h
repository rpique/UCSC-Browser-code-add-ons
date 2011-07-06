/* xbfiles - Interface with binary files. */
//#include "common.h"
//#include "linefile.h"
//#include "hash.h"
//#include "options.h"

//#include "cutMapperCore.h"
#include <stdio.h>
#include <sys/mman.h>

#define XB_MAGIC 0xCA60B175ul
#define XB_MAGIC_REV 0xB175CA60ul
#define XB_VERSION 0x0001

/* xb File structure.
 *     fixedWidthHeader
 *         magic# 		4 bytes
 *         version              2 bytes
 *	   sizeOf(xbVal_t)	2 bytes
 *         xbStranded           2 bytes (strandedness, 0 no strand, 1 two strands)
 *         xbCount              2 bytes (number of xb elements + 1, We start counting from 1 here!)
 *     xbInfo                   one for each xbCount
 *         nameSize             1byte
 *         xbName               1*nameSize
 *         xbSize               4 bytes (size)
 *     xbRecord                 one for each xbCount (and strand if xbStranded=1)
 *         data                 sizeOf(xbVal_t)*xbSize[xbCount]
 * 
 */

typedef unsigned char xbVal_t;

typedef struct{
  size_t n;
  xbVal_t *a;
}xbVec_t;


/* Maybe it should have its own hash?.. so I can retrieve the chromosome, and etc. */

typedef struct 
/* Holds header and index info from .Xb file. */
{
  //xbList_t *next;
  //boolean isSwapped;	   /* Is byte-swapping needed. I don't want to handle it but just to detect it.*/
  boolean isStranded;      /* There is strand specific information */ 
  bits16 version;	   /* Version of .2bit file */
  bits16 count;	          /* Number of sequences. */
  bits16 typeSize;       /* sizeof(xbVal_t); */
  char **names;        /* Names of the chromosomes */
  unsigned *sizes;     /* Lengths of the crhomosomes */
  xbVec_t *vec;          /* Data containers */ 
  xbVal_t *a;           /* The entire concatenated vector, or NULL if not used */
  void *pmmap;          /* Pointer to the mmaped area, or NULL if not mmaped */
  long long int mmapLength; /* direction to the file */
  int fd;  /* File descriptor index */
}xbList_t;



// Add type size??, so we don't need that 
xbList_t *xbInit(boolean isStranded, bits16 count,char **names, unsigned *sizes){
  xbList_t *xbl=malloc(sizeof(xbList_t));
  int j;
  xbl->version = XB_VERSION;
  xbl->typeSize = sizeof(xbVal_t);
  xbl->count = count;
  xbl->isStranded = isStranded;
  xbl->vec = malloc(sizeof(xbVec_t) * (xbl->count << isStranded));
  xbl->names = malloc(sizeof(char *) * (xbl->count));
  xbl->sizes = malloc(sizeof(unsigned) * (xbl->count));
  verbose(2,"before loop %d\n",count);
  for (j=0;j<xbl->count;j++){
    verbose(2,"%d,%s %d\t",j,names[j],sizes[j]);
    xbl->names[j]=cloneString(names[j]);
    xbl->sizes[j]=sizes[j];
    verbose(2,"%s %d\n",xbl->names[j],xbl->sizes[j]);
  }
  if(isStranded){
    for (j=0;j<xbl->count;j++){ 
      int k=j<<1;  // I think I should change k by j on the sizes.... 
      xbl->vec[k].n=sizes[j]; 
      xbl->vec[k].a=malloc(sizeof(xbVal_t)*sizes[j]);
      xbl->vec[k+1].n=sizes[j]; 
      xbl->vec[k+1].a=malloc(sizeof(xbVal_t)*sizes[j]);
    }
  }
  else{
    for (j=0;j<xbl->count;j++){ 
      xbl->vec[j].n=sizes[j]; 
      xbl->vec[j].a=malloc(sizeof(xbVal_t)*sizes[j]);
    }
  }
  xbl->a=NULL;
  xbl->pmmap=NULL;
  xbl->mmapLength=0;
  return xbl;
}

void xbDestroy(xbList_t *xbl){
  int j;
  
  if(xbl->pmmap==NULL){ //Regular free?
    for (j=0;j<(xbl->count << xbl->isStranded);j++)
      free(xbl->vec[j].a);
  }
  else{// mmaped memory
    munmap(xbl->pmmap,xbl->mmapLength);
  }
  free(xbl->vec);
  free(xbl->sizes);
  free(xbl->names);
  
  free(xbl);
}


// Doesn't need typesize
void xbDump(xbList_t *xbl,char *fName)
/* Write out header portion of twoBit file, including initial
 * index */
{
  FILE *f=fopen(fName,"wb");
  bits32 magic = XB_MAGIC;
  bits16 version = xbl->version;
  bits16 typeSize = xbl->typeSize;
  unsigned long long totalSize = 0;
  int j;
  unsigned char nameSize;
  //long long counter = 0; /* check for 32 bit overflow */
  
  /* Write out fixed parts of header. */
  fwrite(&magic,sizeof(bits32),1,f);
  fwrite(&version,sizeof(bits16),1,f);
  fwrite(&typeSize,sizeof(bits16),1,f);
  fwrite(&(xbl->isStranded),sizeof(boolean),1,f);
  fwrite(&(xbl->count),sizeof(bits16),1,f);

  /* we don't need this I think */
  //  fwrite(&offset,sizeof(bits64),1,f); //PlaceHolder offset beggining of data 
  //  fwrite(&offset,sizeof(bits64),1,f); //PlaceHolder total size of the data.
  
  /* Sequence specific info and calculating offsets...*/
  for (j=0;j<xbl->count;j++){ //chr 0 reserved...
    //    verbose(2,"##Dumping %d,%s %d\n",j,xbl->names[j],xbl->sizes[j]);
    nameSize=strlen(xbl->names[j]);
    fwrite(&nameSize,1,1,f);
    fwrite(xbl->names[j],1,nameSize,f);
    fwrite(&xbl->sizes[j],sizeof(bits32),1,f);
    totalSize += (unsigned long long) xbl->sizes[j];
  }
  verbose(2,"##Dumping data...\n"); //Not prepared for stranded data... !!
  for (j=0;j<xbl->count;j++){
    //verbose(2,"#%d) %d %d\n",j,xbl->typeSize,xbl->sizes[j]);
    verbose(2,"##Dumping %d,%s %d...",j,xbl->names[j],xbl->sizes[j]);
    totalSize=fwrite(xbl->vec[j].a,xbl->typeSize,xbl->sizes[j],f);
    verbose(2,"%lld,ftell%ld\n",totalSize,ftell(f));
    assert(totalSize==xbl->sizes[j]);
    //verbose(2," %lld\n",totalSize);
    //verbose(2,"# ftell=%ld\n",ftell(f));
  }
  verbose(2,"...Completed\n");
  carefulClose(&f);
}


//Shouldn't need type size I think
xbList_t *xbLoadMmap(char *fName)
/* Open file, read in header but not index.  
 * Squawk and die if there is a problem. */
{
  FILE *f=fopen(fName,"r+b");
  xbList_t *xbl=malloc(sizeof(xbList_t));

  int j;
  bits32 sig;
  unsigned long long int totalSize = 0;
  long long int offset = 0;
  unsigned char nameSize;
  int fdi;

  verbose(2,"## Mmaping file %s, PageSize %d\n",fName,getpagesize());
  mustReadOne(f, sig);
  assert(sig==XB_MAGIC);
  /*
  if (sig == twoBitSwapSig)
    isSwapped = tbf->isSwapped = TRUE;
    else if (sig != twoBitSig)
    errAbort("%s doesn't have a valid twoBitSig", fileName);
  */

  fread(&(xbl->version),sizeof(bits16),1,f);
  assert(xbl->version == XB_VERSION);
  fread(&(xbl->typeSize),sizeof(bits16),1,f);
  fread(&(xbl->isStranded),sizeof(boolean),1,f);
  fread(&(xbl->count),sizeof(bits16),1,f);

  verbose(2,"##  There are %d vectors, and %d strands, of typesize %d\n",xbl->count,xbl->isStranded+1,xbl->typeSize);

  xbl->vec = malloc(sizeof(xbVec_t) * (xbl->count << xbl->isStranded));
  xbl->names = malloc(sizeof(char *) * xbl->count);
  xbl->sizes = malloc(sizeof(unsigned) * xbl->count);
  
  /* we don't need this I think */
  //  fwrite(&offset,sizeof(bits64),1,f); //PlaceHolder 
  //  fwrite(&offset,sizeof(bits64),1,f); //PlaceHolder 
  
  /* Sequence specific info and calculating offsets...*/
  for (j=0;j<xbl->count;j++){ //chr 0 reserved...
    fread(&nameSize,1,1,f);
    xbl->names[j]=malloc(nameSize);
    fread(xbl->names[j],1,nameSize,f);
    fread(&(xbl->sizes[j]),sizeof(bits32),1,f);
    totalSize += (bits64) xbl->sizes[j];
    verbose(3,"### %d,%s %d\n",j,xbl->names[j],xbl->sizes[j]);
  }
  offset=ftell(f);
  xbl->mmapLength=offset + (((unsigned long long)totalSize << (xbl->isStranded)) * (unsigned long long)xbl->typeSize);  
  verbose(2,"## At position %lld, data size %lld, total size %lld\n",offset,(long long int)totalSize,xbl->mmapLength);
  // Now we do the mmap!

  verbose(2,"##...mmap-ing data...");
  rewind(f);
  fdi=fileno(f);
 
  //  xbl->pmmap=mmap(0,xbl->mmapLength,PROT_READ,MAP_PRIVATE,fdi,0); PRIVATE
  xbl->pmmap=mmap(0,xbl->mmapLength,PROT_READ|PROT_WRITE,MAP_SHARED,fdi,0);
  xbl->a=(xbl->pmmap+offset);


  if(xbl->isStranded){
    totalSize=0;
    for (j=0;j<xbl->count;j++){ //chr 0 reserved...    
      int k=j<<1;      
      xbl->vec[k].a=&(xbl->a[totalSize]);
      xbl->vec[k].n=xbl->sizes[j];
      totalSize += (unsigned long long) xbl->sizes[j];
      xbl->vec[k+1].a=&(xbl->a[totalSize]);
      xbl->vec[k+1].n=xbl->sizes[j];
      totalSize += (unsigned long long) xbl->sizes[j];
      verbose(3,"###%d x 2)%d,%lld\n",j,xbl->sizes[j],totalSize);
    }
    verbose(2,"...%lld,completed! \n",totalSize);
    verbose(2,"Last Element +:%d ",(int)xbl->vec[(j-1)<<1].a[xbl->sizes[j-1]-1]);
    verbose(2,"Last Element -:%d\n",(int)xbl->vec[((j-1)<<1) + 1].a[xbl->sizes[j-1]-1]);
  }
  else{
    totalSize=0;
    for (j=0;j<xbl->count;j++){ //chr 0 reserved...    
      xbl->vec[j].a=&(xbl->a[totalSize]);
      xbl->vec[j].n=xbl->sizes[j];
      totalSize += (unsigned long long) xbl->sizes[j];
      verbose(3,"###%d x 2)%d,%lld\n",j,xbl->sizes[j],totalSize);
    }
    verbose(2,"...%d,completed! ",j);
    verbose(2,"Last Element:%d\n",(int)xbl->vec[j-1].a[xbl->sizes[j-1]-1]);
  }



  //fclose(f);
  
  return xbl;
}
  

// Add type size??, so we don't need that 
xbList_t *xbInitMmap(boolean isStranded, bits16 count,char **names, unsigned *sizes, char *fName){
  //  xbList_t *xbl;//=malloc(sizeof(xbList_t));

  xbVal_t aux=0;
  
  FILE *f=fopen(fName,"wb");
  bits32 magic = XB_MAGIC;
  bits16 version = XB_VERSION;
  bits16 typeSize = sizeof(xbVal_t);
  unsigned long long totalSize = 0;
  int j;
  unsigned char nameSize;
  //long long counter = 0; /* check for 32 bit overflow */
  
  verbose(2,"Initializing header of %s...\n",fName);
  /* Write out fixed parts of header. */
  fwrite(&magic,sizeof(bits32),1,f);
  fwrite(&version,sizeof(bits16),1,f);
  fwrite(&typeSize,sizeof(bits16),1,f);
  fwrite(&isStranded,sizeof(boolean),1,f);
  fwrite(&count,sizeof(bits16),1,f);

  /* we don't need this I think */
  //  fwrite(&offset,sizeof(bits64),1,f); //PlaceHolder offset beggining of data 
  //  fwrite(&offset,sizeof(bits64),1,f); //PlaceHolder total size of the data.
  
  /* Sequence specific info and calculating offsets...*/
  for (j=0;j<count;j++){ //chr 0 reserved...
    //    verbose(2,"##Dumping %d,%s %d\n",j,xbl->names[j],xbl->sizes[j]);
    nameSize=strlen(names[j]);
    fwrite(&nameSize,1,1,f);
    fwrite(names[j],1,nameSize,f);
    fwrite(&sizes[j],sizeof(bits32),1,f);
    totalSize += (unsigned long long) sizes[j];
  }
  verbose(2,"..completed\n");

  verbose(2,"Initializing data space for %lld bytes...\n",(totalSize<<isStranded)*typeSize);
  fseek(f,((totalSize<<isStranded) -1)*typeSize, SEEK_CUR);
  fwrite(&aux,typeSize,1,f);
  fclose(f);
  verbose(2,"..completed\n");
  
  return xbLoadMmap(fName);
}

  

/* struct twoBitFile *xbOpen(char *fileName) */
/* /\* Open file, read in header and index.   */
/*  * Squawk and die if there is a problem. *\/ */
/* { */
/* struct twoBitFile *tbf = twoBitOpenReadHeader(fileName); */
/* struct twoBitIndex *index; */
/* boolean isSwapped = tbf->isSwapped; */
/* int i; */
/* struct hash *hash; */
/* FILE *f = tbf->f; */

/* /\* Read in index. *\/ */
/* hash = tbf->hash = hashNew(digitsBaseTwo(tbf->seqCount)); */
/* for (i=0; i<tbf->seqCount; ++i) */
/*     { */
/*     char name[256]; */
/*     if (!fastReadString(f, name)) */
/*         errAbort("%s is truncated", fileName); */
/*     lmAllocVar(hash->lm, index); */
/*     index->offset = readBits32(f, isSwapped); */
/*     hashAddSaveName(hash, name, index, &index->name); */
/*     slAddHead(&tbf->indexList, index); */
/*     } */
/* slReverse(&tbf->indexList); */
/* return tbf; */
/* } */

