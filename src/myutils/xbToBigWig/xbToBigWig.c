/* xbToBigWig - Converts a xb file to a bw file. */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"

static char const rcsid[] = "$Id: newProg.c,v 1.30 2010/03/24 21:18:33 hiram Exp $";

void usage()
/* Explain usage and exit. */
{
  errAbort(
  "xbToBigWig - Converts a xb file to a bw file\n"
  "usage:\n"
  "   xbToBigWig in.xb out.F.bw out.R.bw\n"
  "in.xb - is the xb binary file \n"
  "out.F.bw - bw file with the forward strand of the xb file \n"
  "out.R.bw - bw file with the reverse strand of the xb file \n"
  "options:\n"
  //  "   -smoothing = 0 (default) no smoothing, 1 - (0.25 0.5 0.25)\n" //, 2 - (1/9,2/9,3/9,2/9,1/9)"
  "   -verbose = verbosity lebel (1 default)\n"
  );
}

static struct optionSpec options[] = {
   {NULL, 0},
};



struct hash *hashChromSizesFromXbList(xbList_t *xbl)
/* Read xbl chrom struct into hash keyed by chrom. */
{
  struct hash *hash = hashNew(0);
  //struct lineFile *lf = lineFileOpen(fileName, TRUE);
  //char *row[2];
  while (lineFileRow(lf, row))
    hashAddInt(hash, row[0], sqlUnsigned(row[1]));

  //lineFileClose(&lf);
  return hash;
}


static void parseBedGraphSection(struct lineFile *lf, boolean clipDontDie, 
				 struct hash *chromSizeHash, struct lm *lm, 
				 int itemsPerSlot, struct bwgSection **pSectionList)
/* Parse out bedGraph section until we get to something that is not in bedGraph format. */
{
  /* Set up hash and list to store chromosomes. */
  struct hash *chromHash = hashNew(0);
  struct bedGraphChrom *chrom, *chromList = NULL;

  /* Collect lines in items on appropriate chromosomes. */
  struct bwgBedGraphItem *item;
  char *line;
  while (lineFileNextReal(lf, &line))
    {
      /* Check for end of section. */
      if (stepTypeLine(line))
        {
	  lineFileReuse(lf);
	  break;
	}

      /* Parse out our line and make sure it has exactly 4 columns. */
      char *words[5];
      int wordCount = chopLine(line, words);
      lineFileExpectWords(lf, 4, wordCount);

      /* Get chromosome. */
      char *chromName = words[0];
      chrom = hashFindVal(chromHash, chromName);
      if (chrom == NULL)
        {
	  lmAllocVar(chromHash->lm, chrom);
	  hashAddSaveName(chromHash, chromName, chrom, &chrom->name);
	  chrom->size = hashIntVal(chromSizeHash, chromName);
	  slAddHead(&chromList, chrom);
	}

      /* Convert to item and add to chromosome list. */
      lmAllocVar(lm, item);
      item->start = lineFileNeedNum(lf, words, 1);
      item->end = lineFileNeedNum(lf, words, 2);
      item->val = lineFileNeedDouble(lf, words, 3);

      /* Do sanity checking on coordinates. */
      if (item->start > item->end)
        errAbort("bedGraph error: start (%u) after end line (%u) %d of %s.", 
		 item->start, item->end, lf->lineIx, lf->fileName);
      if (item->end > chrom->size)
	{
	  warn("bedGraph error line %d of %s: chromosome %s has size %u but item ends at %u",
	       lf->lineIx, lf->fileName, chrom->name, chrom->size, item->end);
	  if (!clipDontDie)
	    noWarnAbort();
	}
      else
	{
	  slAddHead(&chrom->itemList, item);
	}
    }
  slSort(&chromList, bedGraphChromCmpName);
  
  for (chrom = chromList; chrom != NULL; chrom = chrom->next)
    {
      slSort(&chrom->itemList, bwgBedGraphItemCmp);

      /* Break up into sections of no more than items-per-slot size. */
      struct bwgBedGraphItem *startItem, *endItem, *nextStartItem = chrom->itemList;
    for (startItem = chrom->itemList; startItem != NULL; startItem = nextStartItem)
	{
	/* Find end item of this section, and start item for next section.
	 * Terminate list at end item. */
	int sectionSize = 0;
	int i;
	endItem = startItem;
	for (i=0; i<itemsPerSlot; ++i)
	    {
	    if (nextStartItem == NULL)
		break;
	    endItem = nextStartItem;
	    nextStartItem = nextStartItem->next;
	    ++sectionSize;
	    }
	endItem->next = NULL;

	/* Fill in section and add it to section list. */
	struct bwgSection *section;
	lmAllocVar(lm, section);
	section->chrom = cloneString(chrom->name);
	section->start = startItem->start;
	section->end = endItem->end;
	section->type = bwgTypeBedGraph;
	section->items.bedGraphList = startItem;
	section->itemCount = sectionSize;
	slAddHead(pSectionList, section);
	}
    }

  /* Free up hash, no longer needed. Free's chromList as a side effect since chromList is in 
   * hash's memory. */
  hashFree(&chromHash);
  chromList = NULL;
}



struct bwgSection *bwgParseWig(char *fileName, boolean clipDontDie, struct hash *chromSizeHash,
			       int maxSectionSize, struct lm *lm)
/* Parse out ascii wig file - allocating memory in lm. */
{
  struct lineFile *lf = lineFileOpen(fileName, TRUE);
  char *line;
  struct bwgSection *sectionList = NULL;

  //This is prepared for both bedgraph and wig, I just need one...
  while (lineFileNextReal(lf, &line)){
    verbose(2, "processing %s\n", line);
    if (stringIn("chrom=", line))
      parseSteppedSection(lf, clipDontDie, chromSizeHash, line, lm, maxSectionSize, &sectionList);
    else{
      /* Check for bed... */
      char *dupe = cloneString(line);
      char *words[5];
      int wordCount = chopLine(dupe, words);
      if (wordCount != 4)
	errAbort("Unrecognized line %d of %s:\n%s\n", lf->lineIx, lf->fileName, line);
      
      /* Parse out a bed graph line just to check numerical format. */
      char *chrom = words[0];
      int start = lineFileNeedNum(lf, words, 1);
      int end = lineFileNeedNum(lf, words, 2);
      double val = lineFileNeedDouble(lf, words, 3);
      verbose(2, "bedGraph %s:%d-%d@%g\n", chrom, start, end, val);

      /* Push back line and call bed parser. */
      lineFileReuse(lf);
      parseBedGraphSection(lf, clipDontDie, chromSizeHash, lm, maxSectionSize, &sectionList);
    }
  }
  slSort(&sectionList, bwgSectionCmp);
  
  /* Check for overlap. */
  struct bwgSection *section, *nextSection;
  for (section = sectionList; section != NULL; section = nextSection){
    nextSection = section->next;
    if (nextSection != NULL){
      if (sameString(section->chrom, nextSection->chrom)){
	if (section->end > nextSection->start){ 
	  //RPR horrible hack
	  //errAbort("There's more than one value for %s base %d (in coordinates that start with 1).\n",
	  //  section->chrom, nextSection->start+1);
	}
      }
    }
  }
return sectionList;
}



void xbToBigWig(char *xbFile, char *outFrwFile, char *outRevFile)
/* xbToBigWig - Converts a xb file to a bw file. */
{
  xbList_t *xbl=xbLoadMmap(xbFile);
  struct hash *chromSizeHash = hashChromSizesFromXbList(xbl);
  struct lm *lm = lmInit(0);

  struct bwgSection *sectionList = bwgParseWig(inName, clipDontDie, chromSizeHash, itemsPerSlot, lm);
  if (sectionList == NULL)
    errAbort("%s is empty of data", xbFile);

  bwgCreate(sectionList, chromSizeHash, blockSize, itemsPerSlot, compress, outName);

  lmCleanup(&lm);
}

int main(int argc, char *argv[])
/* Process command line. */
{
  optionInit(&argc, argv, options);
  if (argc != 3)
    usage();
  xbToBigWig(argv[1],argv[2]);
  return 0;
}
