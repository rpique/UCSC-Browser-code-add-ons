/* trackDbPatch - Patch files in trackDb with a specification .ra file that has db, track, and file fields that say where to apply the patch, and other fields that are patched in.. */
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "ra.h"
#include "portable.h"

static char const rcsid[] = "$Id: trackDbPatch.c,v 1.2 2009/11/23 07:39:46 kent Exp $";

char *clPatchDir = NULL;
char *clKey = "track";
boolean clFirstFile = FALSE;

void usage()
/* Explain usage and exit. */
{
errAbort(
  "trackDbPatch - Patch files in trackDb with a specification .ra file that has db, track, and file fields that say where to apply the patch, and other fields that are patched in.\n"
  "usage:\n"
  "   trackDbPatch patches.ra backupDir\n"
  "options:\n"
  "   -test=patchDir - rather than doing patches in place, write patched output to this dir\n"
  "   -key=tagName - use tagName as key.  Default '%s'\n"
  "   -firstFile - when a patch can go to multiple files apply it to first rather than last file\n"
  , clKey
  );
}

static struct optionSpec options[] = {
   {"test", OPTION_STRING},
   {"key", OPTION_STRING},
   {"firstFile", OPTION_BOOLEAN},
   {NULL, 0},
};


/* Program first reads in patchs to list of raPatches.  Then it creates a list of filePatches so that it can do all patches
 * in one pile at once. */

struct raPatch
/* A patch record. */
    {
    struct raPatch *next;
    char *db;		/* Database. */
    char *track;	/* Track. */
    struct slName *fileList;	/* List of files. */
    struct slPair *tagList;	/* List of tags to merge in. */
    };

struct fileToPatch 
/* A file needing patching */
    {
    struct fileToPatch *next;
    char *fileName;		/* Name of file. */
    struct slRef *patchList;	/* References to raPatches associated with file.  */
    };

struct fileToPatch *groupPatchesByFiles(struct raPatch *patchList, boolean firstFile)
/* Make fileToPatch list that covers all files in patchList. If lastFile is set will apply patch to first (as 
 * opposed to the usual last) file in list of files  associated with a patch. */
{
struct fileToPatch *fileList = NULL, *file;
struct hash *fileHash = hashNew(0);
struct raPatch *patch;

for (patch = patchList; patch != NULL; patch = patch->next)
    {
    struct slName *fileName = patch->fileList;
    if (!firstFile)
        fileName = slLastEl(patch->fileList);
    assert(fileName);
    file = hashFindVal(fileHash, fileName->name);
    if (file == NULL)
        {
	AllocVar(file);
	file->fileName = cloneString(fileName->name);
	slAddHead(&fileList, file);
	hashAdd(fileHash, file->fileName, file);
	}
    refAdd(&file->patchList, patch);
    }

/* Straighten out all lists. */
slReverse(&fileList);
for (file = fileList; file != NULL; file = file->next)
    slReverse(&file->patchList);
hashFree(&fileHash);
return fileList;
}


static struct slName *makeFileList(char *filesAndPos)
/* Convert something that looks like "file # file #" to a list of files.  This
 * writes zeroes into the filesAndPos input. */
{
struct slName *list = NULL;
char *word;
while ((word = nextWord(&filesAndPos)) != NULL)
    {
    slNameAddTail(&list, word);
    word = nextWord(&filesAndPos);
    if (word == NULL && !isdigit(word[0]))
        errAbort("Expecting number in makeFileList, got %s", word);
    }
return list;
}

struct raPatch *raPatchReadOne(struct lineFile *lf)
/* Read next patch from lineFile. */
{
struct slPair *tagList = raNextRecordAsSlPairList(lf);
if (tagList == NULL)
    return NULL;

/* Go through tag list, diverting some tags to be actual fields in patch. */
struct slPair *newTagList = NULL, *tag, *next;
struct raPatch *patch;
AllocVar(patch);
for (tag = tagList; tag != NULL; tag = next)
    {
    next = tag->next;
    if (sameString(tag->name, "db"))
	{
        patch->db = tag->val;
	freez(&tag);
	}
    else if (sameString(tag->name, clKey))
        {
	patch->track = tag->val;
	freez(&tag);
	}
    else if (sameString(tag->name, "file"))
        {
	patch->fileList = makeFileList(tag->val);
	freeMem(tag->val);
	freez(&tag);
	}
    else
        {
	slAddHead(&newTagList, tag);
	}
    }
slReverse(&newTagList);

if (patch->track == NULL)
    errAbort("Missing %s tag before line %d of %s", clKey, lf->lineIx, lf->fileName);
if (patch->fileList == NULL)
    errAbort("Missing %s tag before line %d of %s", "file", lf->lineIx, lf->fileName);
patch->tagList= newTagList;
return patch;
}

char *cloneFirstWord(char *string)
/* Clone leading word in string. */
{
char *s = skipLeadingSpaces(string);
char *e = skipToSpaces(s);
if (e == NULL)
    return cloneString(s);
int size = e - s;
if (size == 0)
    return NULL;
return cloneStringZ(s, size);
}

char *findLastMatchInString(char *string, char *match)
/* Return last position in string that starts with match.  Return NULL
 * if no match */
{
int matchLen = strlen(match);
int stringLen = strlen(string);
int startIx = stringLen - matchLen + 1;
while (--startIx >= 0)
    {
    if (memcmp(string + startIx, match, matchLen) == 0)
        return string + startIx;
    }
return NULL;
}

char *skipTrackDbPathPrefix(char *full)
/* Get suffix of full that skips the trackDb source code position. */
{
char *pat = "makeDb/trackDb/";
char *start = findLastMatchInString(full, pat);
if (start != NULL)
    start += strlen(pat);
return start;
}


void makeDirForFile(char *file)
/* Create directory that fill will sit in if it doesn't already exist. */
{
char dir[PATH_LEN], name[FILENAME_LEN], extension[FILEEXT_LEN];
splitPath(file, dir, name, extension);
if (dir[0] != 0)
    {
    char *simplePath = simplifyPathToDir(dir);
    if (simplePath[0] != 0)
	{
	uglyf("makeDirsOnPath(%s)\n", simplePath);
	makeDirsOnPath(simplePath);
	}
    freeMem(simplePath);
    }
}

void mustRename(char *oldPath, char *newPath)
/* Rename.  If fail print error message and abort. */
{
int err = rename(oldPath, newPath);
if (err != 0)
    errnoAbort("Couldn't rename %s to %s", oldPath, newPath);
}

struct slName *raNextStanza(struct lineFile *lf)
/* Return list of lines starting from current position, up through last line of next stanza.
 * May return a few blank/comment lines at end with no real stanza. */
{
struct slName *lineList = NULL;
char *line;
while (lineFileNext(lf, &line, NULL))
    {
    slNameAddHead(&lineList, line);
    line = skipLeadingSpaces(line);
    if (line[0] != 0 && line[0] != '#')
	break;
    }
while (lineFileNext(lf, &line, NULL))
    {
    char *s = skipLeadingSpaces(line);
    if (*s == 0)
         {
	 lineFileReuse(lf);
	 break;
	 }
    slNameAddHead(&lineList, line);
    }
slReverse(&lineList);
return lineList;
}

static void applyPatches(char *inName, struct slRef *patchRefList, char *keyField, char *outName)
/* Apply patches in list. */
{
int keyFieldLen = strlen(keyField);

/* Convert list of patch references to hash of patches. */
struct slRef *ref;
struct hash *patchHash = hashNew(0);
for (ref = patchRefList; ref != NULL; ref = ref->next)
    {
    struct raPatch *patch = ref->val;
    char *key = cloneFirstWord(patch->track);
    hashAdd(patchHash, key, patch);
    freeMem(key);
    }
uglyf("%d patches in hash, %d in list\n", patchHash->elCount, slCount(patchRefList));

/* Open input and output files. */
struct lineFile *lf = lineFileOpen(inName, TRUE);
FILE *f = mustOpen(outName, "w");

/* Scan through one stanza at a time. */
for (;;)
    {
    /* First just fetch stanza - including any blank and comment lines before, into a list of strings. */
    struct slName *stanza = raNextStanza(lf);
    if (stanza == NULL)
        break;

    /* Go through stanza once just to see if have any patches to apply. */
    struct raPatch *patch = NULL;
    struct slName *line;
    for (line = stanza; line != NULL; line = line->next)
        {
	char *tagStart = skipLeadingSpaces(line->name);
	if (startsWithWord(keyField, tagStart))
	     {
	     char *valStart = skipLeadingSpaces(tagStart + keyFieldLen);
	     char *key = cloneFirstWord(valStart);
	     patch = hashFindVal(patchHash, key);
	     freeMem(key);
	     break;
	     }
	}

    /* If have patch apply it, otherwise just copy. */
    if (patch)
        {
	uglyf("Got patch %s with %d tags starting %s %s\n", patch->track, slCount(patch->tagList), patch->tagList->name, (char *)patch->tagList->val);
	int indent = 0;
	struct hash *appliedHash = hashNew(0);
	for (line = stanza; line != NULL; line = line->next)
	    {
	    char *lineStart = line->name;
	    char *tagStart = skipLeadingSpaces(lineStart);
	    boolean copyLine = TRUE;
	    if (tagStart[0] != 0 && tagStart[0] != '#')
	        {
		indent = tagStart - lineStart;
		struct slPair *tagPatch;
		for (tagPatch = patch->tagList; tagPatch != NULL; tagPatch = tagPatch->next)
		    {
		    if (startsWithWord(tagPatch->name, tagStart))
		        {
			copyLine = FALSE;
			spaceOut(f, indent);
			fprintf(f, "%s %s\n", tagPatch->name, (char*)tagPatch->val);
			uglyf("Applying patch '%s' to modify %s'\n", (char*)tagPatch->val, tagStart);
			hashAdd(appliedHash, tagPatch->name, NULL);
			break;
			}
		    }
		}
	    if (copyLine)
	        {
		fprintf(f, "%s\n", line->name);
		}
	    }
	struct slPair *tagPatch;
	for (tagPatch = patch->tagList; tagPatch != NULL; tagPatch = tagPatch->next)
	    {
	    if (!hashLookup(appliedHash, tagPatch->name))
	        {
		spaceOut(f, indent);
		uglyf("Applying patch to %s adding %s %s\n", patch->track, tagPatch->name, (char*)tagPatch->val);
		fprintf(f, "%s %s\n", tagPatch->name, (char*)tagPatch->val);
		hashAdd(appliedHash, tagPatch->name, NULL);
		}
	    }
	hashFree(&appliedHash);
	}
    else
        {
	for (line = stanza; line != NULL; line = line->next)
	    {
	    fprintf(f, "%s\n", line->name);
	    }
	}

    slFreeList(&stanza);
    }
lineFileClose(&lf);
carefulClose(&f);
}

void trackDbPatch(char *patchesFile, char *backupDir)
/* trackDbPatch - Patch files in trackDb with a specification .ra file that has db, track, and file fields that say 
 * where to apply the patch, and other fields that are patched in.. */
{
/* Read in patch file. */
struct lineFile *lf = lineFileOpen(patchesFile, TRUE);
struct raPatch *patch, *patchList = NULL;
while ((patch = raPatchReadOne(lf)) != NULL)
    slAddHead(&patchList, patch);
slReverse(&patchList);

/* Group it by file to patch */
struct fileToPatch *file, *fileList = groupPatchesByFiles(patchList, clFirstFile);
verbose(1, "Got %d patches covering %d files in %s\n", slCount(patchList), slCount(fileList), patchesFile);

/* Do some setting up for backups and patches */
makeDirsOnPath(backupDir);
if (clPatchDir != NULL)
    makeDirsOnPath(clPatchDir);

for (file = fileList; file != NULL; file = file->next)
    {
    /* Figure out name that skips most of the long path to the trackDb files */
    char *relName = skipTrackDbPathPrefix(file->fileName);
    if (relName == NULL)
         relName = file->fileName;
    uglyf("full rel %s %s\n", file->fileName, relName);

    /* Create file names for backup file and temp file. */
    char backupPath[PATH_LEN], patchPath[PATH_LEN];
    char tempPath[PATH_LEN];
    safef(backupPath, sizeof(backupPath), "%s/%s", backupDir, relName);
    safef(tempPath, sizeof(tempPath), "%s/%s.tmp", backupDir, relName);
    if (clPatchDir)
        safef(patchPath, sizeof(patchPath), "%s/%s", clPatchDir, relName);
    else
        safef(patchPath, sizeof(patchPath), "%s", file->fileName);


    /* Do patch reading original source and creating temp file. */
    makeDirForFile(backupPath);
    applyPatches(file->fileName, file->patchList, clKey, tempPath);

    /* If testing, move temp to patch */
    if (clPatchDir)
	{
        makeDirForFile(patchPath);
	mustRename(tempPath, patchPath);
	}
    else
    /* If not testing then move original to backup and temp to original location. */
        {
	mustRename(file->fileName, backupPath);
	mustRename(tempPath, file->fileName);
	}
    }

lineFileClose(&lf);
}

int main(int argc, char *argv[])
/* Process command line. */
{
optionInit(&argc, argv, options);
if (argc != 3)
    usage();
clKey = optionVal("key", clKey);
clPatchDir = optionVal("test", clPatchDir);
clFirstFile = optionExists("firstFile");
if (!clPatchDir)
    errAbort("Must specify test option currently");
trackDbPatch(argv[1], argv[2]);
return 0;
}