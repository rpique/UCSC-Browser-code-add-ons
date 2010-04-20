/* mdbUpdate - Adds, updates or removes metadata objects and variables from the 'mdb' metadata table. */
#include "common.h"
#include "linefile.h"
#include "options.h"
#include "mdb.h"
#include "hash.h"

static char const rcsid[] = "$Id: mdbUpdate.c,v 1.2 2010/04/15 19:26:54 tdreszer Exp $";

#define OBJTYPE_DEFAULT "table"

void usage()
/* Explain usage and exit. */
{
errAbort(
  "mdbUpdate - Adds, updates or removes metadata objects and variables from the 'mdb' metadata table.\n"
  "usage:\n"
  "   mdbUpdate {db} [-table= [-force]] [-recreate] [-test [-verbose=2]]\n"
  "                  [-obj= [-type=] [-delete] [-var=] [-binary] [-val=]]\n"
  "                  [-vars=\"var1=val1 var2=val2... [-delete] [-var=] [-binary] [-val=]]\n"
  "                  [{fileName}] [-replace]\n"
  "Options:\n"
  "    {db}     Database to update metadata in.  This argument is required.\n"
  "    -table   Table to update metadata to.  Default is the sandbox version of\n"
  "             '" MDB_DEFAULT_NAME "'.\n"
  "       -recreate   Creates or empties the table.  No further arguemts are needed.\n"
  "       -force      Overrides restrictions placed on shared  '" MDB_DEFAULT_NAME "'.\n"
  "    -test    Does not update but only reports results.  Use with -verbose=2 to see SQL.\n"
  "  if {fileName} argument not provided, then -obj or -vars must be provided\n"
  "    -obj={objName}     Update this single object from the command lineas :\n"
  "    -vars={var=val...} Apply update to group of objects matching these restrictions.\n"
  "                       It is recommended you test the selection criteria with mdbPrint\n"
  "                       to verify which objects meet your selection criteria first.\n"
  "                       Use of 'var!=val', 'var=v%%' and 'var=?' are supported.\n"
  "    These options work on objects selected with -obj or -vars:\n"
  "       -delete         Remove a specific var or entire obj (if -var not provided).\n"
  "       -var={varName}  Provide variable name (if no -var then must be -delete)\n"
  "       -binary         NOT YET IMPLEMENTED.  This var has a binary val and -val={file}\n"
  "       -val={value}    (Enclosed in \"quotes if necessary\".) Not required to delete var\n"
  "       -setVars={var=val...}  Allows setting multiple var=val pairs.\n"
  "    [{fileName}] File containing formatted metadata lines.  (Ignored if -obj or -vars used).\n"
  "      -replace   Means remove all old variables for each object before adding new variables\n\n"
  "There are two ways to call mdbUpdate.  The object (or objects matching vars) and var to update "
  "can be declared on the command line, or a file of formatted metadata lines can be provided. "
  "The file can be the formatted output from mdbPrint or the following special formats:\n"
  "  metadata objName var1=val1 var2=\"val2 with spaces\" var3=...\n"
  "    Adds or updates the specific object and variables\n"
  "  metadata objName delete\n"
  "    delete all metadata for objName\n"
  "  metadata objName delete var=val var2=val2\n"
  "     deletes specifically named variables for objName. (Must provide val but value is ignored).\n"
  "  Special ENCODE format as produced by the doEncodeValidate.pl\n"
  "NOTE: Updates to the shared '" MDB_DEFAULT_NAME "' can only be done by a file written\n"
  "      directly from mdbPrint.  Update sandbox first, then move updates to shared table.\n"
  "HINT: Use '%%' in any command line obj, var or val as a wildcard for selection.\n\n"
  "Examples:\n"
  "  mdbUpdate hg19 -vars=\"grant=Snyder cell=GM12878 antibody=CTCF\" -var=expId -val=1427\n"
  "            Update all objects matcing Snyder/GM12878/CTCF and set the expId=1472.\n"
  "  mdbUpdate hg19 -obj=fredsTable -var=description val=\"Ethyl's husband's clutter\"\n"
  "            Updates fredsTable with a description.\n"
  "  mdbUpdate mm9 -table=mdb_braney -recreate\n"
  "            Creates or empties the named metadata table.\n"
  "  mdbUpdate hg18 -test vars=\"composite=wgEncodeDukeDNase\" -delete\n"
  "            Tests the delete of all objects that have the named composite defined.\n"
  "  mdbUpdate hg18 -table=mdb encBroadHistone.mdb.ra -replace\n"
  "            Replaces all metadata for objects found in the provided file.\n"
  "            File must have been printed with mdbPrint.\n"
  );
}

static struct optionSpec optionSpecs[] = {
    {"table",   OPTION_STRING}, // default "mdb"
    {"obj",     OPTION_STRING}, // objName or objId
    {"vars",    OPTION_STRING}, // Select set of object by vars
    {"var",     OPTION_STRING}, // variable
    {"val",     OPTION_STRING}, // value
    {"setVars", OPTION_STRING}, // Allows setting multiple var=val pairs
    {"delete",  OPTION_BOOLEAN},// delete one obj or obj/var
    {"binary",  OPTION_BOOLEAN},// val is binary (NOT YET IMPLEMENTED) implies -val={file}
    {"replace", OPTION_BOOLEAN},// replace entire obj when loading from file
    {"recreate",OPTION_BOOLEAN},// creates or recreates the table
    {"force",   OPTION_BOOLEAN},// override restrictions on shared table
    {"test",    OPTION_BOOLEAN},// give it a test, will ya?
    {NULL,      0}
};

int main(int argc, char *argv[])
// Process command line.
{
if(argc == 1)
    usage();

struct mdbObj * mdbObjs = NULL;

optionInit(&argc, argv, optionSpecs);

if(argc < 2)
    {
    verbose(1, "REQUIRED 'DB' argument not found:\n");
    usage();
    }

char *db         = argv[1];
char *table      = optionVal("table",NULL);
boolean deleteIt = optionExists("delete");
boolean testIt   = optionExists("test");
boolean recreate = optionExists("recreate");
boolean force    = optionExists("force");
boolean replace  = FALSE;
char *var        = optionVal("var",NULL);
char *val        = optionVal("val",NULL);
char *setVars    = optionVal("setVars",NULL);

struct sqlConnection *conn = sqlConnect(db);

// Find the table if necessary
if(table == NULL)
    {
    table = mdbTableName((recreate?NULL:conn),TRUE); // Look for sandBox name first
    if(table == NULL)
        {
        table = mdbTableName((recreate?NULL:conn),FALSE); // Okay, default then
        if(table == NULL)  // Now you are just getting me angry!
            {
            sqlDisconnect(&conn);
            if(!recreate) // assertable
                errAbort("No '%s.%s' found.  Consider using -recreate flag.\n",db,MDB_DEFAULT_NAME);
            else
                errAbort("No '%s.%s' found.\n",db,MDB_DEFAULT_NAME);
            }
        }
    verbose(1, "Using table named '%s.%s'.\n",db,table);
    }

boolean sharedTbl = sameWord(table,MDB_DEFAULT_NAME);  // Special restrictions apply

// Recreate the table
if(recreate)
    {
    if(sharedTbl && ! force)
        {
        sqlDisconnect(&conn);
        errAbort("NOT SUPPORTED for shared table '%s'.\n",MDB_DEFAULT_NAME);
        }
    boolean recreated = sqlTableExists(conn,table);
    mdbReCreate(conn,table,testIt);
    if(testIt)
        {
        verbose(1, "Would %screate table named '%s'.\n",
                (recreated?"re":""),table);
        if(!recreated)
            {
        sqlDisconnect(&conn);
        if(optionExists("obj") || optionExists("vars") || argc > 2)
                verbose(1, "Can't test further commands.  Consider '-db= [-table=] -recreate' as the only arguments.\n");
            return 0;  // Don't test any update if we haven't actually created the table!
            }
        }
    else
        verbose(1, "%s table named '%s'.\n",(recreated?"Recreated":"Created"),table);
    }

if(argc > 2 && (deleteIt || var != NULL || val != NULL || setVars != NULL))
    {
    verbose(1, "INCONSISTENT REQUEST: can't combine supplied file with -delete, -var, -val or -setVars.\n");
    usage();
    }
if(deleteIt && var != NULL && val != NULL)
    {
    verbose(1, "INCONSISTENT REQUEST: can't combine -delete with -var and -val.\n");
    usage();
    }
if (argc != 3 && !deleteIt)
    {
    if(setVars == NULL && (var == NULL || val == NULL))
        {
        if(recreate) // no problem
            return 0;
        verbose(1, "INCONSISTENT REQUEST: need both -var and -val.\n");
        usage();
        }
    else if (setVars != NULL && (var != NULL || val != NULL))
        {
        if(recreate) // no problem
            return 0;
        verbose(1, "INCONSISTENT REQUEST: can't combin -var or -val with -setVars.\n");
        usage();
        }
    }

// Now get the object list
if(optionExists("obj"))
    {
    if(sharedTbl && !force)
        {
        sqlDisconnect(&conn);
        verbose(1, "NOT SUPPORTED for shared table '%s'.\n",MDB_DEFAULT_NAME);
        usage(); // Must not have submitted formatted file also
        }
    if(argc > 2 || optionExists("vars"))
        {
        sqlDisconnect(&conn);
        verbose(1, "INCONSISTENT REQUEST: can't combine -obj with -vars or a supplied file.\n");
        usage(); // Must not have submitted formatted file also
        }

    mdbObjs = mdbObjCreate(optionVal("obj",  NULL),var,
                            (optionExists("binary") ? "binary" : "txt"), // FIXME: don't know how to deal with binary yet
                            val);
    mdbObjs->deleteThis = deleteIt;

    if(setVars != NULL)
        mdbObjSwapVars(mdbObjs,setVars,deleteIt);

    verbose(2, "metadata %s %s%s%s%s\n",
         mdbObjs->obj,(mdbObjs->deleteThis ? "delete ":""),
        (mdbObjs->vars && mdbObjs->vars->var!=NULL?mdbObjs->vars->var:""),
        (mdbObjs->vars && mdbObjs->vars->val!=NULL?"=":""),
        (mdbObjs->vars && mdbObjs->vars->val!=NULL?mdbObjs->vars->val:""));
    }
else if(optionExists("vars"))
    {
    if(sharedTbl && !force)
        {
        sqlDisconnect(&conn);
        verbose(1, "NOT SUPPORTED for shared table '%s'.\n",MDB_DEFAULT_NAME);
        usage(); // Must not have submitted formatted file also
        }
    if(argc > 2)
        {
        sqlDisconnect(&conn);
        verbose(1, "INCONSISTENT REQUEST: can't combine -vars with a supplied file.\n");
        usage(); // Must not have submitted formatted file also
        }
    struct mdbByVar * mdbByVars = mdbByVarsLineParse(optionVal("vars", NULL));
    mdbObjs = mdbObjsQueryByVars(conn,table,mdbByVars);

    // replace all found vars but update request
    if(setVars != NULL)
        mdbObjSwapVars(mdbObjs,setVars,deleteIt);
    else
        mdbObjTransformToUpdate(mdbObjs,var,
                            (optionExists("binary") ? "binary" : "txt"), // FIXME: don't know how to deal with binary yet
                             val,deleteIt);
    }
else // Must be submitting formatted file
    {
    if(argc != 3)
        {
        sqlDisconnect(&conn);
        if(recreate) // no problem
            return 0;
        verbose(1, "REQUIRED: must declare -obj, -vars or supply a file.\n");
        usage(); // Must not have submitted formatted file also
        }

    replace = optionExists("replace");
    boolean validated = FALSE;
    mdbObjs = mdbObjsLoadFromFormattedFile(argv[2],&validated);
    if(sharedTbl && !force && !validated)
        {
        sqlDisconnect(&conn);
        errAbort("Update to shared table '%s' requires file directly written by mdbPrint from sandbox file.\n", table);
        }
    if(mdbObjs != NULL)
        verbose(1, "Read %d metadata objects from %s\n", slCount(mdbObjs),argv[1]);
    }

int count = 0;

if(mdbObjs != NULL)
    {
    if(testIt && verboseLevel() > 2)
        mdbObjPrint(mdbObjs,FALSE);

    count = mdbObjsSetToDb(conn,table,mdbObjs,replace,testIt);
    }
if(testIt)
    verbose(1, "Command would affected %d row(s) in %s.%s\n", count,db,table);
else
    verbose(1, "Affected %d row(s) in %s.%s\n", count,db,table);

sqlDisconnect(&conn);
mdbObjsFree(&mdbObjs);
return 0;

// TODO:
// 1) Case insensitive hashs?
// 2) -ra by default (-1 line)?
// 3) expId table?  -exp=start requests generating unique ids for selected vars, then updating them. -expTbl generates expTable as id,"var=val var=val var=val"
}