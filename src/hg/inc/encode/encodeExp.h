/* encodeExp.h was originally generated by the autoSql program, which also
 * generated encodeExp.c and encodeExp.sql.  This header links the database and
 * the RAM representation of objects. */

#ifndef ENCODEEXP_H
#define ENCODEEXP_H

#include "jksql.h"

#define ENCODEEXP_NUM_COLS 8

struct encodeExp
/* ENCODE experiments */
    {
    struct encodeExp *next;  /* Next in singly linked list. */
    int ix;	/* auto-increment ID */
    char *organism;	/* human | mouse */
    char *accession;	/* wgEncodeE[H|M]00000N */
    char *lab;	/* lab name from ENCODE cv.ra */
    char *dataType;	/* dataType from ENCODE cv.ra */
    char *cellType;	/* cellType from ENCODE cv.ra */
    char *expVars;	/* typeOfTerm=term list of experiment-defining variables */
    char *lastUpdated;  /* auto-update timestamp */  // WARNING: hand-edit here
    };

void encodeExpStaticLoad(char **row, struct encodeExp *ret);
/* Load a row from encodeExp table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct encodeExp *encodeExpLoadByQuery(struct sqlConnection *conn, char *query);
/* Load all encodeExp from table that satisfy the query given.
 * Where query is of the form 'select * from example where something=something'
 * or 'select example.* from example, anotherTable where example.something =
 * anotherTable.something'.
 * Dispose of this with encodeExpFreeList(). */

void encodeExpSaveToDb(struct sqlConnection *conn, struct encodeExp *el, char *tableName, int updateSize);
/* Save encodeExp as a row to the table specified by tableName.
 * As blob fields may be arbitrary size updateSize specifies the approx size
 * of a string that would contain the entire query. Arrays of native types are
 * converted to comma separated strings and loaded as such, User defined types are
 * inserted as NULL. Note that strings must be escaped to allow insertion into the database.
 * For example "autosql's features include" --> "autosql\'s features include"
 * If worried about this use encodeExpSaveToDbEscaped() */

void encodeExpSaveToDbEscaped(struct sqlConnection *conn, struct encodeExp *el, char *tableName, int updateSize);
/* Save encodeExp as a row to the table specified by tableName.
 * As blob fields may be arbitrary size updateSize specifies the approx size.
 * of a string that would contain the entire query. Automatically
 * escapes all simple strings (not arrays of string) but may be slower than encodeExpSaveToDb().
 * For example automatically copies and converts:
 * "autosql's features include" --> "autosql\'s features include"
 * before inserting into database. */

struct encodeExp *encodeExpLoad(char **row);
/* Load a encodeExp from row fetched with select * from encodeExp
 * from database.  Dispose of this with encodeExpFree(). */

struct encodeExp *encodeExpLoadAll(char *fileName);
/* Load all encodeExp from whitespace-separated file.
 * Dispose of this with encodeExpFreeList(). */

struct encodeExp *encodeExpLoadAllByChar(char *fileName, char chopper);
/* Load all encodeExp from chopper separated file.
 * Dispose of this with encodeExpFreeList(). */

#define encodeExpLoadAllByTab(a) encodeExpLoadAllByChar(a, '\t');
/* Load all encodeExp from tab separated file.
 * Dispose of this with encodeExpFreeList(). */

struct encodeExp *encodeExpCommaIn(char **pS, struct encodeExp *ret);
/* Create a encodeExp out of a comma separated string.
 * This will fill in ret if non-null, otherwise will
 * return a new encodeExp */

void encodeExpFree(struct encodeExp **pEl);
/* Free a single dynamically allocated encodeExp such as created
 * with encodeExpLoad(). */

void encodeExpFreeList(struct encodeExp **pList);
/* Free a list of dynamically allocated encodeExp's */

void encodeExpOutput(struct encodeExp *el, FILE *f, char sep, char lastSep);
/* Print out encodeExp.  Separate fields with sep. Follow last field with lastSep. */

#define encodeExpTabOut(el,f) encodeExpOutput(el,f,'\t','\n');
/* Print out encodeExp as a line in a tab-separated file. */

#define encodeExpCommaOut(el,f) encodeExpOutput(el,f,',',',');
/* Print out encodeExp as a comma separated list including final comma. */

void encodeExpJsonOutput(struct encodeExp *el, FILE *f);
/* Print out encodeExp in JSON format. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

#include "mdb.h"

/* WARNING:  When schema is changed, encodeExp.c will need changes as well */

#define ENCODE_EXP_FIELD_IX             "ix"
#define ENCODE_EXP_FIELD_ORGANISM       "organism"
#define ENCODE_EXP_FIELD_ACCESSION      "accession"
#define ENCODE_EXP_FIELD_LAB            "lab"
#define ENCODE_EXP_FIELD_DATA_TYPE      "dataType"
#define ENCODE_EXP_FIELD_CELL_TYPE      "cellType"
#define ENCODE_EXP_FIELD_FACTORS        "expVars"
#define ENCODE_EXP_FIELD_LAST_UPDATED   "lastUpdated"

#define ENCODE_EXP_NO_CELL              "None"
#define ENCODE_EXP_NO_VAR               "None"

#define ENCODE_EXP_TABLE        "encodeExp"
#define ENCODE_EXP_DATABASE     "hgFixed"
#define ENCODE_EXP_ACC_PREFIX   "wgEncodeE"
#define ENCODE_EXP_TABLE_LOCK   "lock_encodeExp"

void encodeExpFieldIndex(char *fieldName);
/* Get column number of named field in EncodeExp schema */

void encodeExpTableCreate(struct sqlConnection *conn, char *table);
/* Create an encodeExp table */

struct encodeExp *encodeExpLoadAllFromTable(struct sqlConnection *conn, char *table);
/* Load all encodeExp in table */

struct encodeExp *encodeExpFromMdb(struct sqlConnection *conn,char *db, struct mdbObj *mdb);
/* Create an encodeExp from an ENCODE metaDb object */

struct encodeExp *encodeExpFromMdbVars(char *db, struct mdbVar *vars);
// Creates and returns an encodeExp struct from mdbVars, but does not .touch the table
// Only Experiment Defining Variables should be in the list.

struct encodeExp *encodeExpFromRa(struct hash *ra);
/* Load an encodeExp from a Ra hash */

struct hash *encodeExpToRaFile(struct encodeExp *exp, FILE *f);
/* Create a Ra hash from an encodeExp.  Print to file if non NULL */

struct hash *encodeExpToRa(struct encodeExp *exp);
/* Create a Ra hash from an encodeExp */

struct encodeExp *encodeExpGetByIdFromTable(struct sqlConnection *conn, char *tableName, int id);
/* Return experiment specified by id from named table */

struct encodeExp *encodeExpGetById(struct sqlConnection *conn, int id);
/* Return experiment specified by id from default table */

void encodeExpAdd(struct sqlConnection *conn, char *tableName, struct encodeExp *exp);
/* Add encodeExp as a new row to the table specified by tableName.
   Update accession using index assigned with autoincrement */

char *encodeExpAddAccession(struct sqlConnection *conn, char *tableName, int id);
/* Add accession field to an existing "temp" experiment.  This is done
 * after experiment is determined to be valid.
 * Return the accession. */

void encodeExpRemoveAccession(struct sqlConnection *conn, char *tableName, int id);
/* Revoke an experiment by removing the accession. */

char *encodeExpKey(struct encodeExp *exp);
/* Create a hash key from an encodeExp */

char *encodeExpVars(struct encodeExp *exp);
// Create a string of all experiment defining vars and vals as "lab=UW dataType=ChipSeq ..."
// WARNING: May be missing var=None if the var was added after composite had exps.

struct encodeExp *encodeExpGetFromTable(char *organism, char *lab, char *dataType, char *cell,
                                struct slPair *varPairs, char *table);
/* Return experiments matching args in named experiment table.
 * Organism, Lab and DataType must be non-null */

struct encodeExp *encodeExpGet(char *organism, char *lab, char *dataType, char *cell,
                                struct slPair *varPairs);
/* Return experiments matching args in default experiment table.
 * Organism, Lab and DataType must be non-null */

struct encodeExp *encodeExpGetByMdbVarsFromTable(char *db, struct mdbVar *vars, char *table);
/* Return experiments by looking up mdb var list from the named experiment table */

struct encodeExp *encodeExpGetByMdbVars(char *db, struct mdbVar *vars);
/* Return experiments by looking up mdb var list from the default experiment table */

struct encodeExp *encodeExpGetOrCreateByMdbVarsFromTable(char *db, struct mdbVar *vars, char *table);
// Return experiment looked up or created from the mdb var list from the named experiment table.
#define encodeExpGetOrCreateByMdbVars(db,vars) encodeExpGetOrCreateByMdbVarsFromTable((db),(vars),ENCODE_EXP_TABLE)

int encodeExpExists(char *db, struct mdbVar *vars);
/* Return TRUE if at least one experiment exists for these vars */

char *encodeGetExpAccessionByMdbVars(char *db, struct mdbVar *vars);
/* Return accession of (first) experiment matching vars, or NULL if not found */

void encodeExpUpdateField(struct sqlConnection *conn, char *tableName,
                                int id, char *field, char *val);
/* Update field in encodeExp identified by id with value.
   Only supported for a few non-interdependent fields
   and only for non-accessioned experiments */

void encodeExpUpdateExpVars(struct sqlConnection *conn, char *tableName,
                                char *accession, struct slPair *varPairs);
/* Update expVars in encodeExp identified by accession */

#endif /* ENCODEEXP_H */


