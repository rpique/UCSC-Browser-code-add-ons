/* rqlEval - evaluate tree returned by rqlParse given a record and function to lookup fields
 * in the record. . */

#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "dystring.h"
#include "tokenizer.h"
#include "sqlNum.h"
#include "rql.h"

static char const rcsid[] = "$Id: rqlEval.c,v 1.1 2009/12/02 19:09:53 kent Exp $";

struct rqlEval rqlEvalCoerceToBoolean(struct rqlEval r)
/* Return TRUE if it's a nonempty string or a non-zero number. */
{
switch (r.type)
    {
    case rqlTypeBoolean:
	break;	/* It's already done. */
    case rqlTypeString:
        r.val.b = (r.val.s != NULL && r.val.s[0] != 0);
	break;
    case rqlTypeInt:
        r.val.b = (r.val.i != 0);
	break;
    case rqlTypeDouble:
        r.val.b = (r.val.x != 0.0);
	break;
    default:
	internalErr();
	r.val.b = FALSE;
	break;
    }
r.type = rqlTypeBoolean;
return r;
}

static struct rqlEval rqlEvalEq(struct rqlParse *p, void *record, RqlEvalLookup lookup)
/* Return true if two children are equal regardless of children type
 * (which are just gauranteed to be the same). */
{
struct rqlParse *lp = p->children;
struct rqlParse *rp = lp->next;
struct rqlEval lv = rqlEvalOnRecord(lp, record, lookup);
struct rqlEval rv = rqlEvalOnRecord(rp, record, lookup);
struct rqlEval res;
res.type = rqlTypeBoolean;
assert(lv.type == rv.type);
switch (lv.type)
    {
    case rqlTypeBoolean:
        res.val.b = (lv.val.b == rv.val.b);
	break;
    case rqlTypeString:
	res.val.b = sameString(lv.val.s, rv.val.s);
	break;
    case rqlTypeInt:
	res.val.b = (lv.val.i == rv.val.i);
	break;
    case rqlTypeDouble:
	res.val.b = (lv.val.x == rv.val.x);
	break;
    default:
	internalErr();
	res.val.b = FALSE;
	break;
    }
return res;
}

static struct rqlEval rqlEvalLt(struct rqlParse *p, void *record, RqlEvalLookup lookup)
/* Return true if r < l . */
{
struct rqlParse *lp = p->children;
struct rqlParse *rp = lp->next;
struct rqlEval lv = rqlEvalOnRecord(lp, record, lookup);
struct rqlEval rv = rqlEvalOnRecord(rp, record, lookup);
struct rqlEval res;
res.type = rqlTypeBoolean;
switch (lv.type)
    {
    case rqlTypeBoolean:
        res.val.b = (lv.val.b < rv.val.b);
	break;
    case rqlTypeString:
	res.val.b = strcmp(lv.val.s, rv.val.s) < 0;
	break;
    case rqlTypeInt:
	res.val.b = (lv.val.i < rv.val.i);
	break;
    case rqlTypeDouble:
	res.val.b = (lv.val.x < rv.val.x);
	break;
    default:
	internalErr();
	res.val.b = FALSE;
	break;
    }
return res;
}

static struct rqlEval rqlEvalGt(struct rqlParse *p, void *record, RqlEvalLookup lookup)
/* Return true if r > l . */
{
struct rqlParse *lp = p->children;
struct rqlParse *rp = lp->next;
struct rqlEval lv = rqlEvalOnRecord(lp, record, lookup);
struct rqlEval rv = rqlEvalOnRecord(rp, record, lookup);
struct rqlEval res;
res.type = rqlTypeBoolean;
switch (lv.type)
    {
    case rqlTypeBoolean:
        res.val.b = (lv.val.b > rv.val.b);
	break;
    case rqlTypeString:
	res.val.b = strcmp(lv.val.s, rv.val.s) > 0;
	break;
    case rqlTypeInt:
	res.val.b = (lv.val.i > rv.val.i);
	break;
    case rqlTypeDouble:
	res.val.b = (lv.val.x > rv.val.x);
	break;
    default:
	internalErr();
	res.val.b = FALSE;
	break;
    }
return res;
}

static struct rqlEval rqlEvalLike(struct rqlParse *p, void *record, RqlEvalLookup lookup)
/* Return true if r like l . */
{
struct rqlParse *lp = p->children;
struct rqlParse *rp = lp->next;
struct rqlEval lv = rqlEvalOnRecord(lp, record, lookup);
struct rqlEval rv = rqlEvalOnRecord(rp, record, lookup);
struct rqlEval res;
res.type = rqlTypeBoolean;
assert(rv.type == rqlTypeString);
assert(rv.type == lv.type);
res.val.b = sqlMatchLike(rv.val.s, lv.val.s);
return res;
}

struct rqlEval rqlEvalOnRecord(struct rqlParse *p, void *record, RqlEvalLookup lookup)
/* Evaluate self on ra. */
{
struct rqlEval res;
switch (p->op)
    {
    case rqlOpLiteral:
	res.val = p->val;
	res.type = p->type;
	break;
    case rqlOpSymbol:
	res.type = rqlTypeString;
	char *s = lookup(record, p->val.s);
	if (s == NULL)
	    res.val.s = "";
	else
	    res.val.s = s;
	break;
    case rqlOpEq:
	res = rqlEvalEq(p, record, lookup);
	break;
    case rqlOpNe:
	res = rqlEvalEq(p, record, lookup);
	res.val.b = !res.val.b;
	break;

    /* Inequalities. */
    case rqlOpLt:
        res = rqlEvalLt(p, record, lookup);
	break;
    case rqlOpGt:
        res = rqlEvalGt(p, record, lookup);
	break;
    case rqlOpLe:
        res = rqlEvalGt(p, record, lookup);
	res.val.b = !res.val.b;
	break;
    case rqlOpGe:
        res = rqlEvalLt(p, record, lookup);
	res.val.b = !res.val.b;
	break;
    case rqlOpLike:
        res = rqlEvalLike(p,record, lookup);
	break;

    /* Logical ops. */
    case rqlOpAnd:
	{
        res.type = rqlTypeBoolean;
	res.val.b = TRUE;
	struct rqlParse *c;
	for (c = p->children; c != NULL; c= c->next)
	    {
	    struct rqlEval e = rqlEvalOnRecord(c, record, lookup);
	    if (!e.val.b)
		{
	        res.val.b = FALSE;
		break;
		}
	    }
	break;
	}
    case rqlOpOr:
	{
        res.type = rqlTypeBoolean;
	res.val.b = FALSE;
	struct rqlParse *c;
	for (c = p->children; c != NULL; c= c->next)
	    {
	    struct rqlEval e = rqlEvalOnRecord(c, record, lookup);
	    if (e.val.b)
		{
	        res.val.b = TRUE;
		break;
		}
	    }
	break;
	}

    case rqlOpNot:
        res = rqlEvalOnRecord(p->children, record, lookup);
	res.val.b = !res.val.b;
	break;

    /* Type casts. */
    case rqlOpStringToBoolean:
	res = rqlEvalOnRecord(p->children, record, lookup);
	res.type = rqlTypeBoolean;
	res.val.b = (res.val.s[0] != 0);
	break;
    case rqlOpIntToBoolean:
        res = rqlEvalOnRecord(p->children, record, lookup);
	res.type = rqlTypeBoolean;
	res.val.b = (res.val.i != 0);
	break;
    case rqlOpDoubleToBoolean:
        res = rqlEvalOnRecord(p->children, record, lookup);
	res.type = rqlTypeBoolean;
	res.val.b = (res.val.x != 0.0);
	break;
    case rqlOpStringToInt:
	res = rqlEvalOnRecord(p->children, record, lookup);
	res.type = rqlTypeInt;
	res.val.i = atoi(res.val.s);
	break;
    case rqlOpStringToDouble:
	res = rqlEvalOnRecord(p->children, record, lookup);
	res.type = rqlTypeDouble;
	res.val.x = atof(res.val.s);
	break;
    case rqlOpBooleanToInt:
	res = rqlEvalOnRecord(p->children, record, lookup);
	res.type = rqlTypeInt;
	res.val.i = res.val.b;
	break;
    case rqlOpBooleanToDouble:
	res = rqlEvalOnRecord(p->children, record, lookup);
	res.type = rqlTypeDouble;
	res.val.x = res.val.b;
	break;
    case rqlOpIntToDouble:
	res = rqlEvalOnRecord(p->children, record, lookup);
	res.type = rqlTypeDouble;
	res.val.x = res.val.b;
	break;

    case rqlOpUnaryMinusDouble:
        res = rqlEvalOnRecord(p->children, record, lookup);
	res.type = rqlTypeDouble;
	res.val.x = -res.val.x;
	break;

    default:
        errAbort("Unknown op %s\n", rqlOpToString(p->op));
	res.type = rqlTypeInt;	// Keep compiler from complaining.
	res.val.i = 0;	// Keep compiler from complaining.
	break;
    }
return res;
}
