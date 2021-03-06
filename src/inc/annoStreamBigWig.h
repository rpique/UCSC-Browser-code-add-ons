/* annoStreamBigWig -- subclass of annoStreamer for bigWig file or URL */

#ifndef ANNOSTREAMBIGWIG_H
#define ANNOSTREAMBIGWIG_H

#include "annoStreamer.h"

struct annoStreamer *annoStreamBigWigNew(char *fileOrUrl);
/* Create an annoStreamer (subclass) object from a file or URL. */

#endif//ndef ANNOSTREAMBIGWIG_H
