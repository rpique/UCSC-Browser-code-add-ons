/* pngwrite.c - write out a memGfx to a PNG file, using the reference library libpng */
/* (libpng is available from sourceforge and included in many open source distros,
 *  has a wide-open license intended to encourage usage of the PNG format, and the lib
 *  has been under development and testing for 14 years -- http://libpng.org/) */

#ifdef USE_PNG

#include "common.h"
#include "memgfx.h"
#ifdef SETJMP_H
//.crap!;
#endif
#include "png.h"

static char const rcsid[] = "$Id: pngwrite.c,v 1.1 2009/08/05 23:34:31 angie Exp $";

static void pngAbort(png_structp png, png_const_charp errorMessage)
/* type png_error wrapper around errAbort */
{
errAbort((char *)errorMessage);
}

static void pngWarn(png_structp png, png_const_charp warningMessage)
/* type png_error wrapper around warn */
{
warn((char *)warningMessage);
}

boolean mgSaveToPng(FILE *png_file, struct memGfx *mg, boolean useAlpha)
/* Save PNG to an already open file.
 * If useAlpha, then the first color in memgfx's colormap/palette is
 * assumed to be the image background color, and pixels of that color
 * are made transparent. */
/* Reference: http://libpng.org/pub/png/libpng-1.2.5-manual.html */
{
if (!png_file || !mg)
    errAbort("mgSaveToPng: caller screwed up: png_file=[%lld], mg=[%lld]",
	     (long long int)png_file, (long long int)mg);
png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING,
					  NULL, // don't need pointer to data for err/warn handlers
					  pngAbort, pngWarn);
if (!png)
    {
    errAbort("png_write_struct failed");
    return FALSE;
    }
png_infop info = png_create_info_struct(png);
if (!info)
    {
    errAbort("png create_info_struct failed");
    png_destroy_write_struct(&png, NULL);
    return FALSE;
    }

// If setjmp returns nonzero, it means png_error is returning control here.
// But that should not happen because png_error should call pngAbort which calls errAbort.
if (setjmp(png_jmpbuf(png)))
    {
    png_destroy_write_struct(&png, &info);
    fclose(png_file);
    errAbort("pngwrite: setjmp nonzero.  "
	     "why didn't png_error..pngAbort..errAbort stop execution before this errAbort?");
    return FALSE;
    }

// Configure PNG output params:
png_init_io(png, png_file);
png_set_IHDR(png, info, mg->width, mg->height, 8, // 8=bit_depth
	     PNG_COLOR_TYPE_PALETTE, PNG_INTERLACE_NONE,
	     PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
png_set_PLTE(png, info,
	     (png_color *)(mg->colorMap), // png_color is same as struct rgbColor!
	     mg->colorsUsed);
if (useAlpha)
    {
    // First palette color is assumed to be background/transparent, so
    // that's the only one we need in the parallel array opacities[].
    png_byte opacities[] = {0};
    int num_opacities = ArraySize(opacities);
    png_color_16p nonPalette_opacities_values = NULL; // n/a for us, we're using palette
    png_set_tRNS(png, info, opacities, num_opacities, nonPalette_opacities_values);
    }

// Write header/params, write pixels, close and clean up.
// PNG wants a 2D array of pointers to byte offsets into palette/colorMap.
// mg has a 1D array of byte offsets.  Make row pointers for PNG:
png_byte **row_pointers = needMem(mg->height * sizeof(png_byte *));
int i;
for (i = 0;  i < mg->height;  i++)
    row_pointers[i] = &(mg->pixels[i*mg->width]);
png_set_rows(png, info, row_pointers);
png_write_png(png, info, PNG_TRANSFORM_IDENTITY, // no transform
	      NULL); // unused as of PNG 1.2
png_destroy_write_struct(&png, &info);
return TRUE;
}

void mgSavePng(struct memGfx *mg, char *filename, boolean useAlpha)
/* Save memory bitmap to filename as a PNG.
 * If useAlpha, then the first color in memgfx's colormap/palette is
 * assumed to be the image background color, and pixels of that color
 * are made transparent. */
{
FILE *pngFile = mustOpen(filename, "wb");
if (!mgSaveToPng(pngFile, mg, useAlpha))
    {
    remove(filename);
    errAbort("Couldn't save %s", filename);
    }
if (fclose(pngFile) != 0)
    errnoAbort("fclose failed");
}

#endif//def USE_PNG