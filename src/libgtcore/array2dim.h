/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef ARRAY2DIM_H
#define ARRAY2DIM_H

#include <libgtcore/xansi.h>

#define array2dim_malloc(ARRAY2DIM, ROWS, COLUMNS, TYPE, ENV)                 \
        {                                                                     \
          unsigned long a2d_i;                                                \
          ARRAY2DIM = env_ma_malloc(ENV, sizeof (TYPE*) * (ROWS));            \
          (ARRAY2DIM)[0] = env_ma_malloc(ENV,                                 \
                                         sizeof (TYPE) * (ROWS) * (COLUMNS)); \
          for (a2d_i = 1; a2d_i < (ROWS); a2d_i++)                            \
            (ARRAY2DIM)[a2d_i] = (ARRAY2DIM)[a2d_i-1] + (COLUMNS);            \
        }

#define array2dim_calloc(ARRAY2DIM, ROWS, COLUMNS, TYPE, ENV)                 \
        {                                                                     \
          unsigned long a2d_i;                                                \
          ARRAY2DIM = env_ma_malloc(ENV, sizeof (TYPE*) * (ROWS));            \
          (ARRAY2DIM)[0] = env_ma_calloc(ENV, (ROWS) * (COLUMNS),             \
                                         sizeof (TYPE));                      \
          for (a2d_i = 1; a2d_i < (ROWS); a2d_i++)                            \
            (ARRAY2DIM)[a2d_i] = (ARRAY2DIM)[a2d_i-1] + (COLUMNS);            \
        }

#define array2dim_delete(ARRAY2DIM, ENV)                                      \
        env_ma_free((ARRAY2DIM)[0], ENV);                                     \
        env_ma_free(ARRAY2DIM, ENV);

#if 0
  example usage:

  double **a2dim;

  /* create a 10 * 20 double array */
  array2dim_malloc(a2dim, 10, 20, double, env);
  /* ... (use array a2dim in conventional way via a2dim[row][column]) */
  array2dim_delete(a2dim, env);
#endif

#endif
