/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef SPLICE_SITE_INFO_STREAM_H
#define SPLICE_SITE_INFO_STREAM_H

#include <stdio.h>
#include "extended/node_stream.h"
#include "extended/region_mapping.h"

/* implements the ``genome_stream'' interface */
typedef struct SpliceSiteInfoStream SpliceSiteInfoStream;

const GtNodeStreamClass* splice_site_info_stream_class(void);

/* create a SpliceSiteInfoStream, takes ownership of GtRegionMapping  */
GtNodeStream*            splice_site_info_stream_new(GtNodeStream*,
                                                     GtRegionMapping*);
/* returns if an intron has been processed, false otherwise */
bool                     splice_site_info_stream_show(GtNodeStream*);

#endif
