/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#include <assert.h>
#include "extended/feature_node.h"
#include "extended/filter_stream.h"
#include "extended/filter_visitor.h"
#include "extended/genome_stream_rep.h"

struct FilterStream
{
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  GenomeVisitor *filter_visitor; /* the actual work is done in the visitor */
};

#define filter_stream_cast(GS)\
        genome_stream_cast(filter_stream_class(), GS);

static int filter_stream_next_tree(GenomeStream *gs, GtGenomeNode **gn,
                                   GtError *err)
{
  FilterStream *fs;
  int had_err;
  gt_error_check(err);
  fs = filter_stream_cast(gs);

  /* we still have nodes in the buffer */
  if (filter_visitor_node_buffer_size(fs->filter_visitor)) {
    /* return one of them */
    *gn = filter_visitor_get_node(fs->filter_visitor);
    return 0;
  }

  /* no nodes in the buffer -> get new nodes */
  while (!(had_err = genome_stream_next_tree(fs->in_stream, gn, err)) && *gn) {
    assert(*gn && !had_err);
    had_err = gt_genome_node_accept(*gn, fs->filter_visitor, err);
    if (had_err)
      break;
    if (filter_visitor_node_buffer_size(fs->filter_visitor)) {
      *gn = filter_visitor_get_node(fs->filter_visitor);
      return 0;
    }
  }

  /* either we have an error or no new node */
  assert(had_err || !*gn);
  return had_err;
}

static void filter_stream_free(GenomeStream *gs)
{
  FilterStream *fs = filter_stream_cast(gs);
  genome_visitor_delete(fs->filter_visitor);
  genome_stream_delete(fs->in_stream);
}

const GenomeStreamClass* filter_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (FilterStream),
                                         filter_stream_next_tree,
                                         filter_stream_free };
  return &gsc;
}

GenomeStream* filter_stream_new(GenomeStream *in_stream,
                                GtStr *seqid, GtStr *typefilter,
                                GtRange contain_range, GtRange overlap_range,
                                GtStrand strand, GtStrand targetstrand,
                                bool has_CDS, unsigned long max_gene_length,
                                unsigned long max_gene_num,
                                double min_gene_score, double max_gene_score,
                                double min_average_splice_site_prob,
                                unsigned long feature_num)
{
  GenomeStream *gs = genome_stream_create(filter_stream_class(),
                                          genome_stream_is_sorted(in_stream));
  FilterStream *filter_stream = filter_stream_cast(gs);
  assert(in_stream);
  filter_stream->in_stream = genome_stream_ref(in_stream);
  filter_stream->filter_visitor = filter_visitor_new(seqid, typefilter,
                                                     contain_range,
                                                     overlap_range, strand,
                                                     targetstrand, has_CDS,
                                                     max_gene_length,
                                                     max_gene_num,
                                                     min_gene_score,
                                                     max_gene_score,
                                                  min_average_splice_site_prob,
                                                     feature_num);
  return gs;
}
