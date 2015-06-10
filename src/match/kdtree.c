/*
  Copyright (c) 2015 Jörg Winkler <joerg.winkler@studium.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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

// #include <stdbool.h>
#include <limits.h>
#include "match/kdtree.h"
#include "core/arraydef.h"
#include "core/ma.h"

#define UNDEF ULONG_MAX
GT_DECLAREARRAYSTRUCT(GtKdtreeNode);

struct GtKdtree
{
  GtUword dimension;
  GtArrayGtKdtreeNode *nodes;
  int (* cmp)(GtKdtreeNode *, GtKdtreeNode *);
};

struct GtKdtreeNode {
  GtUword id_self;
  GtUword id_low;
  GtUword id_high;
  void *keys;
  void *value;
};

GtKdtree *gt_kdtree_new(void)
{
  GtKdtree *kdtree;
  kdtree = gt_malloc(sizeof (GtKdtree));
  kdtree->dimension = 0;
  GT_INITARRAY(kdtree->nodes, GtKdtreeNode);
  return kdtree;
}

GtKdtreeNode *gt_kdtreenode_new()
{
  GtKdtreeNode *kdtreenode;
  kdtreenode = gt_malloc(sizeof (GtKdtreeNode));
  kdtreenode->id_self = 0;
  kdtreenode->id_high = UNDEF;
  kdtreenode->id_low = UNDEF;
}

void gt_kdtreenode_delete(GtKdtreeNode *kdtreenode)
{
  if (kdtreenode != NULL) {
    gt_free(kdtreenode->keys);
    gt_free(kdtreenode->value);
  }
}

void gt_kdtree_delete(GtKdtree *kdtree)
{
  GT_FREEARRAY(kdtree->nodes, GtKdtreeNode);
  gt_free(kdtree);
}

GtUword gt_kdtree_size(const GtKdtree *kdtree)
{
  return kdtree->nodes->nextfreeGtKdtreeNode;
}

GtKdtreeNode *gt_kdtree_get(const GtKdtree *kdtree, GtUword id)
{
  if (id < gt_kdtree_size(kdtree)) {
    return kdtree->nodes->spaceGtKdtreeNode + id;
  } else {
    return NULL;
  }
}

GtUword gt_kdtree_insert(GtKdtree *kdtree, GtKdtreeNode *new, GtUword currdim)
{
  gt_assert(kdtree && new);
  GtKdtreeNode *currnode = gt_kdtree_get(kdtree, new->id);
  if (currnode == NULL) {
    GT_STOREINARRAY(kdtree->nodes, GtKdtreeNode, ARR_INCR, new);
  } else {
    int cmpval = kdtree->cmp(new, currnode);
    if (cmpval == 0) {
      *(currnode->value) = *(new->value);
      gt_kdtreenode_delete(new);
    } else if (cmpval > 0) {
      if (currnode->id_high == UNDEF) {
        new->id = gt_kdtree_size(kdtree);
        GT_STOREINARRAY(kdtree->nodes, GtKdtreeNode, ARR_INCR, new);
      } else {
        new->id_self = currnode->id_high;
        gt_kdtree_insert(kdtree, new, (currdim+1)%kdtree->dimension);
      }
    } else {
      if (currnode->id_low == UNDEF) {
        new->id = gt_kdtree_size(kdtree);
        GT_STOREINARRAY(kdtree->nodes, GtKdtreeNode, ARR_INCR, new);
      } else {
        new->id_self = currnode->id_low;
        gt_kdtree_insert(kdtree, new, (currdim+1)%kdtree->dimension);
      }
    }
  }
  return new->id_self;

}
