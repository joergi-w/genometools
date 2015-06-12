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

#include <limits.h>
#include "match/kdtree.h"
#include "core/arraydef.h"
#include "core/ma.h"

#define UNDEF ULONG_MAX
#define ARR_INCR 256
GT_DECLAREARRAYSTRUCT(GtKdtreeNode);

struct GtKdtree
{
  GtUword dimension;
  GtArrayGtKdtreeNode *nodes;
  int (* cmp)(const void *, const void *, GtUword);
};

struct GtKdtreeNode {
  GtUword id_self;
  GtUword id_low;
  GtUword id_high;
  void *key;
  void *value;
};

/* Constructors and Destructors for GtKdtree and GtKdtreeNode */

GtKdtreeNode *gt_kdtreenode_new(const void *key, const void *value)
{
  GtKdtreeNode *kdtreenode;
  kdtreenode = gt_malloc(sizeof (GtKdtreeNode));
  kdtreenode->id_self = UNDEF;
  kdtreenode->id_high = UNDEF;
  kdtreenode->id_low = UNDEF;
  kdtreenode->key = key;
  kdtreenode->value = value;
}

void gt_kdtreenode_delete(GtKdtreeNode *kdtreenode)
{
  if (kdtreenode != NULL) {
    gt_free(kdtreenode);
  }
}

GtKdtree *gt_kdtree_new(GtUword dimension,
                        int (* cmp)(const void *, const void *, GtUword))
{
  GtKdtree *kdtree;
  kdtree = gt_malloc(sizeof (GtKdtree));
  kdtree->dimension = dimension;
  GT_INITARRAY(kdtree->nodes, GtKdtreeNode);
  kdtree->cmp = (* cmp)(const void *, const void *, GtUword);
  return kdtree;
}

void gt_kdtree_delete(GtKdtree *kdtree)
{
  GT_FREEARRAY(kdtree->nodes, GtKdtreeNode);
  gt_free(kdtree);
}

/* KdtreeNode methods */

void *gt_kdtreenode_value(const GtKdtreeNode *kdtreenode)
{
  gt_assert(kdtreenode);
  return kdtreenode->value;
}

/* Kdtree methods */

GtUword gt_kdtree_size(const GtKdtree *kdtree)
{
  return kdtree->nodes->nextfreeGtKdtreeNode;
}

GtKdtreeNode *gt_kdtree_get(const GtKdtree *kdtree, GtUword node_id)
{
  gt_assert(kdtree);
  gt_assert(node_id < gt_kdtree_size(kdtree));
  return kdtree->nodes->spaceGtKdtreeNode + node_id;
}

bool gt_kdtree_empty(const GtKdtree *kdtree)
{
  gt_assert(kdtreenode);
  if (gt_kdtree_size(kdtree) == 0) {
    return true;
  } else {
    return false;
  }
}

GtUword gt_kdtree_dimension(GtKdtree *kdtree)
{
  gt_assert(kdtree);
  return kdtree->dimension;
}

/* recursive insertion of a new node (key, value), starting from currid */
void gt_kdtree_insert_rec(GtKdtree *kdtree, const void *key, const void *value,
                          GtUword *currid, GtUword currdim)
{
  if (*currid == UNDEF) { /* insertion site found */
    GtKdtreeNode *new = gt_kdtreenode_new(key, value);
    *currid = new->id_self = gt_kdtree_size(kdtree);
    GT_STOREINARRAY(kdtree->nodes, GtKdtreeNode, ARR_INCR, new);
  } else {
    GtKdtreeNode *currnode = gt_kdtree_get(kdtree, *currid);
    const int cmpval = cmp(key, currnode->key, currdim);
    const GtUword nextdim = (currdim + 1) % kdtree->dimension;
    if (cmpval < 0) {
      gt_kdtree_insert_rec(kdtree, key, value, &currnode->id_low, nextdim);
    } else if (cmpval > 0) {
      gt_kdtree_insert_rec(kdtree, key, value, &currnode->id_high, nextdim);
    } else { /* equal keys: overwrite value */
      gt_kdtreenode_value(currnode) = value;
    }
  }
}

void gt_kdtree_insert(GtKdtree *kdtree, const void *key, const void *value)
{
  gt_assert(kdtree && key && value);
  GtUword currid = (gt_kdtree_empty(kdtree) ? UNDEF : 0);
  gt_kdtree_insert_rec(kdtree, key, value, &currid, 0);
}

GtKdtreeNode *gt_kdtree_find(const GtKdtree *kdtree, const void *key) {
  return NULL;
}

