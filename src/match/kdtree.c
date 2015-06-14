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
#include <math.h>
#include "match/kdtree.h"
#include "core/arraydef.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "extended/prority_queue.h"

#define GT_KDTREE_UNDEF ULONG_MAX
#define GT_KDTREE_ARRAYINCR 256
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

struct GtKdtreeQueueElem {
  GtUword node_id;
  double priority;
};

/* Constructors and Destructors for GtKdtree and GtKdtreeNode */

GtKdtreeNode *gt_kdtreenode_new(const void *key, const void *value)
{
  GtKdtreeNode *kdtreenode;
  kdtreenode = gt_malloc(sizeof (GtKdtreeNode));
  kdtreenode->id_self = GT_KDTREE_UNDEF;
  kdtreenode->id_high = GT_KDTREE_UNDEF;
  kdtreenode->id_low = GT_KDTREE_UNDEF;
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

double square_dist(void *key_a, void* key_b, GtUword length) {
  gt_assert(key_a && key_b);
  GtUword dim;
  double diff, sum;
  sum = 0.0;
  for (dim = 0; dim < length; dim++) {
    diff = (double)(key_a[dim]) - (double)(key_b[dim]);
    sum += diff * diff;
  }
  return sum;
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
  if (*currid == GT_KDTREE_UNDEF) { /* insertion site found */
    GtKdtreeNode *new = gt_kdtreenode_new(key, value);
    *currid = new->id_self = gt_kdtree_size(kdtree);
    GT_STOREINARRAY(kdtree->nodes, GtKdtreeNode, GT_KDTREE_ARRAYINCR +
                    1.2*kdtree->nodes->allocatedGtKdtreeNode, new);
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
  GtUword currid = (gt_kdtree_empty(kdtree) ? GT_KDTREE_UNDEF : 0);
  gt_kdtree_insert_rec(kdtree, key, value, &currid, 0);
}

bool *gt_kdtree_find(const GtKdtree *kdtree, const void *key, GtUword *loc)
{
  GtKdtreeNode *kdtreenode;
  int cmpval;
  GtUword currdim, nextid = 0;
  *loc = GT_KDTREE_UNDEF;

  gt_assert(kdtree && key && loc);
  if (gt_kdtree_empty(kdtree))
    return false;

  for (currdim = 0; nextid != GT_KDTREE_UNDEF; currdim = (currdim+1)%kdtree->dimension) {
    *loc = nextid;
    kdtreenode = gt_kdtree_get(kdtree, *loc);
    cmpval = cmp(key, kdtreenode->key, currdim);
    if (cmpval < 0) {
      nextid = kdtreenode->id_low;
    } else if (cmpval > 0) {
      nextid = kdtreenode->id_high;
    } else { /* exact key match */
      return true;
    }
  }
  return false;
}

/* recursive knn-search, starting from currid */
void gt_kdtree_knn_rec(const GtKdtree *kdtree, const void *key, GtUword currid, GtUword currdim, GtPriorityQueue *queue)
{
  if (currid == GT_KDTREE_UNDEF) { /* no subtree */
    return;
  } else {
    GtUword nextid, nextdim;
    int cmpval;
    GtKdtreeNode *currnode = gt_kdtree_get(kdtree, currid);
    GtKdtreeQueueElem elem = {currid, -square_dist(key, currnode->key)};

    if (gt_priority_queue_is_full(queue)) {
      if (gt_double_larger_double(elem.priority,
                                  gt_priority_queue_find_min(queue))) {
        gt_priority_queue_extract_min(queue);
        gt_priority_queue_add(queue, &elem);
      }
    } else {
      gt_priority_queue_add(queue, &elem);
    }

    cmpval = cmp(key, currnode->key, currdim);
    nextid = (cmpval <= 0) ? currnode->id_low : currnode->id_high;
    nextdim = (currdim + 1) % kdtree->dimension;
    gt_kdtree_knn_rec(kdtree, key, nextid, nextdim, queue);

    if (!gt_priority_queue_is_full(queue) || gt_double_larger_double
        (-gt_priority_queue_find_min(queue),
         square_dist(key+currdim, currnode->key+currdim, 1))) {
      /* search other subtree */
      nextid = (cmpval <= 0) ? currnode->id_high : currnode->id_low;
      gt_kdtree_knn_rec(kdtree, key, nextid, nextdim, queue);
    }
  }
}

void gt_kdtree_knn(const GtKdtree *kdtree, const void *key, GtUint kvalue,
                   void *value)
{
  GtPriorityQueue *queue;
  GtKdtreeQueueElem *elem;
  GtUword currid;
  gt_assert(kdtree && key && value);
  if (kvalue < gt_kdtree_size(kdtree) || kvalue == 0) {
    /* error: not enough entries for k */
    return NULL;
  }
  queue = gt_priority_queue_new(queue_cmp, kvalue);
  currid = 0;
  gt_kdtree_knn_rec(kdtree, key, &currid, 0, queue);
  while (!gt_priority_queue_is_empty(queue)) {
    elem = (GtKdtreeQueueElem *)gt_priority_queue_extract_min(queue);
    value = gt_kdtreenode_value(gt_kdtree_get(kdtree, elem.node_id));
    value++;
  }
  gt_priority_queue_delete(queue);
}
