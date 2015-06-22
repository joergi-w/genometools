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
#include "extended/kdtree.h"
#include "core/arraydef.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "extended/priority_queue.h"

#define GT_KDTREE_UNDEF ULONG_MAX
#define GT_KDTREE_ARRAYINCR 256
GT_DECLAREARRAYSTRUCT(GtKdtreeNode);

struct GtKdtree
{
  GtUword dimension;
  GtArrayGtKdtreeNode *nodes;
  GtKdtreeAccessFunction access;
  GtKdtreeCompare cmp;
};

struct GtKdtreeNode {
  GtUword id_low;
  GtUword id_high;
  const void *key;
  const void *value;
};

typedef struct GtKdtreeQueueElem {
  GtUword node_id;
  double priority;
} GtKdtreeQueueElem;

/* Constructors and Destructors */

GtKdtree *gt_kdtree_new(GtUword dimension, GtKdtreeCompare cmp,
                        GtKdtreeAccessFunction access)
{
  GtKdtree *kdtree = gt_malloc(sizeof *kdtree);
  kdtree->dimension = dimension;
  GT_INITARRAY(kdtree->nodes, GtKdtreeNode);
  kdtree->access = access;
  kdtree->cmp = cmp;
  return kdtree;
}

void gt_kdtree_delete(GtKdtree *kdtree)
{
  gt_assert(kdtree);
  GT_FREEARRAY(kdtree->nodes, GtKdtreeNode);
  gt_free(kdtree);
}

GtKdtreeNode *gt_kdtreenode_new(const void *key, const void *value)
{
  GtKdtreeNode *kdtreenode = gt_malloc(sizeof *kdtreenode);
  kdtreenode->id_high = GT_KDTREE_UNDEF;
  kdtreenode->id_low = GT_KDTREE_UNDEF;
  kdtreenode->key = key;
  kdtreenode->value = value;
  return kdtreenode;
}

void gt_kdtreenode_delete(GtKdtreeNode *kdtreenode)
{
  if (kdtreenode != NULL)
    gt_free(kdtreenode);
}

/* KdtreeNode methods */

const void *gt_kdtreenode_value(const GtKdtreeNode *kdtreenode)
{
  gt_assert(kdtreenode);
  return kdtreenode->value;
}

/* Kdtree methods */

static double gt_kdtree_square_dist(const void *key_a, const void* key_b,
                                    GtKdtreeAccessFunction access,
                                    GtUword dim_start, GtUword dim_length) {
  GtUword dim;
  double sum = 0.0;
  gt_assert(access != NULL && key_a && key_b);
  for (dim = dim_start; dim < dim_start + dim_length; dim++) {
    double diff = access(key_a, dim) - access(key_b, dim);
    sum += diff * diff;
  }
  return sum;
}

GtUword gt_kdtree_size(const GtKdtree *kdtree)
{
  gt_assert(kdtree);
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
  gt_assert(kdtree);
  return gt_kdtree_size(kdtree) == 0 ? true : false;
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
  gt_assert(kdtree && key && value && currid);
  if (*currid == GT_KDTREE_UNDEF) { /* insertion site found */
    GtKdtreeNode *new_node = gt_kdtreenode_new(key, value);
    *currid = gt_kdtree_size(kdtree);
    GT_STOREINARRAY(kdtree->nodes, GtKdtreeNode, GT_KDTREE_ARRAYINCR +
                    1.2*kdtree->nodes->allocatedGtKdtreeNode, *new_node);
  } else {
    GtKdtreeNode *currnode = gt_kdtree_get(kdtree, *currid);
    const int cmpval = kdtree->cmp(key, currnode->key, currdim);
    const GtUword nextdim = (currdim + 1) % kdtree->dimension;
    if (cmpval < 0) {
      gt_kdtree_insert_rec(kdtree, key, value, &currnode->id_low, nextdim);
    } else if (cmpval > 0) {
      gt_kdtree_insert_rec(kdtree, key, value, &currnode->id_high, nextdim);
    } else { /* equal keys: overwrite value */
      currnode->value = value;
    }
  }
}

void gt_kdtree_insert(GtKdtree *kdtree, const void *key, const void *value)
{
  GtUword currid;
  gt_assert(kdtree && key && value);

  currid = (gt_kdtree_empty(kdtree) ? GT_KDTREE_UNDEF : 0);
  gt_kdtree_insert_rec(kdtree, key, value, &currid, 0);
}

bool gt_kdtree_find(const GtKdtree *kdtree, const void *key, GtUword *loc)
{
  GtUword currdim, nextid = 0;
  gt_assert(kdtree && key && loc);

  *loc = GT_KDTREE_UNDEF;
  if (!gt_kdtree_empty(kdtree)) {
    for (currdim = 0; nextid != GT_KDTREE_UNDEF;
         currdim = (currdim+1)%kdtree->dimension) {
      GtKdtreeNode *kdtreenode = gt_kdtree_get(kdtree, nextid);
      const int cmpval = kdtree->cmp(key, kdtreenode->key, currdim);
      *loc = nextid;
      if (cmpval < 0) {
        nextid = kdtreenode->id_low;
      } else if (cmpval > 0) {
        nextid = kdtreenode->id_high;
      } else { /* exact key match */
        return true;
      }
    }
  }
  return false;
}

/* recursive knn-search, starting from currid */
void gt_kdtree_knn_rec(const GtKdtree *kdtree, const void *key, GtUword currid,
                       GtUword currdim, GtPriorityQueue *queue)
{
  gt_assert(kdtree && key && queue);
  if (currid != GT_KDTREE_UNDEF) { /* has subtree */
    GtKdtreeQueueElem elem;
    double maxdist, currdist;
    GtKdtreeNode *currnode = gt_kdtree_get(kdtree, currid);
    const int cmpval = kdtree->cmp(key, currnode->key, currdim);
    const GtUword nextdim = (currdim + 1) % kdtree->dimension;
    GtUword nextid = (cmpval <= 0) ? currnode->id_low : currnode->id_high;

    /* store negative distance, such that find_min gives the maximum */
    elem.node_id = currid;
    elem.priority = -1.0 * gt_kdtree_square_dist(key, currnode->key,
                                                 kdtree->access, 0,
                                                 gt_kdtree_size(kdtree));
    if (gt_priority_queue_is_full(queue)) {
      maxdist = -1.0 * *(double *)gt_priority_queue_find_min(queue);
      if (gt_double_smaller_double(-elem.priority, maxdist)) {
        gt_priority_queue_extract_min(queue);
        gt_priority_queue_add(queue, &elem);
      }
    } else {
      gt_priority_queue_add(queue, &elem);
    }

    gt_kdtree_knn_rec(kdtree, key, nextid, nextdim, queue);

    maxdist = -1.0 * *(double *)gt_priority_queue_find_min(queue);
    currdist = gt_kdtree_square_dist(key, currnode->key, kdtree->access,
                                     currdim, 1);
    if (!gt_priority_queue_is_full(queue) ||
        gt_double_smaller_double(currdist, maxdist)) {
      /* search other subtree */
      nextid = (cmpval <= 0) ? currnode->id_high : currnode->id_low;
      gt_kdtree_knn_rec(kdtree, key, nextid, nextdim, queue);
    }
  }
}

const void **gt_kdtree_knn(const GtKdtree *kdtree, const void *key,
                           GtUword kvalue)
{
  gt_assert(kdtree && key);
  const void **knn = gt_malloc(kvalue * sizeof *knn);
  if (kvalue >= gt_kdtree_size(kdtree) && kvalue != 0) {
    GtPriorityQueue *queue = gt_priority_queue_new((GtCompare)gt_double_compare,
                                                   kvalue);
    gt_kdtree_knn_rec(kdtree, key, 0, 0, queue);
    while (!gt_priority_queue_is_empty(queue)) {
      GtKdtreeQueueElem *elem;
      elem = (GtKdtreeQueueElem *)gt_priority_queue_extract_min(queue);
      *knn = gt_kdtreenode_value(gt_kdtree_get(kdtree, elem->node_id));
      knn++;
    }
    gt_priority_queue_delete(queue);
    return knn;
  }
  return NULL;
  /* else error: not enough entries for k */
}
