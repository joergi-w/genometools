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

#ifndef KDTREE_H
#define KDTREE_H

#include <stdbool.h>

typedef struct GtKdtree GtKdtree;
typedef struct GtKdtreeNode GtKdtreeNode;

GtKdtree *gt_kdtree_new(GtUword dimension,
                        int (* cmp)(const void *, const void *, GtUword));
void gt_kdtree_delete(GtKdtree *kdtree);

/* return a pointer to the value of the specified node */
void *gt_kdtreenode_value(const GtKdtreeNode *kdtreenode);

/* return the number of nodes in kdtree */
GtUword gt_kdtree_size(const GtKdtree *kdtree);

/* return a pointer to the node with specified node_id */
GtKdtreeNode *gt_kdtree_get(const GtKdtree *kdtree, GtUword node_id);

/* return whether the specified kdtree is empty */
bool gt_kdtree_empty(const GtKdtree *kdtree);

/* return the dimension of the keys in kdtree */
GtUword gt_kdtree_dimension(GtKdtree *kdtree);

/* Create and insert a node of key and value into kdtree. */
void gt_kdtree_insert(GtKdtree *kdtree, const void* key, const void* value);
#endif
