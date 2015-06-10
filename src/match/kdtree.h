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

typedef struct GtKdtree GtKdtree;
typedef struct GtKdtreeNode GtKdtreeNode;

GtKdtree *gt_kdtree_new(void);
void gt_kdtree_delete(GtKdtree *kdtree);

//GtKdtreeNode *gt_kdtreenode_new(GtUword ndim);
//void gt_kdtreenode_delete(GtKdtreeNode *kdtreenode);

GtUword gt_kdtree_insert(GtKdtree *kdtree, GtKdtreeNode *new, GtUword currdim);

#endif
