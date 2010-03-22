/*
  Copyright (c) 2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#include <stdio.h>
#include <string.h>
#include "core/chardef.h"
#include "core/divmodmul.h"
#include "core/ma_api.h"

#include "core/encodedsequence.h"

typedef struct
{
  unsigned long lowerbound,
         upperbound,
         rank;
} Rankedbounds;

Rankedbounds *fillrankbounds(const GtEncodedsequence *encseq,
                             GtReadmode readmode)
{
  if (hasspecialranges(encseq))
  {
    Specialrangeiterator *sri;
    GtSequencerange range;
    unsigned long currentrank = 0, realspecialranges;
    Rankedbounds *rankedbounds, *rbptr;

    realspecialranges = getencseqrealspecialranges(encseq);
    rankedbounds = gt_malloc(sizeof (Rankedbounds) * realspecialranges);
    sri = newspecialrangeiterator(encseq,
                                  GT_ISDIRREVERSE(readmode)
                                  ? false : true);
    for (rbptr = rankedbounds; nextspecialrangeiterator(&range,sri); rbptr++)
    {
      rbptr->lowerbound = range.leftpos;
      rbptr->upperbound = range.rightpos;
      rbptr->rank = currentrank;
      currentrank += rbptr->upperbound - rbptr->lowerbound;
    }
    gt_assert(rbptr == rankedbounds + realspecialranges);
    freespecialrangeiterator(&sri);
    return rankedbounds;
  }
  return NULL;
}

unsigned long frompos2rank(const Rankedbounds *leftptr,
                    const Rankedbounds *rightptr,
                    unsigned long specialpos)
{
  const Rankedbounds *midptr;

  while (leftptr <= rightptr)
  {
    midptr = leftptr + GT_DIV2((unsigned long) (rightptr-leftptr));
    if (specialpos < midptr->lowerbound)
    {
      rightptr = midptr-1;
    } else
    {
      if (specialpos >= midptr->upperbound)
      {
        leftptr = midptr + 1;
      } else
      {
        return midptr->rank + specialpos - midptr->lowerbound;
      }
    }
  }
  fprintf(stderr,"frompos2rank: cannot find pos %lu"
                 " in ranges",specialpos);
  exit(GT_EXIT_PROGRAMMING_ERROR);
  /*@ignore@*/
  return 0;
  /*@end@*/
}

unsigned long fromrank2pos(const Rankedbounds *leftptr,
                    const Rankedbounds *rightptr,
                    unsigned long rank)
{
  const Rankedbounds *midptr;

  while (leftptr <= rightptr)
  {
    midptr = leftptr + GT_DIV2((unsigned long) (rightptr-leftptr));
    if (rank < midptr->rank)
    {
      rightptr = midptr-1;
    } else
    {
      if (rank >= midptr->rank + (midptr->upperbound - midptr->lowerbound))
      {
        leftptr = midptr + 1;
      } else
      {
        return midptr->lowerbound + (rank - midptr->rank);
      }
    }
  }
  fprintf(stderr,"fromrank2rank: cannot find rank %lu"
                 " in ranges",rank);
  exit(GT_EXIT_PROGRAMMING_ERROR);
  /*@ignore@*/
  return 0;
  /*@end@*/
}

typedef struct
{
  unsigned long specialrank,
         key;
} Specialrank;

static int compareSpecialrank(const void *a,const void *b)
{
  const Specialrank *aptr = (const Specialrank *) a,
                    *bptr = (const Specialrank *) b;

  if (aptr->key < bptr->key)
  {
    return -1;
  }
  if (aptr->key > bptr->key)
  {
    return 1;
  }
  gt_assert(false);
  return 0;
}

Specialrank *fillspecialranklist(const GtEncodedsequence *encseq,
                                 GtReadmode readmode,
                                 const unsigned long *inversesuftab)
{
  if (hasspecialranges(encseq))
  {
    Specialrangeiterator *sri;
    GtSequencerange range;
    unsigned long realspecialranges, specialrank, totallength;
    Specialrank *specialranklist, *rbptr;

    totallength = gt_encodedsequence_total_length(encseq);
    realspecialranges = getencseqrealspecialranges(encseq);
    specialranklist = gt_malloc(sizeof (Specialrank) * realspecialranges);
    sri = newspecialrangeiterator(encseq,
                                  GT_ISDIRREVERSE(readmode)
                                  ? false : true);
    rbptr = specialranklist;
    specialrank = 0;
    while (nextspecialrangeiterator(&range,sri))
    {
      gt_assert(rbptr < specialranklist + realspecialranges);
      gt_assert(range.rightpos<=totallength);
      specialrank += range.rightpos - range.leftpos;
      rbptr->specialrank = specialrank - 1;
      rbptr->key = inversesuftab[range.rightpos];
      rbptr++;
    }
    gt_assert(rbptr == specialranklist + realspecialranges);
    freespecialrangeiterator(&sri);
    qsort(specialranklist,(size_t) realspecialranges,
          sizeof (Specialrank),compareSpecialrank);
    return specialranklist;
  }
  return NULL;
}
