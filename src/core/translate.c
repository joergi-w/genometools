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

#include <assert.h>
#include "core/codon.h"
#include "core/ma.h"
#include "core/translate.h"

void gt_translate_dna(GT_Str *protein, const char *dnaseq, unsigned long dnalen,
                      unsigned int frame)
{
  const char *dnaptr;
  assert(protein && !gt_str_length(protein) && dnaseq && frame < 3);
  /* translate the DNA in forward direction */
  for (dnaptr = dnaseq + frame; dnaptr < dnaseq + dnalen - 2; dnaptr += 3) {
    gt_str_append_char(protein,
                       gt_codon2amino(*dnaptr, *(dnaptr+1), *(dnaptr+2)));
  }
}

static int encode(char nucleotide)
{
  switch (nucleotide) {
    case 'A':
    case 'a':
      return GT_A_CODE;
    case 'C':
    case 'c':
      return GT_C_CODE;
    case 'G':
    case 'g':
      return GT_G_CODE;
  }
  return GT_T_CODE;
}

void gt_translate_all_frames(char **frame1, char **frame2, char **frame3,
                             const char *dna_sequence, unsigned long seqlen)
{
  unsigned long i, frame1len, frame2len, frame3len;
  int codon;
  assert(frame1 && frame2 && frame3 && dna_sequence);

  if (seqlen < 3)
    return; /* nothing to translate here */

  /* allocate appropriate space for frames (incl. terminal '\0') */
  frame1len = seqlen / 3;
  frame2len = (seqlen - 1) / 3;
  frame3len = (seqlen - 2) / 3;
  (*frame1) = gt_malloc(frame1len + 1);
  (*frame2) = gt_malloc(frame2len + 1);
  (*frame3) = gt_malloc(frame3len + 1);

  /* encode first two nucleotides */
  codon = (encode(dna_sequence[0]) << 2) | encode(dna_sequence[1]);

  /* iterate over DNA sequence */
  for (i = 0; i < seqlen - 2; i++) {
    /* compute next codon code */
    codon = ((codon << 2) | encode(dna_sequence[i+2])) & 0x3f; /* 0..0111111 */

    /* store amino acid in appropriate frame */
    switch (i % GT_CODON_LENGTH) {
      case 0:
        (*frame1)[i/GT_CODON_LENGTH] = gt_aminos[codon];
        break;
      case 1:
        (*frame2)[i/GT_CODON_LENGTH] = gt_aminos[codon];
        break;
      case 2:
        (*frame3)[i/GT_CODON_LENGTH] = gt_aminos[codon];
    }
  }

  /* set terminal '\0's */
  (*frame1)[frame1len] = '\0';
  (*frame2)[frame2len] = '\0';
  (*frame3)[frame3len] = '\0';
}
