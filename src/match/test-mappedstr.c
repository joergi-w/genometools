/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "core/arraydef.h"
#include "core/chardef.h"
#include "core/error.h"
#include "core/unused_api.h"
#include "core/encodedsequence.h"
#include "intcode-def.h"
#include "sfx-nextchar.h"
#include "kmer2string.h"

#include "sfx-mappedstr.pr"

static Codetype qgram2codefillspecial(unsigned int numofchars,
                                      unsigned int kmersize,
                                      const GtEncodedsequence *encseq,
                                      GtReadmode readmode,
                                      unsigned long startpos,
                                      unsigned long totallength)
{
  Codetype integercode;
  unsigned long pos;
  bool foundspecial;
  GtUchar cc;

  if (startpos >= totallength)
  {
    integercode = (Codetype) (numofchars - 1);
    foundspecial = true;
  } else
  {
    /* for testing */
    cc = gt_encodedsequence_getencodedchar(encseq,startpos,readmode);
    if (ISSPECIAL(cc))
    {
      integercode = (Codetype) (numofchars - 1);
      foundspecial = true;
    } else
    {
      integercode = (Codetype) cc;
      foundspecial = false;
    }
  }
  for (pos = startpos + 1; pos < startpos + kmersize; pos++)
  {
    if (foundspecial)
    {
      ADDNEXTCHAR(integercode,numofchars-1,numofchars);
    } else
    {
      if (pos >= totallength)
      {
        ADDNEXTCHAR(integercode,numofchars-1,numofchars);
        foundspecial = true;
      } else
      {
        /* for testing */
        cc = gt_encodedsequence_getencodedchar(encseq,pos,readmode);
        if (ISSPECIAL(cc))
        {
          ADDNEXTCHAR(integercode,numofchars-1,numofchars);
          foundspecial = true;
        } else
        {
          ADDNEXTCHAR(integercode,cc,numofchars);
        }
      }
    }
  }
  return integercode;
}

GT_DECLAREARRAYSTRUCT(Codetype);

static void outkmeroccurrence(void *processinfo,
                              Codetype code,
                              GT_UNUSED unsigned long position,
                              GT_UNUSED const Firstspecialpos
                                           *firstspecialposition)
{
  GtArrayCodetype *codelist = (GtArrayCodetype *) processinfo;

  GT_STOREINARRAY(codelist,Codetype,1024,code);
}

/*
   The function to collect the code from a stream of fasta files
   can only produce the sequence of code in forward mode.
   Hence we compute the corresponding sequence also in GT_READMODE_FORWARD.
   Thus we restrict the call for verifymappedstr to the case where
   the suffix array is in readmode = GT_READMODE_FORWARD.
*/

static void collectkmercode(GtArrayCodetype *codelist,
                            const GtEncodedsequence *encseq,
                            unsigned int kmersize,
                            unsigned int numofchars,
                            unsigned long stringtotallength)
{
  unsigned long offset;
  Codetype code;

  for (offset=0; offset<=stringtotallength; offset++)
  {
    code = qgram2codefillspecial(numofchars,
                                 kmersize,
                                 encseq,
                                 GT_READMODE_FORWARD,
                                 offset,
                                 stringtotallength);
    GT_STOREINARRAY(codelist,Codetype,1024,code);
  }
}

static int comparecodelists(const GtArrayCodetype *codeliststream,
                            const GtArrayCodetype *codeliststring,
                            unsigned int kmersize,
                            unsigned int numofchars,
                            const char *characters,
                            GtError *err)
{
  unsigned long i;
  char buffer1[64+1], buffer2[64+1];

  gt_error_check(err);
  if (codeliststream->nextfreeCodetype != codeliststring->nextfreeCodetype)
  {
    gt_error_set(err,"length codeliststream= %lu != %lu =length codeliststring",
                  (unsigned long) codeliststream->nextfreeCodetype,
                  (unsigned long) codeliststring->nextfreeCodetype);
    return -1;
  }
  for (i=0; i<codeliststream->nextfreeCodetype; i++)
  {
    if (codeliststream->spaceCodetype[i] != codeliststring->spaceCodetype[i])
    {
      fromkmercode2string(buffer1,
                          codeliststream->spaceCodetype[i],
                          numofchars,
                          kmersize,
                          characters);
      fromkmercode2string(buffer2,
                          codeliststring->spaceCodetype[i],
                          numofchars,
                          kmersize,
                          characters);
      gt_error_set(err,"codeliststream[%lu] = " FormatCodetype " != "
                    FormatCodetype " = codeliststring[%lu]\n%s != %s",
                    i,
                    codeliststream->spaceCodetype[i],
                    codeliststring->spaceCodetype[i],
                    i,
                    buffer1,
                    buffer2);
      return -1;
    }
  }
  return 0;
}

static int verifycodelists(const GtEncodedsequence *encseq,
                           unsigned int kmersize,
                           unsigned int numofchars,
                           const GtArrayCodetype *codeliststream,
                           GtError *err)
{
  bool haserr = false;
  GtArrayCodetype codeliststring;
  const GtUchar *characters;
  unsigned long stringtotallength;

  gt_error_check(err);
  stringtotallength = gt_encodedsequence_total_length(encseq);
  characters = gt_encodedsequence_alphabetcharacters(encseq);
  GT_INITARRAY(&codeliststring,Codetype);
  collectkmercode(&codeliststring,
                  encseq,
                  kmersize,
                  numofchars,
                  stringtotallength);
  if (comparecodelists(codeliststream,
                       &codeliststring,
                       kmersize,
                       numofchars,
                       (const char *) characters,
                       err) != 0)
  {
    haserr = true;
  }
  GT_FREEARRAY(&codeliststring,Codetype);
  return haserr ? -1 : 0;
}

int verifymappedstr(const GtEncodedsequence *encseq,unsigned int prefixlength,
                    GtError *err)
{
  unsigned int numofchars;
  GtArrayCodetype codeliststream;
  bool haserr = false;

  gt_error_check(err);
  numofchars = gt_encodedsequence_alphabetnumofchars(encseq);
  GT_INITARRAY(&codeliststream,Codetype);
  if (getfastastreamkmers(gt_encodedsequence_filenames(encseq),
                          outkmeroccurrence,
                          &codeliststream,
                          numofchars,
                          prefixlength,
                          gt_encodedsequence_alphabetsymbolmap(encseq),
                          false,
                          err) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (verifycodelists(encseq,
                        prefixlength,
                        numofchars,
                        &codeliststream,
                        err) != 0)
    {
      haserr = true;
    }
  }
  GT_FREEARRAY(&codeliststream,Codetype);
  return haserr ? -1 : 0;
}
