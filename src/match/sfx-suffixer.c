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
#include <stdio.h>
#include <limits.h>
#include <errno.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include "core/arraydef.h"
#include "core/assert_api.h"
#include "core/error_api.h"
#include "core/unused_api.h"
#include "core/progressbar.h"
#include "core/minmax.h"
#include "core/fa.h"
#include "spacedef.h"
#include "core/progress_timer.h"
#include "intcode-def.h"
#include "core/encodedsequence.h"
#include "core/safecast-gen.h"
#include "esa-fileend.h"
#include "sfx-partssuf-def.h"
#include "sfx-suffixer.h"
#include "sfx-bentsedg.h"
#include "sfx-enumcodes.h"
#include "sfx-strategy.h"
#include "diff-cover.h"
#include "sfx-copysort.h"
#include "stamp.h"

#include "sfx-mappedstr.pr"

static inline void setsortspace(Suftab *suftab,unsigned long idx,
                                unsigned long value)
{
  suftab->sortspace[idx - suftab->offset] = value;
}

struct Sfxiterator
{
  bool storespecials;
  Codetype currentmincode,
           currentmaxcode;
  unsigned long specialcharacters,
         widthofpart,
         totallength;
  Suftab suftab;
  unsigned long nextfreeCodeatposition;
  Codeatposition *spaceCodeatposition;
  Suftabparts *suftabparts;
  const GtEncodedsequence *encseq;
  GtReadmode readmode;
  Outlcpinfo *outlcpinfo;
  unsigned int part,
               numofchars,
               prefixlength;
  GtArrayGtUlong fusp;
  Specialrangeiterator *sri;
  GtSequencerange overhang;
  bool exhausted;
  Bcktab *bcktab;
  Codetype numofallcodes;
  unsigned long *leftborder; /* points to bcktab->leftborder */
  unsigned long long bucketiterstep; /* for progressbar */
  Sfxstrategy sfxstrategy;
  GtLogger *logger;
  GtProgressTimer *sfxprogress;
  Differencecover *dcov;
};

#ifdef SKDEBUG
static unsigned long iterproduceCodeatposition(Codeatposition *codelist,
                                               const  GtEncodedsequence *encseq,
                                               GtReadmode readmode,
                                               unsigned int prefixlength,
                                               unsigned int numofchars)
{
  if (prefixlength > 1U)
  {
    Enumcodeatposition *ecp;
    unsigned long insertindex;
    Specialcontext specialcontext;

    ecp = newEnumcodeatposition(encseq,
                                readmode,
                                prefixlength,
                                numofchars);
    for (insertindex = 0; nextEnumcodeatposition(&specialcontext,ecp);
         insertindex++)
    {
      codelist[insertindex].maxprefixindex = specialcontext.maxprefixindex;
      codelist[insertindex].position = specialcontext.position;
      codelist[insertindex].code
        = computefilledqgramcode(ecp,
                                 specialcontext.maxprefixindex,
                                 specialcontext.position -
                                 specialcontext.maxprefixindex);
    }
    freeEnumcodeatposition(&ecp);
    return insertindex;
  }
  return 0;
}

static void compareCodeatpositionlists(const Codeatposition *codelist1,
                                       unsigned long len1,
                                       const Codeatposition *codelist2,
                                       unsigned long len2)
{
  unsigned long idx;

  if (len1 != len2)
  {
    fprintf(stderr,"len1 = %lu != %lu = len2\n",len1,len2);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  for (idx=0; idx<len1; idx++)
  {
    if (codelist1[idx].position != codelist2[idx].position)
    {
      fprintf(stderr,"idx %lu, codelist1.position = " FormatSeqpos " != "
                      FormatSeqpos " = codelist2.position\n",idx,
                      PRINTSeqposcast(codelist1[idx].position),
                      PRINTSeqposcast(codelist2[idx].position));
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    if (codelist1[idx].maxprefixindex != codelist2[idx].maxprefixindex)
    {
      fprintf(stderr,"idx %lu, codelist1.maxprefixindex = %u != %u = "
                     "codelist2.maxprefixindex\n",idx,
                      codelist1[idx].maxprefixindex,
                      codelist2[idx].maxprefixindex);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    if (codelist1[idx].code != codelist2[idx].code)
    {
      fprintf(stderr,"idx %lu, codelist1.code = %u != %u = "
                     "codelist2.code\n",idx,
                      codelist1[idx].code,
                      codelist2[idx].code);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
}

static void verifycodelistcomputation(
                       const GtEncodedsequence *encseq,
                       GtReadmode readmode,
                       unsigned long realspecialranges,
                       unsigned int prefixlength,
                       unsigned int numofchars,
                       unsigned long nextfreeCodeatposition1,
                       const Codeatposition *spaceCodeatposition1)
{
  unsigned long nextfreeCodeatposition2;
  Codeatposition *spaceCodeatposition2;

  ALLOCASSIGNSPACE(spaceCodeatposition2,NULL,Codeatposition,
                   realspecialranges+1);
  nextfreeCodeatposition2 = iterproduceCodeatposition(spaceCodeatposition2,
                                                      encseq,
                                                      readmode,
                                                      prefixlength,
                                                      numofchars);
  gt_assert(realspecialranges+1 >= (unsigned long) nextfreeCodeatposition2);
  compareCodeatpositionlists(spaceCodeatposition1,
                             nextfreeCodeatposition1,
                             spaceCodeatposition2,
                             nextfreeCodeatposition2);
  FREESPACE(spaceCodeatposition2);
}
#endif

#ifdef SKDEBUG
static Codetype getencseqcode(const GtEncodedsequence *encseq,
                              GtReadmode readmode,
                              unsigned long totallength,
                              const Codetype **multimappower,
                              unsigned int prefixlength,
                              unsigned long pos)
{
  Codetype code = 0;
  unsigned int idx;
  GtUchar cc;

  for (idx=0; idx<prefixlength; idx++)
  {
    gt_assert((unsigned long) (pos + idx) < totallength);
    cc = gt_encodedsequence_getencodedcharnospecial(encseq,pos + idx, readmode);
    gt_assert(ISNOTSPECIAL(cc));
    code += multimappower[idx][cc];
  }
  return code;
}

static Codetype previouscode = 0;
static bool previousfirstspecialdefined = false,
            previousstorespecials = false;
unsigned int previousspecialpos = 0;
#endif

static void updatekmercount(void *processinfo,
                            Codetype code,
                            unsigned long position,
                            const Firstspecialpos *firstspecial)
{
  Sfxiterator *sfi = (Sfxiterator *) processinfo;

  if (firstspecial->defined)
  {
    if (sfi->storespecials)
    {
      if (firstspecial->specialpos > 0)
      {
        if (sfi->sfxstrategy.storespecialcodes)
        {
          Codeatposition *cp;

          cp = sfi->spaceCodeatposition + sfi->nextfreeCodeatposition++;
          gt_assert(code <= (Codetype) MAXCODEVALUE);
          cp->code = (unsigned int) code;
          gt_assert(firstspecial->specialpos <= MAXPREFIXLENGTH);
          cp->maxprefixindex = firstspecial->specialpos;
          cp->position = position + firstspecial->specialpos;
        }
        sfi->storespecials = false;
        gt_assert(code > 0);
        sfi->leftborder[code]++;
      }
    } else
    {
      if (firstspecial->specialpos > 0)
      {
        gt_assert(code > 0);
        sfi->leftborder[code]++;
      } else
      {
        sfi->storespecials = true;
      }
    }
  } else
  {
#ifdef SKDEBUG
    if (code == 0)
    {
      Codetype code2 = getencseqcode(sfi->encseq,
                                   sfi->readmode,
                                   gt_encodedsequence_total_length(sfi->encseq),
                                   bcktab_multimappower(sfi->bcktab),
                                   sfi->prefixlength,
                                   position);
      if (code2 != 0)
      {
        printf("### position " FormatSeqpos ", code2 = %lu != 0\n",
                        PRINTSeqposcast(position),code2);
        printf("previouscode = " FormatCodetype "\n",
                        previouscode);
        if (previousfirstspecialdefined)
        {
          printf("previousfirstspecialdefined = true\n");
          printf("previousstorespecials = %s\n",
                  previousstorespecials ? "true" : "false");
          printf("previousspecialpos = %u\n",previousspecialpos);
        }
        exit(GT_EXIT_PROGRAMMING_ERROR);
      }
    }
#endif
    sfi->leftborder[code]++;
  }
#ifdef SKDEBUG
  previouscode = code;
  previousfirstspecialdefined = firstspecial->defined;
  previousstorespecials = sfi->storespecials;
  previousspecialpos = firstspecial->specialpos;
#endif
}

static void insertwithoutspecial(void *processinfo,
                                 Codetype code,
                                 unsigned long position,
                                 const Firstspecialpos *firstspecial)
{
  if (!firstspecial->defined)
  {
    Sfxiterator *sfi = (Sfxiterator *) processinfo;

    if (code >= sfi->currentmincode && code <= sfi->currentmaxcode)
    {
      setsortspace(&sfi->suftab,--sfi->leftborder[code],position);
      /* from right to left */
    }
  }
}

static void reversespecialcodes(Codeatposition *spaceCodeatposition,
                                unsigned long nextfreeCodeatposition)
{
  Codeatposition *front, *back, tmp;

  for (front = spaceCodeatposition,
       back = spaceCodeatposition + nextfreeCodeatposition - 1;
       front < back; front++, back--)
  {
    tmp = *front;
    *front = *back;
    *back = tmp;
  }
}

static void derivespecialcodesfromtable(Sfxiterator *sfi,bool deletevalues)
{
  Codetype code;
  unsigned int prefixindex;
  unsigned long insertindex, j;
  unsigned long stidx;

  for (prefixindex=1U; prefixindex < sfi->prefixlength; prefixindex++)
  {
    for (j=0, insertindex = 0; j < sfi->nextfreeCodeatposition; j++)
    {
      if (prefixindex <= sfi->spaceCodeatposition[j].maxprefixindex)
      {
        code = codedownscale(sfi->bcktab,
                             (Codetype) sfi->spaceCodeatposition[j].code,
                             prefixindex,
                             sfi->spaceCodeatposition[j].maxprefixindex);
        if (code >= sfi->currentmincode && code <= sfi->currentmaxcode)
        {
          updatebckspecials(sfi->bcktab,code,sfi->numofchars,prefixindex);
          stidx = --sfi->leftborder[code];
          /* from right to left */
          setsortspace(&sfi->suftab,stidx,
                       sfi->spaceCodeatposition[j].position - prefixindex);
        }
      }
      if (deletevalues)
      {
        if (prefixindex < sfi->prefixlength - 1 &&
            prefixindex < sfi->spaceCodeatposition[j].maxprefixindex)
        {
          if (insertindex < j)
          {
            sfi->spaceCodeatposition[insertindex] =
              sfi->spaceCodeatposition[j];
          }
          insertindex++;
        }
      }
    }
    if (deletevalues)
    {
      sfi->nextfreeCodeatposition = insertindex;
    }
  }
}

static void derivespecialcodesonthefly(Sfxiterator *sfi)
{
  Codetype code;
  unsigned int prefixindex;
  unsigned long stidx;
  Enumcodeatposition *ecp;
  Specialcontext specialcontext;

  for (prefixindex=1U; prefixindex < sfi->prefixlength; prefixindex++)
  {
    ecp = newEnumcodeatposition(sfi->encseq,sfi->readmode,
                                sfi->prefixlength,
                                sfi->numofchars);
    while (nextEnumcodeatposition(&specialcontext,ecp))
    {
      if (prefixindex <= specialcontext.maxprefixindex)
      {
        if (computefilledqgramcodestopatmax(&code,
                                            ecp,
                                            prefixindex,
                                            specialcontext.position-prefixindex,
                                            sfi->currentmaxcode))
        {
          gt_assert(code <= sfi->currentmaxcode);
          if (code >= sfi->currentmincode)
          {
            updatebckspecials(sfi->bcktab,code,sfi->numofchars,prefixindex);
            gt_assert(code > 0);
            stidx = --sfi->leftborder[code];
            /* from right to left */
            setsortspace(&sfi->suftab,stidx,
                         specialcontext.position - prefixindex);
          }
        }
      }
    }
    freeEnumcodeatposition(&ecp);
  }
}

void freeSfxiterator(Sfxiterator **sfiptr)
{
  Sfxiterator *sfi = (Sfxiterator *) *sfiptr;
#ifdef SKDEBUG
  if (sfi->bcktab != NULL)
  {
    checkcountspecialcodes(sfi->bcktab);
  }
#endif
  if (sfi->bcktab != NULL)
  {
    addfinalbckspecials(sfi->bcktab,sfi->numofchars,sfi->specialcharacters);
  }
  if (sfi->sri != NULL)
  {
    freespecialrangeiterator(&sfi->sri);
  }
  FREESPACE(sfi->spaceCodeatposition);
  FREESPACE(sfi->suftab.sortspace);
  freesuftabparts(sfi->suftabparts);
  if (sfi->bcktab != NULL)
  {
    bcktab_delete(&sfi->bcktab);
  }
  if (sfi->dcov != NULL)
  {
    differencecover_delete(sfi->dcov);
  }
  FREESPACE(*sfiptr);
}

#ifdef SKDEBUG
static void showleftborder(const unsigned long *leftborder,
                           Codetype numofallcodes)
{
  Codetype i;

  for (i=0; i<MIN(numofallcodes,(Codetype) 1024); i++)
  {
    printf("leftborder[" FormatCodetype "]=" FormatSeqpos "\n",
            i,leftborder[i]);
  }
}
#endif

Sfxiterator *newSfxiterator(const GtEncodedsequence *encseq,
                            GtReadmode readmode,
                            unsigned int prefixlength,
                            unsigned int numofparts,
                            Outlcpinfo *outlcpinfo,
                            const Sfxstrategy *sfxstrategy,
                            GtProgressTimer *sfxprogress,
                            GtLogger *logger,
                            GtError *err)
{
  Sfxiterator *sfi = NULL;
  unsigned long realspecialranges, specialcharacters;
  bool haserr = false;

  gt_error_check(err);

  realspecialranges = getencseqrealspecialranges(encseq);
  specialcharacters = getencseqspecialcharacters(encseq);
  gt_assert(prefixlength > 0);
  if (sfxstrategy != NULL && sfxstrategy->storespecialcodes &&
      prefixlength > MAXPREFIXLENGTH)
  {
    gt_error_set(err,"argument for option -pl must be in the range [1,%u]",
                  MAXPREFIXLENGTH);
    haserr = true;
  }
  if (!haserr)
  {
    ALLOCASSIGNSPACE(sfi,NULL,Sfxiterator,1);
    if (sfxstrategy != NULL && sfxstrategy->storespecialcodes)
    {
      ALLOCASSIGNSPACE(sfi->spaceCodeatposition,NULL,
                       Codeatposition,realspecialranges+1);
      gt_logger_log(logger,"sizeof (spaceCodeatposition)=%lu",
                              (unsigned long) (sizeof (Codeatposition) *
                                               (realspecialranges+1)));
    } else
    {
      sfi->spaceCodeatposition = NULL;
    }
    sfi->bcktab = NULL;
    sfi->nextfreeCodeatposition = 0;
    sfi->suftab.sortspace = NULL;
    sfi->suftabparts = NULL;
    sfi->encseq = encseq;
    sfi->readmode = readmode;
    sfi->numofchars = gt_encodedsequence_alphabetnumofchars(encseq);
    sfi->prefixlength = prefixlength;
    sfi->dcov = NULL;
    if (sfxstrategy != NULL)
    {
       sfi->sfxstrategy = *sfxstrategy;
       if (sfxstrategy->cmpcharbychar || !possibletocmpbitwise(encseq))
       {
         sfi->sfxstrategy.cmpcharbychar = true;
       } else
       {
         sfi->sfxstrategy.cmpcharbychar = false;
       }
    } else
    {
      defaultsfxstrategy(&sfi->sfxstrategy,
                         possibletocmpbitwise(encseq) ? false : true);
    }
    gt_logger_log(logger,"maxinsertionsort=%lu",
                sfi->sfxstrategy.maxinsertionsort);
    gt_logger_log(logger,"maxbltriesort=%lu",
                sfi->sfxstrategy.maxbltriesort);
    gt_logger_log(logger,"maxcountingsort=%lu",
                sfi->sfxstrategy.maxcountingsort);
    gt_logger_log(logger,"storespecialcodes=%s",
                sfi->sfxstrategy.storespecialcodes ? "true" : "false");
    gt_logger_log(logger,"cmpcharbychar=%s",
                sfi->sfxstrategy.cmpcharbychar ? "true" : "false");
    if (sfi->sfxstrategy.ssortmaxdepth.defined)
    {
      gt_logger_log(logger,"ssortmaxdepth=%u",
                               sfi->sfxstrategy.ssortmaxdepth.valueunsignedint);
    } else
    {
      gt_logger_log(logger,"ssortmaxdepth=undefined");
    }
    sfi->totallength = gt_encodedsequence_total_length(encseq);
    gt_logger_log(logger,"totallength=%lu",
                        sfi->totallength);
    sfi->specialcharacters = specialcharacters;
    sfi->outlcpinfo = outlcpinfo;
    sfi->sri = NULL;
    sfi->part = 0;
    sfi->exhausted = false;
    sfi->bucketiterstep = 0;
    sfi->logger = logger;
    sfi->sfxprogress = sfxprogress;

    if (sfi->sfxstrategy.differencecover > 0 &&
        getencseqspecialcharacters(encseq)
          < gt_encodedsequence_total_length(encseq))
    {
      if (sfxprogress != NULL)
      {
        gt_progress_timer_start_new_state(sfxprogress,
                                          "sorting difference cover sample",
                                          stdout);
      }
      sfi->dcov = differencecover_new(sfi->sfxstrategy.differencecover,
                                      encseq,readmode,logger);
      if (sfi->dcov == NULL)
      {
        gt_error_set(err,"no difference cover modulo %u found",
                     sfi->sfxstrategy.differencecover);
        haserr = true;
      } else
      {
        if (differencecover_vparamverify(sfi->dcov,err) != 0)
        {
          haserr = true;
          differencecover_delete(sfi->dcov);
          sfi->dcov = NULL;
        } else
        {
          gt_logger_log(logger,"presorting sample suffixes according to "
                                  "difference cover modulo %u",
                                  sfi->sfxstrategy.differencecover);
          differencecover_sortsample(sfi->dcov,sfi->sfxstrategy.cmpcharbychar,
                                     false);
        }
      }
    }
  }
  if (!haserr)
  {
    gt_assert(sfi != NULL);
    sfi->bcktab = allocBcktab(sfi->numofchars,
                              prefixlength,
                              sfi->sfxstrategy.storespecialcodes,
                              logger,
                              err);
    if (sfi->bcktab == NULL)
    {
      haserr = true;
      sfi->leftborder = NULL;
      sfi->numofallcodes = 0;
    } else
    {
      sfi->leftborder = bcktab_leftborder(sfi->bcktab);
      sfi->numofallcodes = bcktab_numofallcodes(sfi->bcktab);
    }
  }
  if (!haserr)
  {
    gt_assert(sfi != NULL);
    sfi->storespecials = true;
    if (sfxprogress != NULL)
    {
      gt_progress_timer_start_new_state(sfxprogress,
                                        "counting prefix distribution",
                                        stdout);
    }
    getencseqkmers(encseq,
                   readmode,
                   updatekmercount,
                   sfi,
                   prefixlength);
    if (sfi->sfxstrategy.storespecialcodes)
    {
      gt_assert(realspecialranges+1
                  >= (unsigned long) sfi->nextfreeCodeatposition);
      reversespecialcodes(sfi->spaceCodeatposition,sfi->nextfreeCodeatposition);
    }
#ifdef SKDEBUG
    verifycodelistcomputation(encseq,
                              readmode,
                              realspecialranges,
                              prefixlength,
                              sfi->numofchars,
                              sfi->nextfreeCodeatposition,
                              sfi->spaceCodeatposition);
#endif
    gt_assert(sfi->leftborder != NULL);
#ifdef SKDEBUG
    showleftborder(sfi->leftborder,sfi->numofallcodes);
#endif
    bcktab_leftborderpartialsums(sfi->bcktab,
                                 sfi->totallength - specialcharacters);
    sfi->suftabparts = newsuftabparts(numofparts,
                                      sfi->leftborder,
                                      sfi->numofallcodes,
                                      sfi->totallength - specialcharacters,
                                      specialcharacters + 1,
                                      logger);
    gt_assert(sfi->suftabparts != NULL);
    ALLOCASSIGNSPACE(sfi->suftab.sortspace,NULL,GtUlong,
                     stpgetlargestwidth(sfi->suftabparts));
    sfi->suftab.longest.defined = false;
    sfi->suftab.longest.valueunsignedlong = 0;
    if (hasspecialranges(sfi->encseq))
    {
      sfi->sri = newspecialrangeiterator(sfi->encseq,
                                         GT_ISDIRREVERSE(sfi->readmode)
                                           ? false : true);
    } else
    {
      sfi->sri = NULL;
    }
    sfi->fusp.spaceGtUlong = sfi->suftab.sortspace;
    sfi->fusp.allocatedGtUlong = stpgetlargestwidth(sfi->suftabparts);
    sfi->overhang.leftpos = sfi->overhang.rightpos = 0;
  }
  if (haserr)
  {
    if (sfi != NULL)
    {
      freeSfxiterator(&sfi);
    }
    return NULL;
  }
  return sfi;
}

bool sfi2longestsuffixpos(unsigned long *longest,const Sfxiterator *sfi)
{
  if (sfi->suftab.longest.defined)
  {
    *longest = sfi->suftab.longest.valueunsignedlong;
    return true;
  }
  return false;
}

static void preparethispart(Sfxiterator *sfi)
{
  unsigned long partwidth;
  unsigned int numofparts = stpgetnumofparts(sfi->suftabparts);

  if (sfi->part == 0 &&
      sfi->sfxprogress != NULL && gt_progress_timer_use_bar(sfi->sfxprogress))
  {
    gt_progressbar_start(&sfi->bucketiterstep,
                         (unsigned long long) sfi->numofallcodes);
  }
  sfi->currentmincode = stpgetcurrentmincode(sfi->part,sfi->suftabparts);
  sfi->currentmaxcode = stpgetcurrentmaxcode(sfi->part,sfi->suftabparts);
  sfi->widthofpart = stpgetcurrentwidthofpart(sfi->part,sfi->suftabparts);
  /*
  sfi->suftabptr = sfi->suftab -
                   stpgetcurrentsuftaboffset(sfi->part,sfi->suftabparts);
  */
  sfi->suftab.offset = stpgetcurrentsuftaboffset(sfi->part,sfi->suftabparts);
  if (sfi->sfxstrategy.storespecialcodes)
  {
    derivespecialcodesfromtable(sfi,(numofparts == 1U) ? true : false);
  } else
  {
    derivespecialcodesonthefly(sfi);
  }
  if (sfi->sfxprogress != NULL)
  {
    gt_progress_timer_start_new_state(sfi->sfxprogress,
                                      "inserting suffixes into buckets",
                                      stdout);
  }
  getencseqkmers(sfi->encseq,
                 sfi->readmode,
                 insertwithoutspecial,
                 sfi,
                 sfi->prefixlength);
  if (sfi->sfxprogress != NULL)
  {
    gt_progress_timer_start_new_state(sfi->sfxprogress,
                                      "sorting the buckets",
                                      stdout);
  }
  partwidth = stpgetcurrentsumofwdith(sfi->part,sfi->suftabparts);
  if (sfi->sfxstrategy.ssortmaxdepth.defined &&
      sfi->prefixlength == sfi->sfxstrategy.ssortmaxdepth.valueunsignedint)
  {
    if (!sfi->sfxstrategy.streamsuftab)
    {
      qsufsort(sfi->suftab.sortspace,
               -1,
               &sfi->suftab.longest.valueunsignedlong,
               sfi->encseq,
               sfi->readmode,
               sfi->currentmincode,
               sfi->currentmaxcode,
               partwidth,
               sfi->bcktab,
               sfi->numofchars,
               sfi->prefixlength,
               false,
               true,
               sfi->outlcpinfo);
      sfi->suftab.longest.defined = true;
    }
  } else
  {
    GtBucketspec2 *bucketspec2 = NULL;
    gt_assert(!sfi->sfxstrategy.streamsuftab);
    if (numofparts == 1U && sfi->outlcpinfo == NULL && sfi->prefixlength >= 2U)
    {
      bucketspec2 = gt_bucketspec2_new(sfi->bcktab,sfi->encseq,sfi->readmode,
                                       partwidth,sfi->numofchars);
    }
    if (sfi->sfxstrategy.differencecover > 0)
    {
      sortbucketofsuffixes(sfi->suftab.sortspace - sfi->suftab.offset,
                           bucketspec2,
                           (unsigned long) partwidth,
                           sfi->encseq,
                           sfi->readmode,
                           sfi->currentmincode,
                           sfi->currentmaxcode,
                           sfi->bcktab,
                           sfi->numofchars,
                           sfi->prefixlength,
                           &sfi->sfxstrategy,
                           (void *) sfi->dcov,
                           dc_sortunsortedbucket,
                           sfi->logger);
    } else
    {
      sortallbuckets (&sfi->suftab,
                      bucketspec2,
                      sfi->encseq,
                      sfi->readmode,
                      sfi->currentmincode,
                      sfi->currentmaxcode,
                      partwidth,
                      sfi->bcktab,
                      sfi->numofchars,
                      sfi->prefixlength,
                      sfi->outlcpinfo,
                      &sfi->sfxstrategy,
                      &sfi->bucketiterstep,
                      sfi->logger);
    }
    if (bucketspec2 != NULL)
    {
      unsigned long *suftabptr = sfi->suftab.sortspace - sfi->suftab.offset;
      gt_copysortsuffixes(bucketspec2,suftabptr,sfi->logger);
      gt_bucketspec2_delete(bucketspec2);
      bucketspec2 = NULL;
    }
  }
  sfi->part++;
}

int postsortsuffixesfromstream(Sfxiterator *sfi, const GtStr *str_indexname,
                               GtError *err)
{
  int mmapfiledesc = -1;
  GtStr *tmpfilename;
  struct stat sb;
  bool haserr = false;

  if (sfi->totallength == sfi->specialcharacters)
  {
    return 0;
  }
  FREESPACE(sfi->suftab.sortspace);
  if (sfi->sfxstrategy.streamsuftab)
  {
    gt_assert(sfi->sfxstrategy.ssortmaxdepth.defined &&
              sfi->prefixlength ==
              sfi->sfxstrategy.ssortmaxdepth.valueunsignedint);
  }
  tmpfilename = gt_str_clone(str_indexname);
  gt_str_append_cstr(tmpfilename,SUFTABSUFFIX);
  mmapfiledesc = open(gt_str_get(tmpfilename), O_RDWR, 0);
  if (mmapfiledesc == -1)
  {
    gt_error_set(err,"cannot open file \"%s\": %s",gt_str_get(tmpfilename),
                 strerror(errno));
    haserr = true;
  }
  if (!haserr && fstat(mmapfiledesc, &sb) == -1)
  {
    gt_error_set(err,"cannot fstat file \"%s\": %s",gt_str_get(tmpfilename),
                 strerror(errno));
    haserr = true;
  }
  if (!haserr && sizeof (off_t) > sizeof (size_t) && sb.st_size > SIZE_MAX)
  {
    gt_error_set(err,"file \"%s\" of size %llu is too large to map",
                 gt_str_get(tmpfilename),(unsigned long long) sb.st_size);
    haserr = true;
  }
  if (!haserr
        && (size_t) sb.st_size != sizeof (unsigned long) * (sfi->totallength+1))
  {
    gt_error_set(err,"mapping file %s: file size "
                     " = %lu != %lu = expected number of units",
                 gt_str_get(tmpfilename),
                 (unsigned long) sb.st_size,
                 (unsigned long) (sfi->totallength+1) * sizeof (unsigned long));
    haserr = true;
  }
  if (!haserr)
  {
    gt_assert(sfi->totallength >= sfi->specialcharacters);
    qsufsort(NULL,
             mmapfiledesc,
             &sfi->suftab.longest.valueunsignedlong,
             sfi->encseq,
             sfi->readmode,
             0,
             sfi->currentmaxcode,
             sfi->totallength - sfi->specialcharacters,
             sfi->bcktab,
             sfi->numofchars,
             sfi->prefixlength,
             sfi->sfxstrategy.hashexceptions,
             sfi->sfxstrategy.absoluteinversesuftab,
             sfi->outlcpinfo);
    sfi->suftab.longest.defined = true;
  }
  gt_str_delete(tmpfilename);
  if (close(mmapfiledesc) == -1)
  {
    gt_error_set(err,"cannot close file \"%s\": %s",gt_str_get(tmpfilename),
                 strerror(errno));
    haserr = true;
  }
  return haserr ? -1 : 0;
}

static void insertfullspecialrange(Sfxiterator *sfi,
                                   unsigned long leftpos,
                                   unsigned long rightpos)
{
  unsigned long pos;

  gt_assert(leftpos < rightpos);
  if (GT_ISDIRREVERSE(sfi->readmode))
  {
    pos = rightpos - 1;
  } else
  {
    pos = leftpos;
  }
  while (true)
  {
    if (GT_ISDIRREVERSE(sfi->readmode))
    {
      sfi->fusp.spaceGtUlong[sfi->fusp.nextfreeGtUlong++]
        = GT_REVERSEPOS(sfi->totallength,pos);
      if (pos == leftpos)
      {
        break;
      }
      pos--;
    } else
    {
      sfi->fusp.spaceGtUlong[sfi->fusp.nextfreeGtUlong++] = pos;
      if (pos == rightpos-1)
      {
        break;
      }
      pos++;
    }
  }
}

static void fillspecialnextpage(Sfxiterator *sfi)
{
  GtSequencerange range;
  unsigned long width;

  while (true)
  {
    if (sfi->overhang.leftpos < sfi->overhang.rightpos)
    {
      width = sfi->overhang.rightpos - sfi->overhang.leftpos;
      if (sfi->fusp.nextfreeGtUlong + width > sfi->fusp.allocatedGtUlong)
      {
        /* does not fit into the buffer, so only output a part */
        unsigned long rest = sfi->fusp.nextfreeGtUlong +
                             width - sfi->fusp.allocatedGtUlong;
        gt_assert(rest > 0);
        if (GT_ISDIRREVERSE(sfi->readmode))
        {
          insertfullspecialrange(sfi,sfi->overhang.leftpos + rest,
                                 sfi->overhang.rightpos);
          sfi->overhang.rightpos = sfi->overhang.leftpos + rest;
        } else
        {
          insertfullspecialrange(sfi,sfi->overhang.leftpos,
                                     sfi->overhang.rightpos - rest);
          sfi->overhang.leftpos = sfi->overhang.rightpos - rest;
        }
        break;
      }
      if (sfi->fusp.nextfreeGtUlong + width == sfi->fusp.allocatedGtUlong)
      { /* overhang fits into the buffer and buffer is full */
        insertfullspecialrange(sfi,sfi->overhang.leftpos,
                               sfi->overhang.rightpos);
        sfi->overhang.leftpos = sfi->overhang.rightpos = 0;
        break;
      }
      /* overhang fits into the buffer and buffer is not full */
      insertfullspecialrange(sfi,sfi->overhang.leftpos,
                             sfi->overhang.rightpos);
      sfi->overhang.leftpos = sfi->overhang.rightpos = 0;
    } else
    {
      if (sfi->sri != NULL && nextspecialrangeiterator(&range,sfi->sri))
      {
        width = range.rightpos - range.leftpos;
        gt_assert(width > 0);
        if (sfi->fusp.nextfreeGtUlong + width > sfi->fusp.allocatedGtUlong)
        { /* does not fit into the buffer, so only output a part */
          unsigned long rest = sfi->fusp.nextfreeGtUlong +
                               width - sfi->fusp.allocatedGtUlong;
          if (GT_ISDIRREVERSE(sfi->readmode))
          {
            insertfullspecialrange(sfi,range.leftpos + rest,
                                   range.rightpos);
            sfi->overhang.leftpos = range.leftpos;
            sfi->overhang.rightpos = range.leftpos + rest;
          } else
          {
            insertfullspecialrange(sfi,range.leftpos,range.rightpos - rest);
            sfi->overhang.leftpos = range.rightpos - rest;
            sfi->overhang.rightpos = range.rightpos;
          }
          break;
        }
        if (sfi->fusp.nextfreeGtUlong + width == sfi->fusp.allocatedGtUlong)
        { /* overhang fits into the buffer and buffer is full */
          insertfullspecialrange(sfi,range.leftpos,range.rightpos);
          sfi->overhang.leftpos = sfi->overhang.rightpos = 0;
          break;
        }
        insertfullspecialrange(sfi,range.leftpos,range.rightpos);
        sfi->overhang.leftpos = sfi->overhang.rightpos = 0;
      } else
      {
        if (sfi->fusp.nextfreeGtUlong < sfi->fusp.allocatedGtUlong)
        {
          sfi->fusp.spaceGtUlong[sfi->fusp.nextfreeGtUlong++]
                                                             = sfi->totallength;
          sfi->exhausted = true;
        }
        break;
      }
    }
  }
}

const unsigned long *nextSfxiterator(unsigned long *numberofsuffixes,
                                     bool *specialsuffixes,
                                     Sfxiterator *sfi)
{
  if (sfi->part < stpgetnumofparts(sfi->suftabparts))
  {
    preparethispart(sfi);
    *numberofsuffixes = sfi->widthofpart;
    *specialsuffixes = false;
    return sfi->suftab.sortspace;
  }
  if (sfi->exhausted)
  {
    if (sfi->sfxprogress != NULL && gt_progress_timer_use_bar(sfi->sfxprogress))
    {
      gt_progressbar_stop();
    }
    return NULL;
  }
  sfi->fusp.nextfreeGtUlong = 0;
  fillspecialnextpage(sfi);
  gt_assert(sfi->fusp.nextfreeGtUlong > 0);
  *numberofsuffixes = (unsigned long) sfi->fusp.nextfreeGtUlong;
  *specialsuffixes = true;
  return sfi->suftab.sortspace;
}

int sfibcktab2file(FILE *fp,
                   const Sfxiterator *sfi,
                   GtError *err)
{
  gt_error_check(err);
  return bcktab2file(fp,sfi->bcktab,err);
}

unsigned int getprefixlenbits(void)
{
  return (unsigned int) PREFIXLENBITS;
}
