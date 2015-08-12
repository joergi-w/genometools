/*
  Copyright (C) 2015 Annika Seidel, annika.seidel@studium.uni-hamburg.de
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

#include <string.h>
#include <inttypes.h>
#include "core/error.h"
#include "core/ma_api.h"
#include "core/versionfunc.h"
#include "core/option_api.h"
#include "core/str.h"
#include "core/str_array.h"
#include "core/types_api.h"
#include "extended/linearalign_affinegapcost.h"
#include "extended/linearalign.h"
#include "match/test-pairwise.h"
#include "tools/gt_paircmp.h"

typedef struct
{
  GtStr *charlist;
  GtUword len;
} Charlistlen;

typedef struct
{
  GtStrArray *strings,
             *files;
  Charlistlen *charlistlen;
  GtStr *text;
  bool showedist;
  bool print;
} Cmppairwiseopt;

static void showsimpleoptions(const Cmppairwiseopt *opt)
{
  if (gt_str_array_size(opt->strings) > 0)
  {
    if (!opt->showedist)
      printf("# two strings \"%s\" \"%s\"\n", gt_str_array_get(opt->strings,0),
             gt_str_array_get(opt->strings,1UL));
    return;
  }
  if (gt_str_array_size(opt->files) > 0)
  {
    printf("# two files \"%s\" \"%s\"\n", gt_str_array_get(opt->files,0),
           gt_str_array_get(opt->files,1UL));
    return;
  }
  if (opt->charlistlen != NULL)
  {
    printf("# alphalen \"%s\" " GT_WU "\n",
           gt_str_get(opt->charlistlen->charlist),
           opt->charlistlen->len);
    return;
  }
  if (gt_str_length(opt->text) > 0)
  {
    printf("# text \"%s\"\n", gt_str_get(opt->text));
    return;
  }
}

static GtOPrval parse_options(int *parsed_args,
                              Cmppairwiseopt *pw,
                              int argc, const char **argv, GtError *err)
{
  GtOptionParser *op;
  GtOption *optionstrings,
         *optionfiles,
         *optioncharlistlen,
         *optiontext,
         *optionshowedist,
         *optionprint;
  GtStrArray *charlistlen;
  GtOPrval oprval;

  gt_error_check(err);
  charlistlen = gt_str_array_new();
  pw->strings = gt_str_array_new();
  pw->files = gt_str_array_new();
  pw->text = gt_str_new();
  pw->charlistlen = NULL;
  pw->showedist = false;
  pw->print = false;
  op = gt_option_parser_new("options", "Apply function to pairs of strings.");
  gt_option_parser_set_mail_address(op, "<kurtz@zbh.uni-hamburg.de>");

  optionstrings = gt_option_new_string_array("ss", "use two strings",
                                             pw->strings);
  gt_option_parser_add_option(op, optionstrings);

  optionfiles = gt_option_new_filename_array("ff", "use two files",
                                             pw->files);
  gt_option_parser_add_option(op, optionfiles);

  optioncharlistlen = gt_option_new_string_array("a",
                                             "use character list and length",
                                             charlistlen);
  gt_option_parser_add_option(op, optioncharlistlen);

  optiontext = gt_option_new_string("t", "use text", pw->text, NULL);
  gt_option_parser_add_option(op, optiontext);

  optionshowedist = gt_option_new_bool("e", "output unit edit distance",
                      &pw->showedist, false);
  gt_option_parser_add_option(op, optionshowedist);

  optionprint = gt_option_new_bool("p", "print edist alignment",
                      &pw->print, false);
  gt_option_parser_add_option(op, optionprint);

  gt_option_exclude(optionstrings, optionfiles);
  gt_option_exclude(optionstrings, optioncharlistlen);
  gt_option_exclude(optionstrings, optiontext);
  gt_option_exclude(optionfiles, optioncharlistlen);
  gt_option_exclude(optionfiles, optiontext);
  gt_option_exclude(optioncharlistlen, optiontext);
  gt_option_imply(optionshowedist, optionstrings);
  gt_option_imply(optionprint, optionstrings);

  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, gt_versionfunc,
                                  err);
  if (oprval == GT_OPTION_PARSER_OK)
  {
    if (gt_option_is_set(optionstrings))
    {
      if (gt_str_array_size(pw->strings) != 2UL)
      {
        gt_error_set(err, "option -ss requires two string arguments");
        oprval = GT_OPTION_PARSER_ERROR;
      }
    } else
    {
      if (gt_option_is_set(optionfiles))
      {
        if (gt_str_array_size(pw->files) != 2UL)
        {
          gt_error_set(err, "option -ff requires two filename arguments");
          oprval = GT_OPTION_PARSER_ERROR;
        }
      } else
      {
        if (gt_option_is_set(optioncharlistlen))
        {
          GtWord readint;
          if (gt_str_array_size(charlistlen) != 2UL)
          {
            gt_error_set(err,
                         "option -a requires charlist and length argument");
            oprval = GT_OPTION_PARSER_ERROR;
          }else
          {
            pw->charlistlen = gt_malloc(sizeof *pw->charlistlen);
            pw->charlistlen->charlist =
              gt_str_ref(gt_str_array_get_str(charlistlen,
                                                                    0));
            if (sscanf(gt_str_array_get(charlistlen,1UL), GT_WD, &readint) != 1
                || readint < 1L)
            {
              gt_error_set(err,
                           "option -a requires charlist and length argument");
              oprval = GT_OPTION_PARSER_ERROR;
            }
            pw->charlistlen->len = (GtUword) readint;
          }
        } else
        {
          if (!gt_option_is_set(optiontext))
          {
            gt_error_set(err,
                         "use exactly one of the options -ss, -ff, -a, -t");
            oprval = GT_OPTION_PARSER_ERROR;
          }
        }
      }
    }
  }
  gt_option_parser_delete(op);
  if (oprval == GT_OPTION_PARSER_OK && *parsed_args != argc)
  {
    gt_error_set(err, "superfluous program parameters");
    oprval = GT_OPTION_PARSER_ERROR;
  }
  gt_str_array_delete(charlistlen);
  return oprval;
}

static void freesimpleoption(Cmppairwiseopt *cmppairwise)
{
  gt_str_array_delete(cmppairwise->strings);
  gt_str_array_delete(cmppairwise->files);
  gt_str_delete(cmppairwise->text);
  if (cmppairwise->charlistlen != NULL)
  {
    gt_str_delete(cmppairwise->charlistlen->charlist);
    gt_free(cmppairwise->charlistlen);
  }
}

static GtUword applycheckfunctiontosimpleoptions(
                                  Checkcmppairfuntype checkfunction,
                                  const Cmppairwiseopt *opt)
{
  if (gt_str_array_size(opt->strings) > 0)
  {
    bool forward = true;
    while (true)
    {
      checkfunction(forward,
                    (const GtUchar *) gt_str_array_get(opt->strings,0),
                    (GtUword) strlen(gt_str_array_get(opt->strings,0)),
                    (const GtUchar *) gt_str_array_get(opt->strings,1UL),
                    (GtUword) strlen(gt_str_array_get(opt->strings,1UL)));
      if (!forward)
      {
        break;
      }
      forward = false;
    }
    return 2UL; /* number of testcases */
  }
  if (gt_str_array_size(opt->files) > 0)
  {
    gt_runcheckfunctionontwofiles(checkfunction,
                               gt_str_array_get(opt->files,0),
                               gt_str_array_get(opt->files,1UL));
    return 2UL;
  }
  if (opt->charlistlen != NULL)
  {
    return gt_runcheckfunctiononalphalen(checkfunction,
                                      gt_str_get(opt->charlistlen->charlist),
                                      opt->charlistlen->len);
  }
  if (gt_str_length(opt->text) > 0)
  {
    return gt_runcheckfunctionontext(checkfunction, gt_str_get(opt->text));
  }
  gt_assert(false);
  return 0;
}

int gt_paircmp(int argc, const char **argv, GtError *err)
{
  int parsed_args;
  Cmppairwiseopt cmppairwise;
  GtOPrval oprval;

  gt_error_check(err);

  oprval = parse_options(&parsed_args, &cmppairwise, argc, argv, err);
  if (oprval == GT_OPTION_PARSER_OK)
  {
    gt_assert(parsed_args == argc);
    showsimpleoptions(&cmppairwise);
    if (cmppairwise.showedist)
    {
      GtUword edist, len1, len2;
      GtStr *s1, *s2;

      gt_assert(gt_str_array_size(cmppairwise.strings) >= 2);
      s1 = gt_str_array_get_str(cmppairwise.strings,0);
      s2 = gt_str_array_get_str(cmppairwise.strings,1UL);
      len1 = gt_str_length(s1);
      len2 = gt_str_length(s2);
      edist = gt_computegreedyunitedist((const GtUchar *) gt_str_get(s1),
                                        len1,
                                        (const GtUchar *) gt_str_get(s2),
                                        len2);
      printf(GT_WU " " GT_WU " " GT_WU " " GT_WU "%% errors\n", 
             edist, len1,len2,(200 * edist)/(len1+len2));
    }
    else if (cmppairwise.print)
    {
      gt_print_edist_alignment(
        (const GtUchar *) gt_str_array_get(cmppairwise.strings,0),0,
        (GtUword) strlen(gt_str_array_get(cmppairwise.strings,0)),
        (const GtUchar *) gt_str_array_get(cmppairwise.strings,1UL),0,
        (GtUword) strlen(gt_str_array_get(cmppairwise.strings,1UL)));
    } else
    {
      GtUword testcases;
      testcases = applycheckfunctiontosimpleoptions(gt_checkgreedyunitedist,
                                                    &cmppairwise);
      printf("# number of testcases for gt_checkgreedyunitedist: " GT_WU "\n",
              testcases);
      testcases = applycheckfunctiontosimpleoptions(gt_checklinearspace,
                                                    &cmppairwise);
      printf("# number of testcases for gt_checklinearspace: " GT_WU "\n",
              testcases);
      testcases = applycheckfunctiontosimpleoptions(gt_checklinearspace_local,
                                                    &cmppairwise);
      printf("# number of testcases for gt_checklinearspace_local: " GT_WU "\n",
              testcases);
      testcases = applycheckfunctiontosimpleoptions(gt_checkaffinelinearspace,
                                                    &cmppairwise);
      printf("# number of testcases for gt_checkaffinelinearspace: " GT_WU "\n",
              testcases);
      testcases = applycheckfunctiontosimpleoptions(gt_checkaffinelinearspace_local,
                                                    &cmppairwise);
      printf("# number of testcases for gt_checkaffinelinearspace: " GT_WU "\n",
              testcases);
    }
  }
  freesimpleoption(&cmppairwise);
  if (oprval == GT_OPTION_PARSER_REQUESTS_EXIT)
  {
    return 0;
  }
  if (oprval == GT_OPTION_PARSER_ERROR)
  {
    return -1;
  }
  return 0;
}
