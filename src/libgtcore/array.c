/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg

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
#include <limits.h>
#include <string.h>
#include <stdlib.h>
#include <libgtcore/array.h>
#include <libgtcore/dynalloc.h>
#include <libgtcore/ensure.h>
#include <libgtcore/range.h>
#include <libgtcore/xansi.h>

#define NUM_OF_TESTS    100
#define MAX_SIZE        1024

struct Array {
  void *space;
  unsigned long next_free;
  size_t allocated,
         size_of_elem;
};

Array* array_new(size_t size_of_elem, Env *env)
{
  Array *a = env_ma_calloc(env, 1, sizeof (Array));
  assert(size_of_elem);
  a->size_of_elem = size_of_elem;
  return a;
}

void* array_get(const Array *a, unsigned long idx)
{
  assert(a && idx < a->next_free);
  return a->space + idx * a->size_of_elem;
}

void* array_get_first(const Array *a)
{
  return array_get(a, 0);
}

void* array_get_last(const Array *a)
{
  assert(a->next_free);
  return array_get(a, a->next_free-1);
}

void* array_pop(Array *a)
{
  assert(a && a->next_free);
  a->next_free--;
  return a->space + a->next_free * a->size_of_elem;
}

void array_rem(Array *a, unsigned long idx)
{
  unsigned long i;
  assert(a && idx < a->next_free);
  /* move elements */
  for (i = idx+1; i < a->next_free; i++) {
    memcpy(a->space + (i-1) * a->size_of_elem, a->space + i * a->size_of_elem,
           a->size_of_elem);
  }
  /* remove last (now duplicated) element */
  a->next_free--;
}

void array_reverse(Array *a, Env *env)
{
  void *front, *back, *tmp;
  assert(a);
  tmp = env_ma_malloc(env, a->size_of_elem);
  for (front = a->space, back = a->space + (a->next_free-1) * a->size_of_elem;
       front < back;
       front += a->size_of_elem, back -= a->size_of_elem) {
    memcpy(tmp, front, a->size_of_elem);
    memcpy(front, back, a->size_of_elem);
    memcpy(back, tmp, a->size_of_elem);
  }
  env_ma_free(tmp, env);
}

void* array_get_space(const Array *a)
{
  assert(a);
  return a->space;
}

void array_add_elem(Array *a, void *elem, size_t size_of_elem, Env *env)
{
  assert(a && elem);
  assert(a->size_of_elem == size_of_elem);
  assert(a->next_free <= a->allocated);
  /* make sure we have enough space */
  if ((a->next_free + 1) * a->size_of_elem > a->allocated)
    a->space = dynalloc(a->space, &a->allocated,
                        (a->next_free + 1) * a->size_of_elem, env);
  /* add */
  memcpy(a->space + a->next_free * a->size_of_elem, elem, a->size_of_elem);
  a->next_free++;
}

void array_add_array(Array *dest, const Array *src, Env *env)
{
  unsigned long i;
  assert(dest && src && dest->size_of_elem == src->size_of_elem);
  for (i = 0; i < array_size(src); i++)
    array_add_elem(dest, array_get(src, i), src->size_of_elem, env);
}

size_t array_elem_size(const Array *a)
{
  assert(a);
  return a->size_of_elem;
}

unsigned long array_size(const Array *a)
{
  return a ? a->next_free : 0;
}

void array_set_size(Array *a, unsigned long size)
{
  assert(a);
  assert(size <= a->next_free);
  a->next_free = size;
}

void array_reset(Array *a)
{
  assert(a);
  a->next_free = 0;
}

Array* array_clone(const Array *a, Env *env)
{
  Array *a_copy;
  assert(a);
  a_copy = env_ma_malloc(env, sizeof (Array));
  /* XXX: overflow checks -> safemult(next_free, size_of_elem) */
  a_copy->space = env_ma_malloc(env, a->next_free * a->size_of_elem);
  memcpy(a_copy->space, a->space, a->next_free * a->size_of_elem);
  a_copy->next_free = a_copy->allocated = a->next_free;
  a_copy->size_of_elem = a->size_of_elem;
  return a_copy;
}

void array_sort(Array *a,int(*compar)(const void *, const void *))
{
  qsort(a->space,a->next_free,a->size_of_elem,compar);
}

int array_compare(Array *a,Array *b,
                  int(*compar)(const void *, const void *,Env *),
                  Env *env)
{
  unsigned long idx;
  size_t size_a, size_b;
  int cmp;

  size_a = array_size(a);
  size_b = array_size(b);
  if (size_a < size_b)
  {
    env_error_set(env,"array_size(a) = %lu < %lu = array_size(b)",
                  (unsigned long) size_a,
                  (unsigned long) size_b);
    return -1;
  }
  if (size_a > size_b)
  {
    env_error_set(env,"array_size(a) = %lu > %lu = array_size(b)",
                  (unsigned long) size_a,
                  (unsigned long) size_b);
    return 1;
  }
  for (idx=0; idx<(unsigned long) size_a; idx++)
  {
    cmp = compar(array_get(a,idx),array_get(b,idx),env);
    if(cmp != 0)
    {
      return cmp;
    }
  }
  return 0;
}

int array_example(Env *env)
{
  unsigned long i;
  Array *a;

  env_error_check(env);

  /* an example array use case */

  a = array_new(sizeof (unsigned long), env);
  for (i = 0; i < 100; i++) {
    array_add(a, i, env);
    assert(i == *(unsigned long*) array_get(a, i));
  }
  assert(array_size(a) == 100);
  assert(*(unsigned long*) array_pop(a) == 99);
  assert(array_size(a) == 99);
  array_delete(a, env);

  return 0;
}

int array_unit_test(Env *env)
{
  Array *char_array, *int_array, *a = NULL;
  char cc, *char_array_test;
  int ci, *int_array_test;
  unsigned long i, j, size;
  Range range;
  int had_err = 0;
  env_error_check(env);

  /* testing an empty array */
  char_array = array_new(sizeof (char), env);
  array_delete(char_array, env);
  int_array = array_new(sizeof (int), env);
  array_delete(int_array, env);

  char_array = array_new(sizeof (char), env);
  int_array = array_new(sizeof (int), env);
  char_array_test = env_ma_malloc(env, (MAX_SIZE + 1) * sizeof (char));
  int_array_test = env_ma_malloc(env, MAX_SIZE * sizeof (int));

  for (i = 0; i < NUM_OF_TESTS && !had_err; i++) {
    size = ((double) rand() / RAND_MAX) * MAX_SIZE;

    array_reset(char_array);
    array_reset(int_array);

    ensure(had_err, array_size(char_array) == 0);
    ensure(had_err, array_size(int_array) == 0);

    for (i = 0; i < size && !had_err; i++) {
      cc = ((double) rand() / RAND_MAX) * CHAR_MAX;
      ci = ((double) rand() / RAND_MAX) * INT_MAX;

      array_add(char_array, cc, env);
      array_add(int_array, ci, env);

      ensure(had_err, array_size(char_array) == i+1);
      ensure(had_err, array_size(int_array) == i+1);
      ensure(had_err, *((char*) array_get(char_array, i)) == cc);
      ensure(had_err, *((int*) array_get(int_array, i)) == ci);

      array_add_elem(char_array, &cc, sizeof (char), env);
      array_add_elem(int_array, &ci, sizeof (int), env);

      ensure(had_err, array_size(char_array) == i+2);
      ensure(had_err, array_size(int_array) == i+2);
      ensure(had_err, *((char*) array_get(char_array, i+1)) == cc);
      ensure(had_err, *((int*) array_get(int_array, i+1)) == ci);
      ensure(had_err, *((char*) array_pop(char_array)) == cc);
      ensure(had_err, *((int*) array_pop(int_array)) == ci);
      ensure(had_err, array_size(char_array) == i+1);
      ensure(had_err, array_size(int_array) == i+1);
      ensure(had_err, *((char*) array_get(char_array, i)) == cc);
      ensure(had_err, *((int*) array_get(int_array, i)) == ci);

      char_array_test[i] = cc;
      char_array_test[i+1]= '\0';
      int_array_test[i] = ci;

      ensure(had_err, strncmp(array_get_space(char_array), char_array_test,
                              strlen(char_array_test)) == 0);

      for (j = 0; j <= i && !had_err; j++)
        ensure(had_err, *(int*) array_get(int_array, j) == int_array_test[j]);
    }
  }

  /* test array_reverse() */
  if (!had_err) {
    a = array_new(sizeof (Range), env);
    for (i = 0; i < 24; i++) {
      range.start = i + 1;
      range.end   = i + 101;
      array_add(a, range, env);
    }
    array_reverse(a, env);
    for (i = 0; !had_err && i < 24; i++) {
      range.start = i + 1;
      range.end   = i + 101;
      ensure(had_err, !range_compare(range, *(Range*) array_get(a, 23 - i)));
    }
  }
  array_delete(a, env);

  array_delete(char_array, env);
  array_delete(int_array, env);
  env_ma_free(char_array_test, env);
  env_ma_free(int_array_test, env);

  return had_err;
}

void array_delete(Array *a, Env *env)
{
  if (!a) return;
  env_ma_free(a->space, env);
  env_ma_free(a, env);
}
