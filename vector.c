
#include "linear_algebra.h"

#include <stdlib.h>
#include <stdio.h>

// ----------------------------------------------------------------------
// vector_create
//
// create a new vector of length "n"

struct vector *
vector_create(int n)
{
  struct vector *v = malloc(sizeof(*v));
  v->n = n;
  v->vals = calloc(n, sizeof(*v->vals));

  return v;
}

// ----------------------------------------------------------------------
// vector_create_and_set
//
// create a new vector of length "n", and set the values from "vals"

struct vector *
vector_create_and_set(int n, const double *vals)
{
  struct vector *v = vector_create(n);

  for (int i = 0; i < v->n; i++) {
    VEC(v, i) = vals[i];
  }

  return v;
}

// ----------------------------------------------------------------------
// vector_destroy
//
// destroys (deallocates) the vector "v"

void
vector_destroy(struct vector *v)
{
  free(v->vals);
  free(v);
}

// ----------------------------------------------------------------------
// vector_is_equal

bool
vector_is_equal(const struct vector *x, const struct vector *y)
{
  for (int i = 0; i < x->n; i++) {
    if (VEC(x, i) != VEC(y, i)) {
      return false;
    }
  }
  return true;
}

// vector_print
//
// prints vector "v" to stdout

void
vector_print(struct vector *v)
{
  printf("[(#=%d)", v->n);
  for (int i = 0; i < v->n; i++) {
    printf(" %g", VEC(v, i));
  }
  printf("]");
}

// ----------------------------------------------------------------------
// vector_create_crd_cc
//
// create a uniformly spaced cell-centered coordinate vector
// with n elements covering [0:len]

struct vector *
vector_create_crd_cc(int n, double len)
{
  struct vector *v = vector_create(n);

  double dx = len / n;
  for (int i = 0; i < n; i++) {
    VEC(v, i) = (i + .5) * dx;
  }

  return v;
}

// ----------------------------------------------------------------------
// vector_create_crd_nc
//
// create a uniformly spaced node-centered coordinate vector
// with n elements covering [0:len]

struct vector *
vector_create_crd_nc(int n, double len)
{
  struct vector *v = vector_create(n);

  double dx = len / (n-1);
  for (int i = 0; i < n; i++) {
    VEC(v, i) = i * dx;
  }
  
  return v;
}

// ----------------------------------------------------------------------
// vector_write
// 
// write a vector on the given coordinates into an ASCII file

void
vector_write(struct vector *crd, struct vector *v, const char *filename)
{
  assert(crd->n == v->n);

  FILE *file = fopen(filename, "w");
  
  for (int i = 0; i < v->n; i++) {
    fprintf(file, "%g %g\n", VEC(crd, i), VEC(v, i));
  }

  fclose(file);
}
