#ifndef __CONFIG_H__
#define __CONFIG_H__
#include "input.h"

void RandomGeneration(Input *, Crystal *, int, int);
void Crossover(Input *, Crystal *, int, Crystal *, int, int);
void Permutation(Input *, Crystal *, int, Crystal *, int, int);
void LatticeMutation(Input *, Crystal *, int, Crystal *, int, int);
#endif
