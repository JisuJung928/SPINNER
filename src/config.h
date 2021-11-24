#ifndef __CONFIG_H__
#define __CONFIG_H__
#include <mpi.h>
#include "input.h"

void RandomGeneration(Input *, Crystal *, int);
void Crossover(Input *, Crystal *, int, Crystal *, int, MPI_Comm);
void Permutation(Input *, Crystal *, int, Crystal *, int);
void LatticeMutation(Input *, Crystal *, int, Crystal *, int);
#endif
