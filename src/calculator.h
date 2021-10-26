#ifndef __CALCULATOR_H__
#define __CALCULATOR_H__
#include <mpi.h>

#include "randSpg.h"
#include "input.h"

void *LammpsInit(Input *, Crystal *, MPI_Comm *, char *);
double Relax(Input *, Crystal *, MPI_Comm *, int);
#endif
