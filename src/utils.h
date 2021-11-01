#ifndef __UTILS_H__
#define __UTILS_H__
#include <mpi.h>
#include "crystal.h"

void DistributeCrystal(Crystal *, MPI_Datatype, MPI_Datatype, int, MPI_Comm *);
#endif
