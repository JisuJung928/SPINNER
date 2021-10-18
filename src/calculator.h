#ifndef __CALCULATOR_H__
#define __CALCULATOR_H__
#include <mpi.h>
#define LAMMPS_LIB_MPI
#include "library.h"
#include "randSpg.h"

#include "input.h"

void *LammpsInit(Input, Crystal);
void *Relax(Input, Crystal &);
#endif
