#ifndef __UTILS_H__
#define __UTILS_H__
#include <mpi.h>
#include <string>
#include <vector>
#include "crystal.h"

using namespace std;
inline bool CheckLattice(latticeStruct ls)
{
    return ((ls.a > 0) &&
            (ls.b > 0) &&
            (ls.c > 0) && 
            (ls.alpha + ls.beta + ls.gamma < 360) &&
            (ls.alpha + ls.beta > ls.gamma) &&
            (ls.beta + ls.gamma > ls.alpha) &&
            (ls.gamma + ls.alpha > ls.beta));
}

void SortCrystal(Crystal *, vector<unsigned int>);
void WriteCrystal(Crystal *, string);
void GathervCrystal(Crystal *, Crystal *, MPI_Datatype, MPI_Datatype, int, MPI_Comm);
void BcastCrystal(Crystal *, MPI_Datatype, MPI_Datatype, int, MPI_Comm);
double GetMassFromSymbol(string);
#endif
