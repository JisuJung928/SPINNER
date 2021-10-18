#include "calculator.h"

void *LammpsInit(Input input, Crystal crystal)
{
    /* split comm
    int lammps;
    if (me < nprocs_lammps) lammps = 1;
    else lammps = MPI_UNDEFINED;
    MPI_Comm comm_lammps;
    MPI_Comm_split(MPI_COMM_WORLD,lammps,0,&comm_lammps);
    */

    void *lmp = nullptr;
    lmp = lammps_open(0, nullptr, MPI_COMM_WORLD, nullptr);
    // TODO: set cell

    return lmp;
}
