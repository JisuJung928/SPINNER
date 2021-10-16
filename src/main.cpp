#include <iostream>
#include <mpi.h>

#define LAMMPS_LIB_MPI
#include "library.h"

#include "config.h"
#include "input.h"

using namespace std;
int main(int argc, char* argv[])
{
    Input input;
    read_input(input, "./INPUT");

    Crystal c = generate(input);
    string config_name = "./POSCAR";
    string title = "./TEST_CONFIG";
    c.writePOSCAR(config_name, title);

    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* split comm
    int lammps;
    if (me < nprocs_lammps) lammps = 1;
    else lammps = MPI_UNDEFINED;
    MPI_Comm comm_lammps;
    MPI_Comm_split(MPI_COMM_WORLD,lammps,0,&comm_lammps);
    */

    void *lmp = NULL;
    lmp = lammps_open(0, NULL, MPI_COMM_WORLD, NULL);

    lammps_close(lmp);
    MPI_Finalize();
  
    return 0;
}
