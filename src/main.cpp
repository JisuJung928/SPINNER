#include <iostream>
#include <mpi.h>

#define LAMMPS_LIB_MPI
#include "library.h"

#include "config.h"
#include "input.h"

using namespace std;
int main(int argc, char* argv[])
{
    Input input = read_input("./INPUT");

    vector<Crystal> c_vector = generate(input);
    /* test */
    string config_name = "./POSCAR";
    string title = "./TEST_CONFIG";
    for (unsigned int i = 0; i < c_vector.size(); ++i) {
        c_vector[i].writePOSCAR(config_name + to_string(i), title);
    }

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

    // TODO: resort the Crystal
    void *lmp = NULL;
    lmp = lammps_open(0, NULL, MPI_COMM_WORLD, NULL);

    lammps_close(lmp);
    MPI_Finalize();
  
    return 0;
}
