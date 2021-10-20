#include <cstddef>
#include <iostream>
#include <mpi.h>

#include "calculator.h"
#include "config.h"
#include "input.h"

using namespace std;
int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Input *input = ReadInput("./INPUT");

    vector<Crystal> crystal_vector;
    if (rank == 0) {
        crystal_vector = GenerateCrystal(input);
    }

    /* atomStruct type for MPI */
    MPI_Datatype as_type[4] =
    {
        MPI_UNSIGNED,
        MPI_DOUBLE,
        MPI_DOUBLE,
        MPI_DOUBLE
    };
    int as_block[4] = {1, 1, 1, 1};
    MPI_Aint as_disp[4];

    as_disp[0] = offsetof(atomStruct, atomicNum);
    as_disp[1] = offsetof(atomStruct, x);
    as_disp[2] = offsetof(atomStruct, y);
    as_disp[3] = offsetof(atomStruct, z);

    MPI_Datatype mpi_atomStruct;
    MPI_Type_create_struct(4, as_block, as_disp, as_type, &mpi_atomStruct);
    MPI_Type_commit(&mpi_atomStruct);

    unsigned int n_atoms;
    if (rank == 0) {
        n_atoms = crystal_vector[0].numAtoms();
    }
    MPI_Bcast(&n_atoms, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    vector<atomStruct> as;
    if (rank == 0) {
        as = crystal_vector[0].getAtoms();
    } else {
        as.resize(n_atoms);
    }
    MPI_Bcast(&as[0], n_atoms, mpi_atomStruct, 0, MPI_COMM_WORLD);

    cout << "rank: " << rank << "\n" << as[0].x << endl;
    cout << "rank: " << rank << "\n" << as[0].y << endl;
    cout << "rank: " << rank << "\n" << as[0].z << endl;

    /* latticeStruct type for MPI */
    MPI_Datatype ls_type[6] =
    {
        MPI_DOUBLE,
        MPI_DOUBLE,
        MPI_DOUBLE,
        MPI_DOUBLE,
        MPI_DOUBLE,
        MPI_DOUBLE
    };
    int ls_block[6] = {1, 1, 1, 1, 1, 1};
    MPI_Aint ls_disp[6];

    ls_disp[0] = offsetof(latticeStruct, a);
    ls_disp[1] = offsetof(latticeStruct, b);
    ls_disp[2] = offsetof(latticeStruct, c);
    ls_disp[3] = offsetof(latticeStruct, alpha);
    ls_disp[4] = offsetof(latticeStruct, beta);
    ls_disp[5] = offsetof(latticeStruct, gamma);

    MPI_Datatype mpi_latticeStruct;
    MPI_Type_create_struct(6, ls_block, ls_disp, ls_type, &mpi_latticeStruct);
    MPI_Type_commit(&mpi_latticeStruct);

    latticeStruct ls;
    if (rank == 0) {
        ls = crystal_vector[0].getLattice();
    }
    MPI_Bcast(&ls, 1, mpi_latticeStruct, 0, MPI_COMM_WORLD);

    cout << "rank: " << rank << "\n" << ls.a << endl;
    cout << "rank: " << rank << "\n" << ls.a << endl;
    cout << "rank: " << rank << "\n" << ls.a << endl;

    Crystal crystal = Crystal(ls, as);
    
    Relax(input, &crystal);

    delete input;

    return 0;
}
