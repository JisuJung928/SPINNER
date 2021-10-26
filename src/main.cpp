#include <cstddef>
#include <iostream>
#include <mpi.h>
#include <unistd.h>

#include "calculator.h"
#include "config.h"
#include "input.h"

using namespace std;
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

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

    /* read input */
    Input *input = ReadInput("./INPUT");

    /* local MPI */
    int group_index = rank / input->GetNpar();
    int local_rank = rank % input->GetNpar();
    MPI_Comm lammps_comm;
    MPI_Comm_split(MPI_COMM_WORLD, group_index, rank, &lammps_comm);

    /* One-sided communication */
    int header;
    if (local_rank == 0) {
        header = 1;
    } else {
        header = 0;
    }
    MPI_Comm header_comm;
    MPI_Comm_split(MPI_COMM_WORLD, header, rank, &header_comm);

    MPI_Win win;
    int crystal_index = -1;
    MPI_Win_create(&crystal_index, (MPI_Aint)sizeof(int), sizeof(int),
                   MPI_INFO_NULL, header_comm, &win);

    /* generate crystal */
    vector<Crystal> crystal_vector;
    vector<double> energy_vector;
    for (int i = 0; i < input->GetGeneration(); ++i) {
        if (i == 0) {
            if (rank == 0) {
                crystal_vector = GenerateCrystal(input, input->GetPopulation());
            }
        } else {
            // TODO: evolution
            if (rank == 0) {
                crystal_vector = GenerateCrystal(input, input->GetPopulation());
            }
        }
        /* broadcast # of population from root */
        int population;
        if (rank == 0) {
            population = (int)crystal_vector.size();
        }
        MPI_Bcast(&population, 1, MPI_INT, 0, MPI_COMM_WORLD);

        /* broadcast all crystals from root to headers */
        unsigned int n_atoms;
        vector<atomStruct> as;
        latticeStruct ls;
        for (int j = 0; j < population; ++j) {
            if (rank == 0) {
                n_atoms = crystal_vector[j].numAtoms();
            }
            MPI_Bcast(&n_atoms, 1, MPI_UNSIGNED, 0, header_comm);

            if (local_rank == 0) {
                if (rank == 0) {
                    as = crystal_vector[j].getAtoms();
                } else {
                    as.resize(n_atoms);
                }
                MPI_Bcast(&as[0], n_atoms, mpi_atomStruct, 0, header_comm);
            }

            if (local_rank == 0) {
                if (rank == 0) {
                    ls = crystal_vector[j].getLattice();
                }
                MPI_Bcast(&ls, 1, mpi_latticeStruct, 0, header_comm);
            }

            if ((local_rank == 0) && (rank != 0)) {
                crystal_vector.push_back(Crystal(ls, as));
            }
        }

        int add = 1;
        Crystal crystal;
        while (crystal_index < population - 1) {
            /* epoch start */
            MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win);
            /* flush the value in rank 0 */
            MPI_Win_flush(0, win);
            /* add 1 to the value in rank 0 */
            MPI_Fetch_and_op(&add, &crystal_index, MPI_INT, 0, 0, MPI_SUM, win);
            /* get the value from rank 0 */
            MPI_Get(&crystal_index, 1, MPI_INT, 0, 0, 1, MPI_INT, win);
            /* epoch end */
            MPI_Win_unlock(0, win);

            /* Caution: local rank */
            MPI_Bcast(&crystal_index, 1, MPI_INT, 0, lammps_comm);

            /* broadcast the crystal from root of lammps_comm to others */
            if (local_rank == 0) {
                n_atoms = crystal_vector[crystal_index].numAtoms();
            }
            MPI_Bcast(&n_atoms, 1, MPI_UNSIGNED, 0, lammps_comm);

            if (local_rank == 0) {
                as = crystal_vector[crystal_index].getAtoms();
            } else {
                as.resize(n_atoms);
            }
            MPI_Bcast(&as[0], n_atoms, mpi_atomStruct, 0, lammps_comm);

            if (local_rank == 0) {
                ls = crystal_vector[crystal_index].getLattice();
            }
            MPI_Bcast(&ls, 1, mpi_latticeStruct, 0, lammps_comm);

            if (local_rank == 0) {
                crystal = crystal_vector[crystal_index];
            } else {
                crystal = Crystal(ls, as);
            }

            /* relax */
            double energy = Relax(input, &crystal, &lammps_comm, group_index);
            energy_vector.push_back(energy);

            /* update crystal_index */
            MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win);
            MPI_Win_flush(0, win);
            MPI_Get(&crystal_index, 1, MPI_INT, 0, 0, 1, MPI_INT, win);
            MPI_Win_unlock(0, win);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        /* initialize crystal_index */
        crystal_index = -1;

        /* Allgather */
    }

    delete input;

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Comm_free(&lammps_comm);
    MPI_Comm_free(&header_comm);
    MPI_Win_free(&win);
    MPI_Finalize();
}
